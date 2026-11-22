"""Known passage graphs rasterized as walls. No trajectory is supplied."""

import unittest
from itertools import pairwise

import numpy as np
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from lidar_svx.topology import extract_network


def walls(solid):
    return np.argwhere(solid & ~ndimage.binary_erosion(solid))


def summary(xyz, edges):
    graph = coo_matrix((np.ones(len(edges)), edges.T), shape=(len(xyz), len(xyz)))
    components, labels = connected_components(graph, directed=False)
    degree = np.bincount(edges.ravel(), minlength=len(xyz))
    return components, len(edges) - len(xyz) + components, degree, labels


class TopologyTest(unittest.TestCase):
    def test_band_interface_retains_widest_contacts_and_resolved_loop(self):
        from lidar_svx.topology import _interface_connections

        # Two distinct arms split at 0 and rejoin at 5 four bands away.
        # Regions 5,6,7,8 create an additional K2,2 at just one interface.
        band = np.array([0, 1, 1, 2, 2, 3, 3, 4, 4])
        pairs = np.array(
            [
                (0, 1),
                (0, 2),
                (1, 3),
                (2, 4),
                (3, 5),
                (4, 5),
                (5, 7),
                (5, 8),
                (6, 7),
                (6, 8),
            ]
        )
        clearance = np.full(len(pairs), 2.0)
        clearance[6] = 1.0
        keep = _interface_connections(pairs, band, clearance)
        assert np.array_equal(np.flatnonzero(~keep), [6])
        count, rank, _, _ = summary(np.zeros((9, 3)), pairs[keep])
        self.assertEqual((count, rank), (1, 1))

    def test_ring_without_endpoints(self):
        x, y, z = np.indices((49, 49, 17)) - np.array([24, 24, 8])[:, None, None, None]
        tube = (np.hypot(x, y) - 17) ** 2 + z**2 <= 4**2
        xyz, edges, _ = extract_network(walls(tube), 1)
        components, rank, degree, _ = summary(xyz, edges)
        self.assertEqual((components, rank, set(degree)), (1, 1, {2}))
        self.assertTrue(tube[tuple(np.floor(xyz).astype(int).T)].all())

    def test_open_ended_tee(self):
        x, y, z = np.indices((65, 49, 15))
        main = (y - 10) ** 2 + (z - 7) ** 2 <= 5**2
        branch = ((x - 32) ** 2 + (z - 7) ** 2 <= 5**2) & (y >= 10)
        occupied = walls(main | branch)
        occupied = occupied[
            (occupied[:, 0] > 0) & (occupied[:, 0] < 64) & (occupied[:, 1] < 48)
        ]
        xyz, edges, _ = extract_network(occupied, 1, min_branch_length=4)
        components, rank, degree, _ = summary(xyz, edges)
        self.assertEqual((components, rank, sum(degree == 1)), (1, 0, 3))
        self.assertLess(xyz[:, 0].min(), 8)
        self.assertGreater(xyz[:, 0].max(), 57)
        self.assertGreater(
            xyz[:, 1].max(), 47 - 2 * 5
        )  # Open ends can retract by a tube diameter.

    def test_even_cross_section_preserves_long_branch(self):
        x, y, z = np.indices((65, 49, 19))
        main = (y - 10.5) ** 2 + (z - 8.5) ** 2 <= 4**2
        branch = ((x - 32.5) ** 2 + (z - 8.5) ** 2 <= 4**2) & (y >= 10)
        xyz, edges, _ = extract_network(walls(main | branch), 1, min_branch_length=4)
        components, rank, degree, _ = summary(xyz, edges)
        self.assertEqual((components, rank, sum(degree == 1)), (1, 0, 3))
        self.assertGreater(np.ptp(xyz[:, 0]), 50)
        self.assertGreater(np.ptp(xyz[:, 1]), 25)

    def test_even_cross_section_single_corridor(self):
        _, y, z = np.indices((57, 32, 32))
        tube = (y - 12.5) ** 2 + (z - 12.5) ** 2 <= 4**2
        xyz, edges, _ = extract_network(walls(tube), 1, min_branch_length=3)
        components, rank, degree, _ = summary(xyz, edges)
        self.assertEqual((components, rank, sum(degree == 1)), (1, 0, 2))
        self.assertGreater(np.ptp(xyz[:, 0]), 44)

    def test_crossing_on_two_levels(self):
        x, y, z = np.indices((49, 51, 29))
        lower = (y - 25) ** 2 + (z - 7) ** 2 <= 4**2
        upper = (x - 24) ** 2 + (z - 20) ** 2 <= 4**2
        xyz, edges, _ = extract_network(walls(lower | upper), 0.4, min_branch_length=2)
        components, rank, degree, labels = summary(xyz, edges)
        self.assertEqual((components, rank, sum(degree == 1)), (2, 0, 4))
        for component in range(components):
            self.assertLess(np.ptp(xyz[labels == component, 2]), 1.6)

    def test_open_ends_separated_by_missing_span(self):
        x, y, z = np.indices((69, 17, 17))
        shell = walls((y - 8) ** 2 + (z - 8) ** 2 <= 5**2)
        shell = shell[(shell[:, 0] > 0) & (shell[:, 0] < 68)]
        shell = shell[(shell[:, 0] < 27) | (shell[:, 0] > 39)]
        xyz, edges, _ = extract_network(shell, 1, min_branch_length=3)
        components, rank, degree, _ = summary(xyz, edges)
        self.assertEqual((components, rank, sum(degree == 1)), (2, 0, 4))
        self.assertFalse(np.any((xyz[:, 0] > 28) & (xyz[:, 0] < 39)))

    def test_rotated_ring_with_a_single_outlet(self):
        from scipy.spatial.transform import Rotation

        points = np.moveaxis(np.indices((81, 75, 57)), 0, -1)
        corners = np.array(
            [[-16, -14, 0], [16, -14, 0], [16, 14, 0], [-16, 14, 0], [-16, -14, 0]]
        )
        outlet = np.array([[16, 2, 0], [34, 2, 0]])
        rotation = Rotation.from_euler("zyx", [37, 18.5, 37], degrees=True)
        solid = np.zeros(points.shape[:3], dtype=bool)
        for curve in (corners, outlet):
            curve = rotation.apply(curve) + [40, 36, 26]
            for start, end in pairwise(curve):
                d = end - start
                t = np.clip(((points - start) @ d) / (d @ d), 0, 1)
                solid |= np.sum((points - start - t[..., None] * d) ** 2, axis=-1) <= 16
        xyz, edges, _ = extract_network(walls(solid), 1, min_branch_length=5)
        components, rank, degree, _ = summary(xyz, edges)
        self.assertEqual((components, rank, sum(degree == 1)), (1, 1, 1))
        self.assertTrue(solid[tuple(np.floor(xyz).astype(int).T)].all())

    def test_two_cycles_in_one_component(self):
        x, y, z = np.indices((89, 55, 17))
        left = (np.hypot(x - 25, y - 27) - 20) ** 2 + (z - 8) ** 2 <= 4**2
        right = (np.hypot(x - 63, y - 27) - 20) ** 2 + (z - 8) ** 2 <= 4**2
        xyz, edges, _ = extract_network(walls(left | right), 1, band_width=5)
        components, rank, degree, _ = summary(xyz, edges)
        self.assertEqual((components, rank, sum(degree == 1)), (1, 2, 0))

    def test_one_wall_gives_no_axis(self):
        coords = np.argwhere(np.ones((32, 19, 1), bool))
        with self.assertRaisesRegex(ValueError, "interior|network|axis"):
            extract_network(coords, 1)
