import numpy as np
from scipy import ndimage

from lidar_svx.generalize import generalize
from lidar_svx.grid import SparseMask
from lidar_svx.topology import _band_network


def test_sparse_corridor_matches_dense_clearance_routing_and_chords():
    x, y, z = np.indices((45, 19, 19))
    tube = (y - 9) ** 2 + (z - 9) ** 2 <= 16
    tube &= ~((x > 15) & (x < 20) & (y > 8) & (z < 10))
    pts = np.argwhere(tube)
    mask = SparseMask(pts, tube.shape)
    expected = ndimage.distance_transform_edt(np.pad(tube, 1))[tuple((pts + 1).T)]
    np.testing.assert_allclose(mask.clearance(), expected)
    masks = []
    for variant in (tube, mask):
        stations, legs = _band_network(variant, 4)
        masks.append(generalize(stations, legs, variant, 1.5))
    for a, b in zip(*masks, strict=True):
        np.testing.assert_array_equal(a, b)
