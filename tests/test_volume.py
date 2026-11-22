import runpy
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage

from lidar_svx.volume import DIRECTIONS, interior_components


def test_sparse_volume_matches_exhaustive_ray_count():
    rng = np.random.default_rng(653)
    wall = rng.random((12, 11, 14)) < 0.18
    score = np.zeros(wall.shape, dtype=int)
    cardinal = np.zeros_like(wall)
    for pos in np.argwhere(~wall):
        for direction in DIRECTIONS:
            odd = []
            for sign in (1, -1):
                runs = 0
                previous = False
                p = pos.copy()
                while True:
                    p += np.array(direction) * sign
                    if np.any(p < 0) or np.any(p >= wall.shape):
                        break
                    hit = wall[tuple(p)]
                    runs += int(hit and not previous)
                    previous = hit
                odd.append(runs % 2)
            if all(odd):
                score[tuple(pos)] += 1
                if sum(abs(v) for v in direction) == 1:
                    cardinal[tuple(pos)] = True
    strong = (score >= 6) & cardinal & ~wall
    majority = ndimage.convolve(
        strong.astype("u1"), np.ones((3, 3, 3), "u1"), mode="constant"
    )
    expected = (score >= 3) & cardinal & ~wall & (majority >= 14)
    boxes, info = interior_components(np.argwhere(wall), wall.shape, max_cells=10000)
    actual = np.zeros_like(wall)
    for origin, mask in boxes:
        p = np.argwhere(mask) + origin
        actual[tuple(p.T)] = True
    np.testing.assert_array_equal(actual, expected)
    assert info["interior_cells"] == actual.sum()


def test_large_empty_extent_does_not_allocate_its_bounding_box():
    # Two closed chambers with a huge diagonal gap.
    wall = np.argwhere(
        np.ones((9, 9, 9), dtype=bool) ^ np.pad(np.ones((7, 7, 7), bool), 1)
    )
    both = np.vstack((wall, wall + [3000, 2000, 1000]))
    boxes, _ = interior_components(both, tuple(both.max(axis=0) + 1), max_cells=20000)
    assert len(list(boxes)) == 2


def test_resource_limit_has_actionable_error():
    cells = np.array([[0, 0, 0], [100_000, 0, 0]])
    with pytest.raises(ValueError, match="limit"):
        interior_components(cells, (100_001, 1, 1), max_cells=2000)


def test_bounded_brackets_recover_passage_past_exterior_returns():
    x, y, z = np.indices((81, 55, 57))
    tube = (y - 27) ** 2 + (z - 27) ** 2 <= 5**2
    walls = tube & ~ndimage.binary_erosion(tube)
    # Two observed exterior surfaces add far-away runs along some rays.
    walls |= (x >= 20) & (x <= 61) & ((z == 43) | (y == 43))
    checkpoints = np.array([[8, 27, 27], [40, 27, 27], [70, 27, 27]])

    def centerline_covered(span):
        boxes, _ = interior_components(
            np.argwhere(walls), walls.shape, max_cells=500_000, local_span=span
        )
        for origin, mask in boxes:
            local = checkpoints - origin
            if (
                np.all((local >= 0) & (local < mask.shape))
                and mask[tuple(local.T)].all()
            ):
                return True
        return False

    assert not centerline_covered(None)
    assert centerline_covered(15)


def test_bounded_brackets_reject_rock_separating_parallel_tubes():
    _, y, z = np.indices((65, 41, 31))
    tubes = ((y - 12) ** 2 + (z - 15) ** 2 <= 4.5**2) | (
        (y - 23) ** 2 + (z - 15) ** 2 <= 4.5**2
    )
    occupied = np.argwhere(tubes & ~ndimage.binary_erosion(tubes))
    boxes, _ = interior_components(occupied, tubes.shape, 500_000, local_span=50)
    long_bodies = []
    for offset, mask in boxes:
        pts = offset + np.argwhere(mask)
        if np.ptp(pts[:, 0]) >= 20:
            long_bodies.append(pts)
    assert len(long_bodies) == 2
    for pts in long_bodies:
        assert tubes[tuple(pts.T)].all()


def test_supported_gap_requires_two_distinct_confident_cores():
    from lidar_svx.volume import _recover_gaps

    shape = (29, 5, 5)
    strides = np.array([25, 5, 1])

    def cells(start, end):
        x, y, z = np.indices((end - start, 2, 2))
        return np.ravel_multi_index(
            (x.ravel() + start, y.ravel() + 1, z.ravel() + 1), shape
        )

    core = np.r_[cells(2, 5), cells(9, 12)]
    # The other weak components either touch only one core or are isolated.
    weak = np.r_[cells(1, 2), cells(5, 9), cells(12, 14), cells(24, 27)]
    healed = _recover_gaps(core, weak, strides)
    np.testing.assert_array_equal(healed, cells(5, 9))


def test_low_count_returns_veto_weak_link_across_two_scanned_sidewalls():
    fixture = runpy.run_path(
        str(Path(__file__).resolve().parents[1] / "scripts" / "benchmark_continuity.py")
    )
    tubes, _, missing = fixture["cases"]()["parallel_passages"]
    xyz, _ = fixture["wall_returns"](tubes, missing, seed=40320)
    xyz = xyz.round(5)
    indices = np.floor((xyz - xyz.min(axis=0)) / 0.12).astype(int)
    shape = tuple(indices.max(axis=0) + 1)
    ids, counts = np.unique(np.ravel_multi_index(indices.T, shape), return_counts=True)
    occupied = np.column_stack(np.unravel_index(ids, shape))
    walls, low_count = occupied[counts >= 3], occupied[counts < 3]

    def long_bodies(barriers):
        boxes, _ = interior_components(
            walls, shape, 1_000_000, local_span=11 / 0.12, barrier_cells=barriers
        )
        return sum(max(mask.shape) > 100 for _, mask in boxes)

    assert long_bodies(None) == 1
    assert long_bodies(low_count) == 2
