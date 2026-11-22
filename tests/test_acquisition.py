import laspy
import numpy as np
import pytest

from lidar_svx.acquisition import acquisition_report


def write_visits(path, *, break_time=False):
    rng = np.random.default_rng(121)
    radii = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    xyz: list[np.ndarray] = []
    times: list[np.float64] = []
    ids: list[int] = []
    y = np.r_[np.linspace(-8, 8, 45), np.linspace(8, -8, 45)]
    for block, centers in enumerate(y.reshape(-1, 5)):
        points = np.array([[np.cos(a), p, np.sin(a)] for p in centers for a in radii])
        stamps = block * 200 + np.arange(len(points), dtype=float)
        if break_time and block >= 4:
            stamps -= 300
        order = rng.permutation(len(points))
        xyz.extend(points[order])
        times.extend(stamps[order])
        ids.extend([block] * len(order))
    cloud = laspy.create(point_format=1, file_version="1.2")
    cloud.header.scales = [0.0001] * 3
    cloud.x, cloud.y, cloud.z = np.array(xyz).T
    cloud.point_source_id = ids
    cloud.gps_time = times
    cloud.write(path)


def test_return_cloud_crosses_section_twice_without_record_order_assumption(tmp_path):
    path = tmp_path / "there-and-back.laz"
    write_visits(path)
    report = acquisition_report(path, [0, 0.1, 0], axis=1, radius=2, slab_width=1)
    assert report["source"]["record_count"] == 2880
    assert report["window_count"] == 90
    assert [r["direction"] for r in report["crossings"]] == [1, -1]
    for visit in report["crossings"]:
        assert np.linalg.norm(np.array(visit["xyz"]) - [0, 0.1, 0]) < 0.02
        assert min(visit["slab_sector_returns"]) > 0
        assert visit["nearest_record_0based"] < 2880
    assert report["source"]["backward_record_steps_within_blocks"] > 100
    assert "trajectory" in report["interpretation"]


def test_blocks_with_overlapping_timestamps_cannot_validate_continuity(tmp_path):
    path = tmp_path / "overlap.laz"
    write_visits(path, break_time=True)
    with pytest.raises(ValueError, match="overlap|order"):
        acquisition_report(path, [0, 0, 0], axis=1, radius=2, slab_width=1)


def test_touching_section_at_turnaround_does_not_count_as_two_crossings(tmp_path):
    path = tmp_path / "turnaround.las"
    write_visits(path)
    # Width-five median at the doubled turnaround takes the preceding Y.
    report = acquisition_report(path, [0, 7.6364, 0], axis=1, radius=2, slab_width=1)
    assert report["crossings"] == []
