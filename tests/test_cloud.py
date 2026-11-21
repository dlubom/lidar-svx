from pathlib import Path

import laspy
import numpy as np
import pytest

from lidar_svx.cloud import read_cloud


@pytest.mark.parametrize("suffix", [".xyz", ".las", ".laz"])
def test_all_chunks_contribute_and_witness_is_original_record(tmp_path: Path, suffix):
    xyz = np.array([[400100.25, -3, 9], [400102.5, -3, 11], [400100.5, -3, 9]])
    source = tmp_path / ("scan" + suffix)
    if suffix == ".xyz":
        np.savetxt(source, xyz, header="X Y Z")
    else:
        las = laspy.create(point_format=1, file_version="1.2")
        las.header.offsets = [400000, 0, 0]
        las.header.scales = [0.001] * 3
        las.x, las.y, las.z = xyz.T
        las.write(source)
    cloud = read_cloud(source, voxel=1.0, min_returns=2, chunk_size=1)
    assert cloud.metadata["record_count"] == 3
    assert cloud.metadata["occupied_before_filter"] == 2
    np.testing.assert_array_equal(cloud.indices, [[0, 0, 0]])
    np.testing.assert_array_equal(cloud.counts, [2])
    np.testing.assert_array_equal(cloud.barrier_cells, [[2, 0, 2]])
    np.testing.assert_array_equal(cloud.records, [0])
    np.testing.assert_allclose(cloud.witness_xyz, xyz[:1])
    np.testing.assert_allclose(cloud.origin, xyz.min(axis=0))


@pytest.mark.parametrize("content", ["", "1 2 3\nnan 2 4\n", "1 2\n"])
def test_invalid_text_stops_before_extraction(tmp_path, content):
    source = tmp_path / "broken.xyz"
    source.write_text(content)
    with pytest.raises(ValueError):
        read_cloud(source, voxel=1.0, min_returns=1)


def test_geographic_las_requires_reprojection(tmp_path):
    from pyproj import CRS

    las = laspy.create(point_format=3, file_version="1.2")
    las.header.add_crs(CRS.from_epsg(4326))
    las.x, las.y, las.z = [20.1], [50.1], [110]
    source = tmp_path / "degrees.las"
    las.write(source)
    with pytest.raises(ValueError, match="reproject"):
        read_cloud(source, voxel=1, min_returns=1)
