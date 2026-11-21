import csv
import hashlib
import json
import subprocess
import sys

import laspy
import numpy as np
from scipy import ndimage


def test_laz_to_survex_preserves_ring_and_source_witnesses(tmp_path):
    x, y, z = np.ogrid[-21:22, -21:22, -6:7]
    solid = (np.hypot(x, y) - 14) ** 2 + z**2 <= 4.5**2
    wall = solid ^ ndimage.binary_erosion(solid)
    xyz = np.argwhere(wall) + [419000.125, 620000.25, -70.5]
    las = laspy.create()
    las.header.offsets = [419000, 620000, -70]
    las.header.scales = [0.001] * 3
    las.x, las.y, las.z = xyz.T
    source = tmp_path / "ring.laz"
    las.write(source)
    output = tmp_path / "result"
    command = [
        sys.executable,
        "-m",
        "lidar_svx",
        str(source),
        str(output),
        "--voxel",
        "1",
    ]
    run = subprocess.run(command, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["graph"]["cycle_rank"] == 1
    assert manifest["graph"]["components"] == 1
    assert manifest["graph"]["length_m"] is None
    assert (
        manifest["source"]["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    )
    with (output / "witnesses.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == manifest["graph"]["stations"]
    for row in rows:
        np.testing.assert_allclose(
            [float(row[k]) for k in ("x", "y", "z")],
            xyz[int(row["record_0based"])],
            rtol=0,
            atol=1e-6,
        )
    assert (output / "network.svx").is_file()
    first = (output / "manifest.json").read_bytes()
    assert subprocess.run(command, capture_output=True).returncode != 0
    assert (output / "manifest.json").read_bytes() == first
