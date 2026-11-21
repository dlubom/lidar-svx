import hashlib
import json
import subprocess
import sys
from pathlib import Path

import laspy
import numpy as np


def compare_visits(tmp_path, points_a, points_b):
    cloud = laspy.create(point_format=1, file_version="1.2")
    cloud.header.scales = [0.0001] * 3
    cloud.x, cloud.y, cloud.z = np.vstack((points_a, points_b)).T
    source = tmp_path / "visits.laz"
    cloud.write(source)
    split = len(points_a)
    report = {
        "source": {"sha256": hashlib.sha256(source.read_bytes()).hexdigest()},
        "probe_xyz": [0, 0, 0],
        "plane_axis": "Y",
        "crossings": [
            {"source_record_range_half_open": [0, split]},
            {"source_record_range_half_open": [split, split + len(points_b)]},
        ],
    }
    section = tmp_path / "section.json"
    section.write_text(json.dumps(report))
    output = tmp_path / "fit.json"
    script = (
        Path(__file__).resolve().parents[1] / "scripts/compare_acquisition_visits.py"
    )
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            str(source),
            str(section),
            str(output),
            "--visits",
            "0",
            "1",
            "--extent",
            ".8",
            "--spacing",
            ".06",
        ],
        capture_output=True,
        text=True,
    )
    return result, output


def test_fit_recovers_known_offset_on_held_out_points_from_compound_surfaces(tmp_path):
    rng = np.random.default_rng(514)
    faces = []
    for i in range(3):
        face = rng.uniform(-1, 1, size=(5000, 3))
        face[:, i] = 0.5 + rng.normal(0, 0.001, len(face))
        faces.append(face)
    a = np.vstack(faces)
    shift = np.array([0.08, -0.06, 0.12])
    b = a + shift + rng.normal(0, 0.001, a.shape)
    result, output = compare_visits(tmp_path, a, b)
    assert result.returncode == 0, result.stderr
    r = json.loads(output.read_text())["forward"]
    np.testing.assert_allclose(r["fit_shift_xyz"], shift, atol=0.008)
    assert (
        r["holdout_distance_after_p50_p90"][0] < r["holdout_distance_before_p50_p90"][0]
    )
    assert r["source_holdout_cells"] > 150


def test_single_plane_does_not_determine_three_dimensional_offset(tmp_path):
    rng = np.random.default_rng(43)
    a = rng.uniform(-1, 1, size=(5000, 3))
    a[:, 2] = 0
    result, output = compare_visits(tmp_path, a, a + [0.1, 0, 0.15])
    assert result.returncode != 0 and "constrain" in result.stderr
    assert not output.exists()
