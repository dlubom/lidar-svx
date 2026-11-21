import csv
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

from lidar_svx.export import write_network


class NetworkExportTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.output = Path(temp.name)
        # Triangle with a spur, and a second corridor above the same XY.
        offset = np.array([450000000.0, 5430000000.0, 25000.125])
        self.points = offset + np.array(
            [
                [0, 0, 0],
                [1000, 0, 0],
                [1000, 1000, 0],
                [1000, 2000, 0],
                [0, 0, 10000],
                [1000, 0, 11000],
            ]
        )
        self.edges = [(0, 1), (1, 2), (2, 0), (2, 3), (4, 5)]
        self.scale = 0.001

    def export(self):
        return write_network(
            self.output,
            self.points,
            iter(self.edges),
            self.scale,
            source_name="test.las",
            source_sha256="a" * 64,
        )

    def test_exports_source_frame_and_component_topology(self):
        stats = self.export()
        self.assertEqual(stats["stations"], 6)
        self.assertEqual(stats["legs"], 5)
        self.assertEqual(stats["cycle_rank"], 1)
        self.assertEqual(stats["components"], 2)
        self.assertEqual(stats["junctions"], 1)
        self.assertEqual(stats["endpoints"], 3)
        self.assertAlmostEqual(stats["length_m"], 3 + 2 * np.sqrt(2))
        self.assertEqual(stats["source_origin_xyz"], self.points[0].tolist())

        with (self.output / "stations.csv").open() as f:
            rows = list(csv.DictReader(f))
        for i, row in enumerate(rows):
            self.assertEqual(row["station"], f"s{i}")
            self.assertEqual(row["degree"], str(sum(i in e for e in self.edges)))
            np.testing.assert_array_equal(
                [float(row[k]) for k in "xyz"], self.points[i]
            )
        self.assertNotEqual(rows[2]["component"], rows[4]["component"])
        with (self.output / "legs.csv").open() as f:
            leg_rows = list(csv.DictReader(f))
        overlay = (self.output / "network.poly").read_text().strip().split("\n\n")
        self.assertEqual(len(overlay), len(self.edges))
        for row, part in zip(leg_rows, overlay, strict=True):
            ends = [int(row[k][1:]) for k in ("from", "to")]
            xyz = [[float(v) for v in line.split()] for line in part.splitlines()]
            np.testing.assert_array_equal(self.points[ends], xyz)
            self.assertAlmostEqual(
                float(row["length_m"]),
                np.linalg.norm(np.diff(xyz, axis=0)) * self.scale,
            )
        svx = (self.output / "network.svx").read_text()
        self.assertEqual(svx.count("*fix "), 2)
        self.assertIn("test.las", svx)
        self.assertIn("a" * 64, svx)

    @unittest.skipUnless(
        shutil.which("cavern") and shutil.which("dump3d"), "requires Survex"
    )
    def test_native_3d_preserves_every_station_leg_and_disconnected_component(self):
        self.export()
        subprocess.run(
            ["cavern", "network.svx"], cwd=self.output, check=True, capture_output=True
        )
        readback = subprocess.check_output(
            ["dump3d", "-l", "network.3d"], cwd=self.output, text=True
        )
        expected = (self.points - self.points[0]) * self.scale
        nodes = {}
        legs = []
        for line in readback.splitlines():
            tokens = line.split()
            if tokens[:1] == ["LEG"]:
                start, end = (
                    tuple(map(float, tokens[1:4])),
                    tuple(map(float, tokens[4:7])),
                )
                legs.append(frozenset((start, end)))
            elif tokens[:1] == ["NODE"]:
                nodes[tokens[4].strip("[]")] = tuple(map(float, tokens[1:4]))
        self.assertEqual(
            nodes, {f"lidar.s{i}": tuple(xyz) for i, xyz in enumerate(expected)}
        )
        self.assertCountEqual(
            legs,
            [
                frozenset((tuple(expected[i]), tuple(expected[j])))
                for i, j in self.edges
            ],
        )

    def test_rejects_invalid_or_ambiguous_graph_before_writing(self):
        for edges in ([], [(0, -1)], [(0, 6)], [(0, 0)], [(0, 1), (1, 0)], [(0.1, 1)]):
            with self.subTest(edges=edges), self.assertRaises(ValueError):
                write_network(self.output, self.points, edges, self.scale)
        self.assertEqual(list(self.output.iterdir()), [])

    def test_model_scale_is_explicit_in_standalone_outputs(self):
        stats = write_network(self.output, self.points, self.edges, 1, model_units=True)
        self.assertIsNone(stats["length_m"])
        self.assertAlmostEqual(
            stats["length_source_units"], 1000 * (3 + 2 * np.sqrt(2))
        )
        with (self.output / "legs.csv").open() as f:
            for row in csv.DictReader(f):
                self.assertEqual(row["length_m"], "")
                self.assertGreater(float(row["length_source_units"]), 0)
        self.assertIn("nominal", (self.output / "network.svx").read_text())
        title = next(
            line
            for line in (self.output / "network.svx").read_text().splitlines()
            if line.startswith("*title")
        )
        self.assertIn("nominal", title)

    def test_fundamental_cycles_list_only_existing_legs_and_include_closure(self):
        # Two independent cycles sharing an edge, plus the separate upper corridor.
        self.edges.append((3, 0))
        stats = self.export()
        with (self.output / "cycles.csv").open() as f:
            cycles = list(csv.DictReader(f))
        self.assertEqual(len(cycles), stats["cycle_rank"])
        self.assertEqual(len(cycles), 2)
        existing = {frozenset(e) for e in self.edges}
        used_walks = []
        for row in cycles:
            self.assertEqual(row["component"], "0")
            walk = [int(s[1:]) for s in row["stations"].split()]
            self.assertEqual(walk[0], walk[-1])
            self.assertEqual(len(walk) - 1, len(set(walk)))
            walked = {frozenset((a, b)) for a, b in zip(walk, walk[1:])}
            self.assertLessEqual(walked, existing)
            self.assertTrue(walked)
            used_walks.append(walked)
            length = np.linalg.norm(np.diff(self.points[walk], axis=0), axis=1).sum()
            self.assertAlmostEqual(length, float(row["perimeter_source_units"]))
        self.assertNotEqual(*used_walks)
        other = self.output / "other"
        write_network(
            other, self.points, ((b, a) for a, b in reversed(self.edges)), self.scale
        )
        self.assertEqual(
            (self.output / "cycles.csv").read_bytes(),
            (other / "cycles.csv").read_bytes(),
        )

    def test_tree_has_empty_cycle_table(self):
        write_network(self.output, self.points[:2], [(0, 1)], 1)
        with (self.output / "cycles.csv").open() as f:
            self.assertEqual(list(csv.DictReader(f)), [])

    def test_principal_survey_retains_station_names_and_other_component_in_full_export(
        self,
    ):
        stats = self.export()
        main = stats["principal_component"]
        svx = (self.output / main["survey"]).read_text()
        self.assertEqual(svx.count("*fix "), 1)
        selected = [
            line.split()[:2] for line in svx.splitlines() if line.startswith("s")
        ]
        with (self.output / "stations.csv").open() as stream:
            stations = {
                row["station"]: int(row["component"]) for row in csv.DictReader(stream)
            }
        with (self.output / "legs.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        expected = [
            [r["from"], r["to"]] for r in rows if stations[r["from"]] == main["id"]
        ]
        self.assertEqual(selected, expected)
        self.assertEqual(len(rows), len(self.edges))
        self.assertEqual(len(selected), main["legs"])
        self.assertGreater(main["length_fraction"], 0.5)

    def test_rejects_coordinates_beyond_native_3d_range_after_local_shift(self):
        with self.assertRaisesRegex(ValueError, "3d coordinate range"):
            write_network(self.output, [[0, 0, 0], [25000000, 0, 0]], [(0, 1)], 1)
        self.assertEqual(list(self.output.iterdir()), [])
