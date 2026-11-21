from pathlib import Path
import runpy

import numpy as np


_CHECKS = runpy.run_path(
    str(Path(__file__).resolve().parents[1] / "scripts" / "benchmark_continuity.py")
)


def test_reference_projection_preserves_fractional_legs_and_separate_components():
    points = np.array([[0.0, 0, 0], [0, 2, 0], [0, 4, 0], [2, 0, 0], [2, 2, 0]])
    edges = np.array([[0, 1], [1, 2], [3, 4]])
    result = _CHECKS["check_routes"](
        points,
        list("abcde"),
        edges,
        np.array([0, 0, 1]),
        [10, 11, 12, 13],
        np.array([[0.1, 0.25, 0], [0.1, 3, 0], [2.1, 1, 0], [0.1, 1, 0]]),
        [(10, 11), (10, 12), (10, 13)],
        2,
    )
    assert result[0]["connected"] and result[0]["shortest_graph_length"] == 2.75
    np.testing.assert_allclose(result[0]["crossings_at_y"], [[0, 2, 0]])
    assert result[1]["shortest_graph_length"] is None
    assert result[2]["shortest_graph_length"] == 0.75
    np.testing.assert_allclose(
        [a["distance"] for a in result[0]["anchors"]], [0.1, 0.1]
    )
