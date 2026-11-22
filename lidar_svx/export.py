"""Export an inferred undirected station graph in source and local coordinates."""

import csv
import json
import math
import operator
from collections import deque
from pathlib import Path

import numpy as np


def _spanning_forest(station_count, edges):
    adjacency: list[list[int]] = [[] for _ in range(station_count)]
    for start, end in edges:
        adjacency[start].append(end)
        adjacency[end].append(start)
    labels = [-1] * station_count
    parent = [-1] * station_count
    depth = [0] * station_count
    anchors: list[int] = []
    for station, neighbours in enumerate(adjacency):
        if not neighbours:
            raise ValueError(f"station s{station} has no leg")
        if labels[station] != -1:
            continue
        label = len(anchors)
        anchors.append(station)
        pending = deque([station])
        labels[station] = label
        while pending:
            current = pending.popleft()
            for neighbour in adjacency[current]:
                if labels[neighbour] == -1:
                    labels[neighbour] = label
                    parent[neighbour] = current
                    depth[neighbour] = depth[current] + 1
                    pending.append(neighbour)
    return labels, anchors, [len(neighbours) for neighbours in adjacency], parent, depth


def _fundamental_cycles(legs, parent, depth):
    """Close each non-tree leg with the unique forest path between its endpoints."""
    for a, b in legs:
        if parent[a] == b or parent[b] == a:
            continue
        left, right = [a], [b]
        while left[-1] != right[-1]:
            i, j = left[-1], right[-1]
            if depth[i] >= depth[j]:
                left.append(parent[i])
            if depth[j] >= depth[i]:
                right.append(parent[j])
        yield left + right[-2::-1] + [a]


def _write_svx(path, header, local, anchors, legs):
    with path.open("x", encoding="ascii") as stream:
        stream.write(header)
        for station in anchors:
            print(
                f"*fix s{station}", *(f"{v:.9f}" for v in local[station]), file=stream
            )
        for a, b in legs:
            print(f"s{a} s{b}", *(f"{v:.9f}" for v in local[b] - local[a]), file=stream)
        print("*end lidar", file=stream)


def write_network(
    directory: Path,
    points,
    edges,
    metres_per_unit: float,
    *,
    source_name: str = "",
    source_sha256: str = "",
    model_units: bool = False,
) -> dict:
    """Write into a fresh output directory. Station indices refer to original XYZ.

    One local fix per component preserves its location relative to the rest of
    the scan. The common origin keeps even large projected XYZ within the .3d
    coordinate range. No rotation or north reference is inferred.
    """
    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("expected finite source XYZ coordinates, shape (N, 3)")
    if not math.isfinite(metres_per_unit) or metres_per_unit <= 0:
        raise ValueError("metres_per_unit must be positive and finite")
    if model_units and metres_per_unit != 1:
        raise ValueError("model units require a nominal scale of 1")

    pairs = set()
    for x, y in edges:
        try:
            a, b = sorted((operator.index(x), operator.index(y)))
        except TypeError as error:
            raise ValueError("leg endpoints must be integer station indices") from error
        if not 0 <= a < b < len(points) or (a, b) in pairs:
            raise ValueError(f"invalid or repeated leg s{a} to s{b}")
        pairs.add((a, b))
    if not pairs:
        raise ValueError("the inferred graph contains no survey legs")
    legs = sorted(pairs)
    labels, anchors, degrees, parent, depth = _spanning_forest(len(points), legs)
    local = (points - points[0]) * metres_per_unit
    index = np.asarray(legs)
    lengths = np.linalg.norm(points[index[:, 1]] - points[index[:, 0]], axis=1)
    if (
        np.any(lengths == 0)
        or not np.isfinite(lengths).all()
        or not np.isfinite(local).all()
    ):
        raise ValueError(
            "graph contains coincident leg endpoints or unrepresentable coordinates"
        )
    # Survex .3d stores signed 32-bit centimetres.
    coordinate_limit = np.iinfo(np.int32)
    if (
        local.min() < coordinate_limit.min / 100
        or local.max() > coordinate_limit.max / 100
    ):
        raise ValueError(
            "local station exceeds the native Survex .3d coordinate range; check source scale"
        )

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    scale_note = (
        "; Unknown physical scale. One source unit is plotted as one nominal metre."
        if model_units
        else f"; Supplied metres per source unit {metres_per_unit}"
    )
    title = "LiDAR centreline hypothesis - local axes" + (
        ", nominal scale" if model_units else ""
    )
    header = "\n".join(
        (
            "; Automatic LiDAR centreline hypothesis in the local axes of the scan.",
            f"; Source file {json.dumps(str(source_name))} SHA-256 {json.dumps(source_sha256)}",
            "; Source XYZ origin " + " ".join(map(str, points[0].tolist())),
            scale_note,
            "; dx, dy, dz follow source +X, +Y, +Z. North and vertical require independent control.",
            "; Each local fix positions a model component. All legs are derived from station coordinates.",
            f'*title "{title}"',
            "*begin lidar",
            "*units dx dy dz metres",
            "*data cartesian from to dx dy dz",
            "",
        )
    )
    _write_svx(directory / "network.svx", header, local, anchors, legs)
    component_ids = np.asarray(labels)[index[:, 0]]
    component_lengths = np.bincount(component_ids, weights=lengths)
    principal = int(np.argmax(component_lengths))
    selected = index[component_ids == principal]
    principal_file = "network.svx"
    if len(anchors) > 1:
        principal_file = "principal.svx"
        _write_svx(
            directory / principal_file,
            "; Connected component with the greatest sum of leg lengths.\n" + header,
            local,
            [anchors[principal]],
            selected,
        )

    with (directory / "stations.csv").open("x", newline="", encoding="ascii") as stream:
        writer = csv.writer(stream)
        writer.writerow(("station", "component", "degree", "x", "y", "z"))
        for i, xyz in enumerate(points):
            writer.writerow((f"s{i}", labels[i], degrees[i], *xyz.tolist()))

    with (
        (directory / "legs.csv").open("x", newline="", encoding="ascii") as table,
        (directory / "network.poly").open("x", encoding="ascii") as poly,
    ):
        writer = csv.writer(table)
        writer.writerow(("from", "to", "length_source_units", "length_m"))
        for (a, b), length in zip(legs, lengths, strict=True):
            writer.writerow(
                (
                    f"s{a}",
                    f"s{b}",
                    float(length),
                    "" if model_units else float(length * metres_per_unit),
                )
            )
            # Each blank line starts a new CloudCompare polyline.
            print(*points[a].tolist(), file=poly)
            print(*points[b].tolist(), file=poly)
            print(file=poly)

    with (directory / "cycles.csv").open("x", newline="", encoding="ascii") as stream:
        writer = csv.writer(stream)
        writer.writerow(("cycle", "component", "stations", "perimeter_source_units"))
        for i, walk in enumerate(_fundamental_cycles(legs, parent, depth), start=1):
            perimeter = math.fsum(np.linalg.norm(np.diff(points[walk], axis=0), axis=1))
            writer.writerow(
                (i, labels[walk[0]], " ".join(f"s{s}" for s in walk), perimeter)
            )

    total = math.fsum(lengths)
    return {
        "principal_component": {
            "id": principal,
            "survey": principal_file,
            "stations": labels.count(principal),
            "legs": len(selected),
            "cycle_rank": len(selected) - labels.count(principal) + 1,
            "length_source_units": float(component_lengths[principal]),
            "length_fraction": float(component_lengths[principal] / total),
        },
        "stations": len(points),
        "legs": len(legs),
        "components": len(anchors),
        "cycle_rank": len(legs) - len(points) + len(anchors),
        "endpoints": degrees.count(1),
        "junctions": sum(degree >= 3 for degree in degrees),
        "source_origin_xyz": points[0].tolist(),
        "length_m": None if model_units else total * metres_per_unit,
        "length_source_units": total,
    }
