"""Measure connected coverage on scans and on deterministic, independently generated tubes.

Synthetic truths are supplied only to the scorer. The extractor receives an XYZ
file in shuffled order, without timestamps, returns labelled by corridor, or guides.
"""

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation


@dataclass
class Tube:
    axis: np.ndarray
    radius: np.ndarray
    component: int = 0
    cap: int | None = None


def tube(axis, radius=1.25, component=0, cap=None):
    axis = np.asarray(axis)
    return Tube(axis, np.broadcast_to(radius, (len(axis),)), component, cap)


def line(start, end, component=0, cap=None):
    length = np.linalg.norm(np.subtract(end, start))
    return tube(
        np.linspace(start, end, int(np.ceil(length / 0.04)) + 1),
        component=component,
        cap=cap,
    )


def cases():
    x = np.linspace(0, 38, 1001)
    axis = np.column_stack((x, 2.6 * np.sin(x / 6), 1.4 * np.cos(x / 7)))
    narrow = tube(axis, 1.4 - 0.67 * np.exp(-(((x - 19) / 3) ** 2)))
    phase = np.linspace(0, 2 * np.pi, 1201)
    ring = tube(
        np.column_stack((7 * np.cos(phase), 7 * np.sin(phase), 0.4 * np.sin(2 * phase)))
    )
    # The two lobes cross at the same open junction, with oblique tangents.
    # Separate swept tubes remove the observed wall samples inside that junction.
    eight = [
        tube(
            np.column_stack((11 * np.sin(t), 6.5 * np.sin(2 * t), 0.35 * np.sin(3 * t)))
        )
        for t in (phase[:601], phase[600:])
    ]
    return {
        **{
            f"bend_missing_{angle}_degrees": ([narrow], 0, angle)
            for angle in (0, 115, 150, 180)
        },
        "ring_with_outlet": ([ring, line((7, 0, 0), (19, 0, 0))], 1, 0),
        "stacked_crossing": (
            [line((-15, 0, 0), (15, 0, 0)), line((0, -15, 3.5), (0, 15, 3.5), 1)],
            0,
            0,
        ),
        "parallel_passages": (
            [line((-15, 0, 0), (15, 0, 0)), line((-15, 3.5, 0), (15, 3.5, 0), 1)],
            0,
            0,
        ),
        "unobserved_span": (
            [line((-19, 0, 0), (-4, 0, 0)), line((4, 0, 0), (19, 0, 0), 1)],
            0,
            0,
        ),
        "solid_plug": (
            [
                line((-16, 0, 0), (-0.5, 0, 0), cap=-1),
                line((0.5, 0, 0), (16, 0, 0), 1, cap=0),
            ],
            0,
            0,
        ),
        "bend_below_surface_scan": ([narrow], 0, 0),
        "oblique_ring_with_outlet": ([ring, line((7, 0, 0), (19, 0, 0))], 1, 0),
        "figure_eight": (eight, 2, 0),
    }


def wall_returns(tubes, missing_angle_degrees, seed, terrain=False, holdout=False):
    rng = np.random.default_rng(seed)
    angle = np.linspace(0, 2 * np.pi, 144, endpoint=False)
    xyz = []
    for i, t in enumerate(tubes):
        tangent = np.gradient(t.axis, axis=0)
        tangent /= np.linalg.norm(tangent, axis=1, keepdims=True)
        first = np.cross(tangent, [0, 0, 1])
        first /= np.linalg.norm(first, axis=1, keepdims=True)
        second = np.cross(tangent, first)
        cos, sin = np.cos(angle)[None, :, None], np.sin(angle)[None, :, None]
        surface = t.axis[:, None] + t.radius[:, None, None] * (
            first[:, None] * cos + second[:, None] * sin
        )
        # Vary the count of observations with the look direction.
        surface += rng.normal(0, 0.014, surface.shape)
        keep = rng.random(surface.shape[:2]) < (0.76 + 0.20 * np.cos(angle))
        if missing_angle_degrees:
            middle = (t.axis[:, 0] > 13) & (t.axis[:, 0] < 25)
            keep &= ~(
                middle[:, None]
                & (angle > 0.85)
                & (angle < 0.85 + np.deg2rad(missing_angle_degrees))
            )
            keep &= ~middle[:, None] | (rng.random(keep.shape) < 0.68)
        surface = surface[keep]
        for other in tubes[:i] + tubes[i + 1 :]:
            distances, ids = cKDTree(other.axis).query(surface)
            surface = surface[distances >= other.radius[ids] - 0.015]
        xyz.append(surface)
        if t.cap is not None:
            rows = 3200
            u = t.radius[t.cap] * np.sqrt(rng.random(rows))
            theta = rng.uniform(0, 2 * np.pi, rows)
            disk = (
                t.axis[t.cap]
                + np.outer(u * np.cos(theta), first[t.cap])
                + np.outer(u * np.sin(theta), second[t.cap])
            )
            xyz.append(disk + rng.normal(0, 0.014, disk.shape))
    if terrain:
        # A registered above-ground patch above the tunnel. It is observed
        # from a different viewpoint and supplies no walls for a new passage.
        x, y = np.meshgrid(np.arange(-21, 60, 0.08), np.arange(-23, 24, 0.08))
        surface = np.column_stack(
            (x.ravel(), y.ravel(), (16 + 0.7 * np.sin(x / 7)).ravel())
        )
        xyz.append(surface + rng.normal(0, 0.014, surface.shape))
    # Fixed non-cardinal orientation; apply it after constructing the junction.
    angles = [-37, 27, 49] if holdout else [23, -12, 17]
    transform = Rotation.from_euler("zyx", angles, degrees=True)
    points = transform.apply(np.concatenate(xyz))
    shift = np.array([0.071, 0.043, 0.089]) if holdout else np.zeros(3)
    if holdout:
        # XYZ uses the cloud minimum as its grid origin. An isolated reference
        # target fixes this origin independently of the tube's subvoxel shift.
        fixed_target = np.floor(points.min(axis=0)) - 2
        points = np.vstack((points + shift, fixed_target))
    points = np.round(points, 5)
    rng.shuffle(points)
    truth = [
        tube(transform.apply(t.axis) + shift, t.radius, t.component, t.cap)
        for t in tubes
    ]
    return points, truth


def read_graph(directory):
    with (directory / "stations.csv").open(newline="") as f:
        stations = list(csv.DictReader(f))
    with (directory / "legs.csv").open(newline="") as f:
        legs = list(csv.DictReader(f))
    points = np.array([[float(row[k]) for k in "xyz"] for row in stations])
    names = [row["station"] for row in stations]
    ids = {name: i for i, name in enumerate(names)}
    edges = np.array([(ids[r["from"]], ids[r["to"]]) for r in legs])
    adjacency = coo_matrix(
        (np.ones(len(edges), bool), edges.T), shape=(len(points),) * 2
    )
    count, component = connected_components(adjacency, directed=False)
    return names, points, edges, count, component[edges[:, 0]]


def component_distances(samples, segments, components):
    """Exact distances from each sample to the finite legs of each component."""
    vec = segments[:, 1] - segments[:, 0]
    norm_sq = np.einsum("ij,ij->i", vec, vec)
    distances = np.full((len(samples), components.max() + 1), np.inf)
    for i, p in enumerate(samples):
        rel = p - segments[:, 0]
        t = np.clip(np.einsum("ij,ij->i", rel, vec) / norm_sq, 0, 1)
        off = rel - vec * t[:, None]
        np.minimum.at(distances[i], components, np.linalg.norm(off, axis=1))
    return distances


def resample(path, step):
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(path, axis=0), axis=1))]
    positions = np.linspace(0, s[-1], int(np.ceil(s[-1] / step)) + 1)
    return np.column_stack([np.interp(positions, s, column) for column in path.T])


def true_tips(truth):
    tips = []
    for i, t in enumerate(truth):
        if np.allclose(t.axis[0], t.axis[-1]):
            continue
        same_cave = [
            other
            for j, other in enumerate(truth)
            if i != j and other.component == t.component
        ]
        for end in (t.axis[0], t.axis[-1]):
            if not any(
                (np.linalg.norm(other.axis - end, axis=1) <= other.radius).any()
                for other in same_cave
            ):
                tips.append(end)
    return tips


def score_synthetic(directory, truth, rank):
    _, points, edges, count, components = read_graph(directory)
    segments = points[edges]
    groups = sorted({t.component for t in truth})
    covered = []
    sample_count = []
    curve_coverage = []
    curve_distances = []
    for group in groups:
        curves = [resample(t.axis, 0.20) for t in truth if t.component == group]
        samples = np.concatenate(curves)
        d = component_distances(samples, segments, components)
        covered.append(np.count_nonzero(d <= 0.8, axis=0))
        sample_count.append(len(samples))
        best = int(covered[-1].argmax())
        for part in np.split(d[:, best], np.cumsum([len(c) for c in curves])[:-1]):
            curve_coverage.append(float((part <= 0.8).mean()))
            curve_distances.append(np.round(np.percentile(part, [50, 90]), 4).tolist())
    covered = np.array(covered)
    share = covered / np.array(sample_count)[:, None]
    tips = true_tips(truth)
    ends = points[np.bincount(edges.ravel(), minlength=len(points)) == 1]
    tip_distance = component_distances(tips, segments, components).min(axis=1)
    tip_to_end = cKDTree(ends).query(tips)[0] if len(ends) and tips else []
    # Test finite legs inside the reference swept tubes. Weight each midpoint
    # by its own interval length instead of counting long and short legs equally.
    lengths = np.linalg.norm(segments[:, 1] - segments[:, 0], axis=1)
    probes, weights = [], []
    for segment, length in zip(segments, lengths):
        n = max(1, int(np.ceil(length / 0.08)))
        fraction = (0.5 + np.arange(n)) / n
        probes.append(segment[0] + fraction[:, None] * (segment[1] - segment[0]))
        weights.append(np.full(n, length / n))
    probes, weights = np.concatenate(probes), np.concatenate(weights)
    inside = np.zeros(len(probes), bool)
    for t in truth:
        distance, nearest = cKDTree(t.axis).query(probes)
        supported = distance <= t.radius[nearest]
        if not np.allclose(t.axis[0], t.axis[-1]):
            for end, neighbor in ((0, 1), (len(t.axis) - 1, -2)):
                outside = (probes - t.axis[end]) @ (t.axis[neighbor] - t.axis[end]) < 0
                supported &= ~((nearest == end) & outside)
        inside |= supported
    return {
        "components": int(count),
        "expected_components": len(groups),
        "cycle_rank": int(len(edges) - len(points) + count),
        "expected_cycle_rank": rank,
        "endpoints": int(np.count_nonzero(np.bincount(edges.ravel()) == 1)),
        "expected_endpoints": len(tips),
        "true_tips_distance_to_leg": np.round(tip_distance, 4).tolist(),
        "true_tips_distance_to_generated_endpoint": np.round(tip_to_end, 4).tolist(),
        "reference_samples": sample_count,
        "covered_reference_samples": covered.max(axis=1).tolist(),
        "connected_coverage_fraction": np.round(share.max(axis=1), 4).tolist(),
        "curve_coverage_in_best_component": np.round(curve_coverage, 4).tolist(),
        "curve_distance_to_best_component_p50_p90": curve_distances,
        "false_merge_components": int(
            np.count_nonzero((share >= 0.20).sum(axis=0) >= 2)
        ),
        "sampled_length_inside_true_tubes_fraction": float(
            np.round(weights[inside].sum() / weights.sum(), 5)
        ),
    }


def summarize_cases(scores):
    disjoint = [name for name, s in scores.items() if s["expected_components"] > 1]
    failed_topology = [
        name
        for name, s in scores.items()
        if s["components"] != s["expected_components"]
        or s["cycle_rank"] != s["expected_cycle_rank"]
        or s["false_merge_components"]
    ]
    failed_coverage = [
        name
        for name, s in scores.items()
        if min(s["curve_coverage_in_best_component"]) < 0.80
    ]
    failed_containment = [
        name
        for name, s in scores.items()
        if (s["sampled_length_inside_true_tubes_fraction"] or 0) < 0.99
    ]
    failed_endpoints = [
        name for name, s in scores.items() if s["endpoints"] != s["expected_endpoints"]
    ]
    return {
        "cases": len(scores),
        "component_count_matches": sum(
            s["components"] == s["expected_components"] for s in scores.values()
        ),
        "cycle_rank_matches": sum(
            s["cycle_rank"] == s["expected_cycle_rank"] for s in scores.values()
        ),
        "disjoint_control_cases": len(disjoint),
        "disjoint_controls_with_false_merges": [
            name for name in disjoint if scores[name]["false_merge_components"]
        ],
        "topology_matches_and_no_false_merge": len(scores) - len(failed_topology),
        "failed_topology": failed_topology,
        "endpoint_count_matches": len(scores) - len(failed_endpoints),
        "failed_endpoint_count": failed_endpoints,
        "component_cycle_endpoint_counts_and_separation_match": (
            len(scores) - len(set(failed_topology) | set(failed_endpoints))
        ),
        "at_least_80pct_of_every_curve_within_0_8_of_its_best_component": len(scores)
        - len(failed_coverage),
        "failed_curve_coverage": failed_coverage,
        "at_least_99pct_of_exported_length_inside_true_tubes": len(scores)
        - len(failed_containment),
        "failed_length_containment": failed_containment,
    }


def synthetic(generator_root, options, base_seed, holdout, selected_cases=None):
    scores = {}
    code_sha = {
        p.name: digest(p) for p in sorted((generator_root / "lidar_svx").glob("*.py"))
    }
    for index, (name, (tubes, rank, partial)) in enumerate(cases().items()):
        if selected_cases and name not in selected_cases:
            continue
        # The complete and obstructed views are subsets of the same noisy surface.
        terrain = name == "bend_below_surface_scan"
        seed = base_seed if terrain else base_seed + max(0, index - 3)
        if name == "oblique_ring_with_outlet" and not holdout:
            seed = base_seed + 18988  # Primary draw repeats seed 59305.
        points, truth = wall_returns(
            tubes, partial, seed=seed, terrain=terrain, holdout=holdout
        )
        if name == "oblique_ring_with_outlet":
            rotation = Rotation.from_euler("xyz", [16, 37, -22], degrees=True)
            points = np.round(rotation.apply(points), 5)
            truth = [
                tube(rotation.apply(t.axis), t.radius, t.component, t.cap)
                for t in truth
            ]
        with tempfile.TemporaryDirectory(prefix="svx-continuity-") as temp:
            source, result = Path(temp) / "cloud.xyz", Path(temp) / "result"
            np.savetxt(source, points, fmt="%.5f")
            process = subprocess.run(
                [sys.executable, "-m", "lidar_svx", str(source), str(result), *options],
                cwd=generator_root,
                capture_output=True,
                text=True,
            )
            scores[name] = {
                "source_sha256": digest(source),
                "returns": len(points),
                "seed": seed,
            }
            if process.returncode == 0:
                manifest = json.loads((result / "manifest.json").read_text())
                if manifest["program_sha256"] != code_sha:
                    raise RuntimeError(
                        "generator changed during benchmark; rerun with one version"
                    )
                scores[name].update(score_synthetic(result, truth, rank))
            else:
                message = process.stderr.strip()
                if not any(
                    text in message
                    for text in (
                        "No supported network edges survived",
                        "Insufficient enclosing wall evidence",
                    )
                ):
                    raise RuntimeError(f"{name}: {message}")
                groups = sorted({t.component for t in truth})
                sample_count = [
                    sum(len(resample(t.axis, 0.20)) for t in truth if t.component == g)
                    for g in groups
                ]
                scores[name].update(
                    no_network=message,
                    components=0,
                    expected_components=len(groups),
                    cycle_rank=0,
                    expected_cycle_rank=rank,
                    endpoints=0,
                    expected_endpoints=len(true_tips(truth)),
                    reference_samples=sample_count,
                    covered_reference_samples=[0] * len(groups),
                    connected_coverage_fraction=[0.0] * len(groups),
                    curve_coverage_in_best_component=[0.0] * len(tubes),
                    curve_distance_to_best_component_p50_p90=[None] * len(tubes),
                    false_merge_components=0,
                    sampled_length_inside_true_tubes_fraction=None,
                )
    return {
        "method": "XYZ wall observations shuffled with fixed seeds. Reference samples at <=0.20 spacing. "
        "Connected coverage is the maximum count within 0.8 of one exported component per truth component. "
        "False merge counts components covering >=20% of each of multiple disjoint truths. "
        "Tube containment is length-weighted on intervals <=0.08, using distance to a densely sampled axis.",
        "seed": base_seed,
        "orientation_zyx_degrees": [-37, 27, 49] if holdout else [23, -12, 17],
        "shift_xyz": [0.071, 0.043, 0.089] if holdout else [0, 0, 0],
        "independent_grid_target": holdout,
        "generator_source_sha256": code_sha,
        "generator_args": options,
        "summary": summarize_cases(scores),
        "cases": scores,
    }


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def native_check(directory, names, points, edges, manifest):
    model = directory / "network.3d"
    if not model.is_file():
        return None
    dump = subprocess.check_output(["dump3d", "-l", model], text=True)
    native, legs = {}, []
    for line in dump.splitlines():
        record = line.split()
        if record[:1] == ["NODE"]:
            native[record[4][1:-1]] = tuple(map(float, record[1:4]))
        elif record[:1] == ["LEG"]:
            legs.append(
                tuple(
                    sorted(
                        (tuple(map(float, record[1:4])), tuple(map(float, record[4:7])))
                    )
                )
            )
    absolute = (points - manifest["graph"]["source_origin_xyz"]) * (
        manifest["parameters"]["metres_per_unit"] or 1
    )
    exact_names = set(native) == {"lidar." + name for name in names}
    coords = [native["lidar." + name] for name in names] if exact_names else []
    expected = (
        [tuple(sorted((coords[a], coords[b]))) for a, b in edges] if exact_names else []
    )
    return {
        "sha256": digest(model),
        "stations": len(native),
        "legs": len(legs),
        "names_match": exact_names,
        "max_axis_difference": float(np.max(np.abs(absolute - coords)))
        if exact_names
        else None,
        "all_native_leg_geometries_match": Counter(legs) == Counter(expected),
    }


def linked_reference(axis_csv):
    with axis_csv.open(newline="") as f:
        axis = list(csv.DictReader(f))
    selected = set()
    for i, row in enumerate(axis[:-1]):
        if int(row["edge_exported"]):
            selected.update((i, i + 1))
    indices = sorted(selected)
    truth = np.array([[float(axis[i][key]) for key in "xyz"] for i in indices])
    return np.asarray([int(axis[i]["station"]) for i in indices]), truth


def missing_station_ranges(stations, hits):
    missing = stations[~hits]
    if not len(missing):
        return []
    breaks = np.flatnonzero(np.diff(missing) != 1) + 1
    return [[int(part[0]), int(part[-1])] for part in np.split(missing, breaks)]


def check_routes(points, names, edges, components, ids, truth, pairs, cross_y):
    """Project reference stations onto finite legs; walk the exported graph between projections."""
    segments = points[edges]
    direction = segments[:, 1] - segments[:, 0]
    length = np.linalg.norm(direction, axis=1)
    graph = coo_matrix((length, edges.T), shape=(len(points),) * 2).tocsr()
    reference = dict(zip(ids, truth, strict=True))
    result = []
    for pair in pairs:
        anchors, projections, legs = [], [], []
        for station in pair:
            if station not in reference:
                raise ValueError(
                    f"reference station {station} has no exported reference leg"
                )
            rel = reference[station] - segments[:, 0]
            t = np.clip(np.einsum("ij,ij->i", rel, direction) / length**2, 0, 1)
            projections_all = segments[:, 0] + t[:, None] * direction
            distances = np.linalg.norm(reference[station] - projections_all, axis=1)
            i = int(np.argmin(distances))
            anchors.append(
                {
                    "reference_station": station,
                    "reference_xyz": reference[station].tolist(),
                    "leg": [names[node] for node in edges[i]],
                    "component": int(components[i]),
                    "projection_xyz": projections_all[i].tolist(),
                    "distance": float(distances[i]),
                }
            )
            projections.append(projections_all[i])
            legs.append(i)
        row = {
            "anchors": anchors,
            "connected": bool(components[legs[0]] == components[legs[1]]),
            "shortest_graph_length": None,
            "crossings_at_y": [],
        }
        if row["connected"]:
            p, q = projections
            origin, target = edges[legs]
            if legs[0] == legs[1]:
                path = np.array([p, q])
            else:
                d, pred = dijkstra(
                    graph, directed=False, indices=origin, return_predecessors=True
                )
                costs = (
                    np.linalg.norm(points[origin] - p, axis=1)[:, None]
                    + d[:, target]
                    + np.linalg.norm(points[target] - q, axis=1)
                )
                a, b = np.unravel_index(np.argmin(costs), costs.shape)
                node = target[b]
                path_nodes = [node]
                while node != origin[a]:
                    node = pred[a, node]
                    path_nodes.append(node)
                path = np.vstack((p, points[path_nodes[::-1]], q))
            row["shortest_graph_length"] = float(
                np.linalg.norm(np.diff(path, axis=0), axis=1).sum()
            )
            if cross_y is not None:
                for start, end in zip(path, path[1:]):
                    if min(start[1], end[1]) <= cross_y < max(start[1], end[1]):
                        crossing = start + (
                            (cross_y - start[1]) / (end[1] - start[1])
                        ) * (end - start)
                        row["crossings_at_y"].append(crossing.tolist())
        result.append(row)
    return result


def real(directory, axis_csv, radius, pairs, cross_y):
    names, points, edges, count, components = read_graph(directory)
    manifest = json.loads((directory / "manifest.json").read_text())
    segments = points[edges]
    lengths = np.linalg.norm(segments[:, 1] - segments[:, 0], axis=1)
    weight = np.bincount(components, weights=lengths)
    main = int(np.argmax(weight))
    station_ids, truth = linked_reference(axis_csv)
    distances = component_distances(truth, segments, components)
    hits = distances <= radius
    audit_path = directory / "clearance.json"
    audit = None
    if audit_path.exists():
        c = json.loads(audit_path.read_text())
        for field, path in (
            ("source_sha256", None),
            ("stations_sha256", "stations.csv"),
            ("legs_sha256", "legs.csv"),
        ):
            if c[field] != (
                digest(directory / path) if path else manifest["source"]["sha256"]
            ):
                raise ValueError(
                    "clearance audit does not match input graph and source"
                )
        audit = {
            k: c[k]
            for k in (
                "radius_las_units",
                "sampled_legs",
                "total_legs",
                "legs_with_a_return_in_radius",
                "source_points_checked",
            )
        }
    return {
        "source_sha256": manifest["source"]["sha256"],
        "reference_sha256": digest(axis_csv),
        "stations_sha256": digest(directory / "stations.csv"),
        "legs_sha256": digest(directory / "legs.csv"),
        "components": int(count),
        "cycle_rank": int(len(edges) - len(points) + count),
        "length": float(lengths.sum()),
        "longest_component": main,
        "longest_component_length": float(weight[main]),
        "longest_component_share": float(weight[main] / weight.sum()),
        "reference_stations": len(truth),
        "radius": radius,
        "longest_component_coverage": int(hits[:, main].sum()),
        "any_component_coverage": int(hits.any(axis=1).sum()),
        "missing_station_id_ranges_longest_component": missing_station_ranges(
            station_ids, hits[:, main]
        ),
        "missing_station_id_ranges_any_component": missing_station_ranges(
            station_ids, hits.any(axis=1)
        ),
        "median_source_axis_distance_to_longest_component": float(
            np.median(distances[:, main])
        ),
        "median_source_axis_distance_to_any_component": float(
            np.median(distances.min(axis=1))
        ),
        "route_projection_method": "Nearest finite 3D leg to each reference station, ties resolved by CSV order. "
        "Shortest path between projections, including fractional endpoint legs. "
        "Component continuity and distances to the probes are reported independently.",
        "cross_y": cross_y,
        "route_checks": check_routes(
            points, names, edges, components, station_ids, truth, pairs, cross_y
        ),
        "component_details": [
            {
                "id": i,
                "legs": int((components == i).sum()),
                "length": float(weight[i]),
                "reference_coverage": int(hits[:, i].sum()),
            }
            for i in range(count)
        ],
        "raw_clearance": audit,
        "native_readback": native_check(directory, names, points, edges, manifest),
    }


def comparison_figure(source, before, after, axis_csv, radius, output):
    import laspy
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Patch

    ref_id, ref = linked_reference(axis_csv)
    graph = []
    source_digest = digest(source)
    for directory in (before, after):
        manifest = json.loads((directory / "manifest.json").read_text())
        if manifest["source"]["sha256"] != source_digest:
            raise ValueError("source differs from network manifest")
        _, xyz, edges, count, labels = read_graph(directory)
        segments = xyz[edges]
        lengths = np.linalg.norm(segments[:, 1] - segments[:, 0], axis=1)
        largest = np.bincount(labels, weights=lengths).argmax()
        hits = component_distances(ref, segments, labels) <= radius
        graph.append(
            (segments, labels == largest, count, hits[:, largest], hits.any(axis=1))
        )
    all_xy = np.concatenate(
        [ref[:, :2], *(g[0][:, :, :2].reshape(-1, 2) for g in graph)]
    )
    low, high = all_xy.min(axis=0) - 8, all_xy.max(axis=0) + 8
    step = max(high - low) / 1000
    nx, ny = np.ceil((high - low) / step).astype(int)
    density = np.zeros((ny, nx), np.int64)
    with laspy.open(source) as stream:
        for chunk in stream.chunk_iterator(500000):
            xy = np.column_stack((chunk.x, chunk.y))
            ij = np.floor((xy - low) / step).astype(int)
            i, j = ij[((ij >= 0) & (ij < [nx, ny])).all(axis=1)].T
            density += np.bincount(j * nx + i, minlength=nx * ny).reshape(density.shape)
    colors = {"main": "#007269", "other": "#bc5b29", "far": "#cdd2d5"}
    fig = plt.figure(figsize=(11.5, 7.2), layout="constrained")
    fig.get_layout_engine().set(h_pad=0.10, hspace=0.18)
    grid = fig.add_gridspec(2, 2, height_ratios=[4.6, 1.1])
    strip = fig.add_subplot(grid[1, :])
    for i, (segments, main, count, hit_main, hit_any) in enumerate(graph):
        ax = fig.add_subplot(grid[0, i])
        ax.imshow(
            np.log1p(density),
            origin="lower",
            extent=(low[0], low[0] + nx * step, low[1], low[1] + ny * step),
            cmap=LinearSegmentedColormap.from_list(
                "density", ["white", "#c0c3c9", "#555f69"]
            ),
            interpolation="nearest",
            aspect="equal",
        )
        for mask, color in ((~main, colors["other"]), (main, colors["main"])):
            ax.add_collection(
                LineCollection(segments[mask, :, :2], colors=color, linewidths=1.35)
            )
        stage = ("Poprzedni wynik", "Nowy wynik")[i]
        ax.set(
            title=f"{stage}   •   {count} skł.   •   {hit_main.sum()}/{len(ref)} stacji w progu",
            xlabel="X [j. LAS]",
            ylabel="Y [j. LAS]",
        )
        ax.spines[["top", "right"]].set_visible(False)
        # Each one-unit-wide bar represents one historical station, including blank missing IDs.
        color = np.where(
            hit_main, colors["main"], np.where(hit_any, colors["other"], colors["far"])
        )
        strip.bar(
            ref_id,
            np.full(len(ref), 0.58),
            bottom=1 - i - 0.29,
            width=1.01,
            color=color,
            linewidth=0,
        )
    strip.set(
        yticks=[0, 1],
        yticklabels=["Nowy", "Poprzedni"],
        ylim=(-0.5, 1.5),
        xlim=(ref_id.min() - 1, ref_id.max() + 1),
        xlabel="Numer stacji wcześniejszej osi wzdłuż skanu",
        title=f"Odległość od każdej stacji do wygenerowanego odcinka 3D, próg {radius:g} j. LAS",
    )
    strip.spines[["left", "top", "right"]].set_visible(False)
    strip.tick_params(axis="y", length=0)
    strip.legend(
        handles=[
            Patch(color=colors["main"], label="najdłuższa spójna składowa"),
            Patch(color=colors["other"], label="pozostałe odcinki"),
            Patch(color=colors["far"], label="brak odcinka w progu"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.63),
        ncol=3,
        frameon=False,
    )
    fig.suptitle(
        f"{source.stem}   •   ciągłość sieci na tych samych odbiciach", fontsize=15
    )
    fig.savefig(output, dpi=145)
    plt.close(fig)


def main():
    source_hash = digest(Path(__file__))
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    syn = sub.add_parser("synthetic")
    syn.add_argument(
        "--generator-root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    syn.add_argument("--output", type=Path)
    syn.add_argument("--case", action="append", choices=tuple(cases()))
    syn.add_argument(
        "--seed", type=int, help="seed for density, noise and record shuffling"
    )
    syn.add_argument(
        "--holdout",
        action="store_true",
        help="independent rotation, noise, subvoxel shift and fixed grid target",
    )
    syn.add_argument("generator_args", nargs=argparse.REMAINDER)
    survey = sub.add_parser("real")
    survey.add_argument("network_directory", type=Path)
    survey.add_argument("axis_csv", type=Path)
    survey.add_argument("--radius", type=float, default=2.0)
    survey.add_argument("--route-pair", nargs=2, type=int, action="append", default=[])
    survey.add_argument("--cross-y", type=float)
    survey.add_argument("--output", type=Path)
    compare = sub.add_parser("figure")
    for key in ("source", "before", "after", "axis_csv", "output"):
        compare.add_argument(key, type=Path)
    compare.add_argument("--radius", type=float, default=2.0)
    args = parser.parse_args()
    if args.action == "real":
        report = real(
            args.network_directory,
            args.axis_csv,
            args.radius,
            args.route_pair,
            args.cross_y,
        )
    elif args.action == "figure":
        comparison_figure(
            args.source,
            args.before,
            args.after,
            args.axis_csv,
            args.radius,
            args.output,
        )
        return
    else:
        flags = args.generator_args
        flags = flags[1:] if flags[:1] == ["--"] else flags
        report = synthetic(
            args.generator_root,
            flags
            or [
                "--voxel",
                "0.22",
                "--min-returns",
                "5",
                "--band-width",
                "2.2",
                "--min-branch-length",
                "4",
                "--tolerance",
                "0.3",
            ],
            args.seed if args.seed is not None else (95431 if args.holdout else 40317),
            args.holdout,
            args.case,
        )
    if source_hash != digest(Path(__file__)):
        raise RuntimeError("benchmark changed during this run")
    report["benchmark_sha256"] = source_hash
    output = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
