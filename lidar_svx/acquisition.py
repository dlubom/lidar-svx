"""Check a section against time-separated neighbourhoods of LAS returns.

This optional diagnostic applies to exporters with contiguous source-ID blocks
whose timestamp ranges are strictly ordered. Coordinates are return centroids
(robust spatial medians). They never enter the survey graph.
"""

import argparse
import hashlib
import json
from dataclasses import dataclass
from importlib.metadata import version
from itertools import pairwise
from pathlib import Path

import laspy
import numpy as np
from scipy.ndimage import median_filter


@dataclass
class Block:
    source_id: int
    first: int
    count: int


@dataclass
class Window:
    block: Block
    index: int
    first_tick: float
    last_tick: float
    center: np.ndarray
    split_distance: float


def _blocks(reader):
    """Stream one contiguous source block at a time; reject oversized blocks."""
    current = None
    samples: list[np.ndarray] = []
    stamps: list[np.ndarray] = []
    offset = 0
    for batch in reader.chunk_iterator(125_000):
        ids = np.asarray(batch.point_source_id)
        cuts = np.r_[0, np.flatnonzero(ids[:-1] != ids[1:]) + 1, len(ids)]
        for lo, hi in pairwise(cuts):
            source_id = int(ids[lo])
            if current is None or source_id != current.source_id:
                if current is not None:
                    yield current, np.concatenate(samples), np.concatenate(stamps)
                current = Block(source_id, int(offset + lo), 0)
                samples, stamps = [], []
            current.count += int(hi - lo)
            if current.count > 250_000:
                raise ValueError(
                    "source-ID block exceeds diagnostic limit of 250000 records"
                )
            points = batch[lo:hi]
            samples.append(np.column_stack((points.x, points.y, points.z)))
            stamps.append(np.asarray(points.gps_time))
        offset += len(ids)
    if current is not None:
        yield current, np.concatenate(samples), np.concatenate(stamps)


def _center(xyz):
    center = np.median(xyz, axis=0)
    # Scaled floor stabilises coincident returns in Weiszfeld's iteration.
    eps = max(0.001 * np.median(np.linalg.norm(xyz - center, axis=1)), 1e-12)
    for _ in range(5):
        weights = 1 / np.maximum(np.linalg.norm(xyz - center, axis=1), eps)
        center = np.average(xyz, weights=weights, axis=0)
    return center


def _read_windows(path):
    windows = []
    used = set()
    previous_max = -np.inf
    count = 0
    reversals = 0
    with laspy.open(path) as reader:
        if "gps_time" not in reader.header.point_format.dimension_names:
            raise ValueError("LAS has no per-return gps_time")
        for block, points, ticks in _blocks(reader):
            count += block.count
            if (
                block.source_id in used
                or not np.isfinite(ticks).all()
                or ticks.min() <= previous_max
            ):
                raise ValueError(
                    "repeated source ID or overlapping/unordered block timestamps"
                )
            if block.count < 40:
                raise ValueError(
                    "each source block needs at least 40 returns for five split checks"
                )
            used.add(block.source_id)
            previous_max = ticks.max()
            reversals += int(np.count_nonzero(np.diff(ticks) < 0))
            for i, selection in enumerate(
                np.array_split(np.argsort(ticks, kind="stable"), 5)
            ):
                xyz = points[selection]
                lo, hi = float(ticks[selection[0]]), float(ticks[selection[-1]])
                if lo == hi:
                    raise ValueError("constant timestamps in acquisition window")
                split = np.linalg.norm(_center(xyz[::2]) - _center(xyz[1::2]))
                windows.append(Window(block, i, lo, hi, _center(xyz), float(split)))
        if count != reader.header.point_count or len(used) < 2:
            raise ValueError(
                "need at least two ordered LAS source blocks and a complete scan"
            )
    return windows, len(used), count, reversals


def _raw_support(reader, first, end, center, axis, radius, slab_width):
    reader.seek(first)
    other = [i for i in range(3) if i != axis]
    sectors = np.zeros(8, dtype=int)
    best = np.inf
    for offset in range(first, end, 125_000):
        raw = reader.read_points(min(end - offset, 125_000))
        xyz = np.column_stack((raw.x, raw.y, raw.z))
        diff = xyz - center
        sq = np.einsum("ij,ij->i", diff, diff)
        i = int(sq.argmin())
        if sq[i] < best:
            best, record, witness = float(sq[i]), offset + i, xyz[i].tolist()
        near = (sq <= radius**2) & (np.abs(diff[:, axis]) <= slab_width / 2)
        angles = np.arctan2(diff[near, other[1]], diff[near, other[0]])
        sectors += np.histogram(angles, bins=np.linspace(-np.pi, np.pi, 9))[0]
    return {
        "source_record_range_half_open": [first, end],
        "nearest_record_0based": record,
        "nearest_xyz": witness,
        "nearest_record_distance": float(np.sqrt(best)),
        "slab_sector_returns": sectors.tolist(),
    }


def acquisition_report(path, probe, *, axis, radius, slab_width):
    probe = np.asarray(probe, dtype=float)
    if (
        probe.shape != (3,)
        or axis not in (0, 1, 2)
        or not np.isfinite([*probe, radius, slab_width]).all()
        or min(radius, slab_width) <= 0
    ):
        raise ValueError(
            "require finite XYZ probe, an axis and positive radius and slab width"
        )
    windows, blocks, count, reversals = _read_windows(path)
    ticks = np.array([[w.first_tick, w.last_tick] for w in windows])
    positions = np.array([w.center for w in windows])
    # Break before a timestamp discontinuity. Never smooth across a missing visit.
    durations = ticks[:, 1] - ticks[:, 0]
    breaks = (ticks[1:, 0] - ticks[:-1, 1]) > 2 * np.maximum(
        durations[1:], durations[:-1]
    )
    segments = np.split(np.arange(len(windows)), np.flatnonzero(breaks) + 1)
    crossings = []
    with laspy.open(path) as reader:
        for run in segments:
            centers = median_filter(positions[run], size=(5, 1), mode="nearest")
            delta = centers[:, axis] - probe[axis]
            epsilon = 16 * np.finfo(float).eps * max(1, abs(probe[axis]))
            off_plane = np.flatnonzero(abs(delta) > epsilon)
            for k, j in pairwise(off_plane):
                if (delta[k] < 0) == (delta[j] < 0) or j - k > 2:
                    continue  # Tangency, or an extended coplanar run without a unique crossing.
                start, end = centers[k], centers[j]
                if (
                    max(np.linalg.norm(start - probe), np.linalg.norm(end - probe))
                    > radius
                ):
                    continue
                t = -delta[k] / (delta[j] - delta[k])
                center = centers[k + 1] if j == k + 2 else start + t * (end - start)
                if np.linalg.norm(center - probe) > radius:
                    continue
                u, v = windows[run[k]], windows[run[j]]
                crossings.append(
                    {
                        "xyz": center.tolist(),
                        "direction": int(np.sign(end[axis] - start[axis])),
                        "window_indices_0based": [int(run[k]), int(run[j])],
                        "source_ids": [u.block.source_id, v.block.source_id],
                        "within_block_windows_0based": [u.index, v.index],
                        "bracket_proxy_xyz": [start.tolist(), end.tolist()],
                        "bracket_window_raw_gps_time_ranges": [
                            [u.first_tick, u.last_tick],
                            [v.first_tick, v.last_tick],
                        ],
                        **_raw_support(
                            reader,
                            u.block.first,
                            v.block.first + v.block.count,
                            center,
                            axis,
                            radius,
                            slab_width,
                        ),
                    }
                )
    with path.open("rb") as stream:
        sha = hashlib.file_digest(stream, "sha256").hexdigest()
    return {
        "source": {
            "name": path.name,
            "sha256": sha,
            "record_count": count,
            "ordered_blocks": blocks,
            "backward_record_steps_within_blocks": reversals,
        },
        "estimator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "window_count": len(windows),
        "time_discontinuities": int(breaks.sum()),
        "window_odd_even_disagreement_median_p95": np.quantile(
            [w.split_distance for w in windows], [0.5, 0.95]
        ).tolist(),
        "probe_xyz": probe.tolist(),
        "plane_axis": "XYZ"[axis],
        "radius_source_units": radius,
        "slab_width_source_units": slab_width,
        "crossings": crossings,
        "method": "Five time-ranked equal-count windows per LAS source block. Five Weiszfeld iterations per window. Coordinatewise rolling median of five centres. Section crossings require both bracketing centres within the probe radius. Raw return checks use those two windows' entire source blocks and eight equal sectors in the other two axes.",
        "interpretation": "Return-distribution proxies for visits to a section. A repeated crossing is indirect evidence; a measured scanner trajectory requires independent acquisition metadata. Split-half agreement measures sampling repeatability, not position accuracy. No survey legs are generated here.",
        "versions": {
            name: version(name) for name in ("numpy", "scipy", "laspy", "lazrs")
        },
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("output", type=Path, help="new JSON report")
    p.add_argument(
        "--probe", type=float, nargs=3, required=True, metavar=("X", "Y", "Z")
    )
    p.add_argument(
        "--axis",
        choices=list("XYZ"),
        required=True,
        help="section normal in source axes",
    )
    p.add_argument("--radius", type=float, required=True)
    p.add_argument("--slab-width", type=float, required=True)
    args = p.parse_args()
    if args.output.exists():
        p.error("output already exists")
    try:
        report = acquisition_report(
            args.source,
            args.probe,
            axis="XYZ".index(args.axis),
            radius=args.radius,
            slab_width=args.slab_width,
        )
        with args.output.open("x", encoding="utf-8") as f:
            json.dump(report, f, indent=2, allow_nan=False)
            f.write("\n")
    except (ValueError, OSError, laspy.errors.LaspyException) as error:
        p.exit(1, f"{error}\n")
    print(
        f"{len(report['crossings'])} section crossings from {report['source']['ordered_blocks']} ordered blocks"
    )


if __name__ == "__main__":
    main()
