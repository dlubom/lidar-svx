import argparse
import csv
import hashlib
import json
import platform
import sys
import tempfile
from importlib import metadata
from pathlib import Path

import numpy as np
from laspy.errors import LaspyException
from pyproj.exceptions import CRSError
from scipy.spatial import cKDTree

from .cloud import read_cloud
from .export import write_network
from .topology import extract_network


def _parser():
    parser = argparse.ArgumentParser(
        description="Build a candidate cave passage graph from registered Cartesian wall returns."
    )
    parser.add_argument(
        "source", type=Path, help="LAS, LAZ or whitespace XYZ in one linear unit"
    )
    parser.add_argument("output", type=Path, help="new result directory")
    parser.add_argument(
        "--voxel",
        type=float,
        required=True,
        metavar="SIZE",
        help="voxel side length in source units",
    )
    evidence = parser.add_mutually_exclusive_group()
    evidence.add_argument(
        "--wall-span",
        type=float,
        metavar="WIDTH",
        help="maximum separation of opposing local wall returns in source units (default 50 voxels)",
    )
    evidence.add_argument(
        "--global-parity",
        action="store_true",
        help="use the original two-sided global wall-run parity test",
    )
    parser.add_argument(
        "--min-returns",
        type=int,
        default=1,
        metavar="N",
        help="minimum number of input records per wall voxel (default 1)",
    )
    parser.add_argument(
        "--min-branch-length",
        type=float,
        default=0.0,
        metavar="L",
        help="terminal branch pruning threshold in source units (default 0)",
    )
    parser.add_argument(
        "--band-width",
        type=float,
        metavar="WIDTH",
        help="width of geodesic distance bands in source units (default 5 voxels)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.0,
        metavar="DIST",
        help="maximum deviation when simplifying station chains in source units (default 0)",
    )
    parser.add_argument(
        "--metres-per-unit",
        type=float,
        metavar="SCALE",
        help="verified metres per XYZ unit; omit for nominal display scale",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=48_000_000,
        metavar="N",
        help="limit on candidate cells for interior reconstruction (default 48000000)",
    )
    return parser


def _sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _witnesses(path, points, cloud, voxel):
    centres = cloud.origin + (cloud.indices + 0.5) * voxel
    _, cells = cKDTree(centres).query(points)
    with path.open("x", encoding="ascii", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "station",
                "record_0based",
                "x",
                "y",
                "z",
                "voxel_returns",
                "station_to_witness_source_units",
            )
        )
        for i, cell in enumerate(cells):
            witness = cloud.witness_xyz[cell]
            writer.writerow(
                (
                    f"s{i}",
                    cloud.records[cell],
                    *witness,
                    cloud.counts[cell],
                    np.linalg.norm(points[i] - witness),
                )
            )


def _run(args, argv):
    source = args.source.resolve(strict=True)
    if args.output.exists() or args.output.is_symlink():
        raise ValueError(
            f"output already exists: {args.output}; choose a new directory"
        )
    output = args.output.resolve()
    protected = Path(__file__).resolve().parents[1] / "data" / "raw"
    if output.is_relative_to(protected) or (
        source.parent.name == "raw" and output.is_relative_to(source.parent)
    ):
        raise ValueError("choose an output directory outside the raw data directory")
    if (
        not np.isfinite(args.min_branch_length)
        or args.min_branch_length < 0
        or args.max_cells < 1
    ):
        raise ValueError("invalid branch length or grid size limit")
    if not np.isfinite(args.tolerance) or args.tolerance < 0:
        raise ValueError(
            "station simplification tolerance must be finite and non-negative"
        )
    if args.band_width is not None and (
        not np.isfinite(args.band_width) or args.band_width <= 0
    ):
        raise ValueError("band width must be finite and positive")
    wall_span = args.wall_span if args.wall_span is not None else 50 * args.voxel
    if args.global_parity:
        wall_span = None
    elif not np.isfinite(wall_span) or wall_span < args.voxel:
        raise ValueError("wall span must be finite and at least one voxel")
    scale = args.metres_per_unit
    if scale is not None and (not np.isfinite(scale) or scale <= 0):
        raise ValueError("metres per unit must be positive and finite")

    initial = source.stat()
    cloud = read_cloud(
        source, args.voxel, args.min_returns, keep_barriers=not args.global_parity
    )
    crs_factor = cloud.metadata.get("crs_unit_factor")
    if crs_factor and scale and not np.isclose(scale, crs_factor, rtol=1e-7, atol=0):
        raise ValueError("supplied scale differs from the unit of the declared LAS CRS")
    source_digest = _sha256(source)
    final = source.stat()
    if (initial.st_ino, initial.st_size, initial.st_mtime_ns) != (
        final.st_ino,
        final.st_size,
        final.st_mtime_ns,
    ):
        raise ValueError("input changed during reading")

    local_points, edges, diagnostics = extract_network(
        cloud.indices,
        args.voxel,
        min_branch_length=args.min_branch_length,
        max_cells=args.max_cells,
        band_width=args.band_width,
        tolerance=args.tolerance,
        wall_span=wall_span,
        global_parity=args.global_parity,
        barrier_cells=cloud.barrier_cells,
    )
    if not len(edges):
        raise ValueError(
            "no supported survey legs at this resolution; inspect coverage and voxel size"
        )
    points = local_points + cloud.origin
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".lidar-svx-", dir=output.parent
    ) as staging:
        directory = Path(staging) / "result"
        graph = write_network(
            directory,
            points,
            edges,
            scale if scale is not None else 1.0,
            source_name=source.name,
            source_sha256=source_digest,
            model_units=scale is None,
        )
        _witnesses(directory / "witnesses.csv", points, cloud, args.voxel)
        package = Path(__file__).parent
        manifest = {
            "command": ["lidar-svx", *argv],
            "source": {"name": source.name, "sha256": source_digest, **cloud.metadata},
            "parameters": {
                "voxel": args.voxel,
                "wall_span": wall_span,
                "min_returns": args.min_returns,
                "tolerance": args.tolerance,
                "band_width": args.band_width
                if args.band_width is not None
                else 5 * args.voxel,
                "min_branch_length": args.min_branch_length,
                "max_cells": args.max_cells,
                "grid_origin_xyz": cloud.origin.tolist(),
                "metres_per_unit": scale,
            },
            "coordinates": {
                "orientation": "Input +X,+Y,+Z without rotation. No geographic north or gravity inferred.",
                "scale": "nominal display only"
                if scale is None
                else "supplied explicitly",
                "survex_local": "(source XYZ - graph.source_origin_xyz) * display scale",
                "cloudcompare": "Original source XYZ, no shift or scale applied in network.poly.",
            },
            "extraction": diagnostics,
            "graph": graph,
            "software": {
                "python": platform.python_version(),
                **{
                    name: metadata.version(name)
                    for name in ("numpy", "scipy", "laspy", "lazrs", "pyproj")
                },
            },
            "program_sha256": {
                file.name: _sha256(file) for file in sorted(package.glob("*.py"))
            },
            "interpretation": "Automatically inferred geometry. Endpoints, junctions and cycles require inspection against the scan.",
        }
        (directory / "manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        directory.rename(output)
    return graph


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        graph = _run(args, argv)
    except (
        ValueError,
        OSError,
        OverflowError,
        MemoryError,
        LaspyException,
        CRSError,
    ) as error:
        parser.exit(1, f"lidar-svx: {error}\n")
    print(
        f"{args.output / 'network.svx'}  ({graph['stations']} stations, {graph['legs']} legs, "
        f"{graph['components']} components, cycle rank {graph['cycle_rank']})"
    )
    return 0
