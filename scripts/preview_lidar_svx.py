"""Plot a Survex network on LAS returns and check selected full 3D legs."""

import argparse
import csv
import hashlib
import json
import sys
from importlib.metadata import version
from pathlib import Path

import laspy
import matplotlib
import numpy as np
from scipy.spatial import cKDTree

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

PROJECTIONS = ((0, 1), (0, 2), (1, 2))


def digest(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def read_legs(directory):
    with (directory / "stations.csv").open(newline="") as source:
        stations = {
            row["station"]: tuple(float(row[axis]) for axis in "xyz")
            for row in csv.DictReader(source)
        }
    with (directory / "legs.csv").open(newline="") as source:
        names = [(row["from"], row["to"]) for row in csv.DictReader(source)]
    if not names:
        raise ValueError("The network contains no legs.")
    return names, np.array([[stations[a], stations[b]] for a, b in names])


def compare_reference(paths, source_digest, names, legs):
    axis_path, review_path = paths
    review = json.loads(review_path.read_text())
    if review["source_sha256"] != source_digest:
        raise ValueError("reference sections use a different source LAZ")
    with axis_path.open(newline="") as stream:
        axis = list(csv.DictReader(stream))
    endpoints: set[int] = set()
    for i, row in enumerate(axis[:-1]):
        if int(row["edge_exported"]):
            endpoints.update((i, i + 1))
    reviewed = set(review["reviewed_stations"])
    start = legs[:, 0]
    direction = legs[:, 1] - start
    length_sq = np.einsum("ij,ij->i", direction, direction)
    distances, reviewed_distances = [], []
    furthest = {}
    furthest_distance = -1.0
    for i in sorted(endpoints):
        row = axis[i]
        point = np.array([float(row[k]) for k in "xyz"])
        rel = point - start
        t = np.clip(np.einsum("ij,ij->i", rel, direction) / length_sq, 0, 1)
        offset = rel - t[:, None] * direction
        squared = np.einsum("ij,ij->i", offset, offset)
        j = int(np.argmin(squared))
        distance = float(np.sqrt(squared[j]))
        distances.append(distance)
        if int(row["station"]) in reviewed:
            reviewed_distances.append(distance)
        if distance > furthest_distance:
            furthest_distance = distance
            furthest = {
                "station": int(row["station"]),
                "xyz": point.tolist(),
                "distance": distance,
                "closest_generated_leg": names[j],
            }

    def summary(values):
        if not values:
            raise ValueError("reference contains no eligible stations")
        values = np.asarray(values)
        return {
            "count": len(values),
            "median": float(np.median(values)),
            "p90": float(np.quantile(values, 0.9)),
            "max": float(values.max()),
            "within": {
                str(t): int(np.count_nonzero(values <= t)) for t in (0.5, 1, 2, 5)
            },
        }

    return {
        "source_sha256": source_digest,
        "axis": {"path": str(axis_path), "sha256": digest(axis_path)},
        "review": {"path": str(review_path), "sha256": digest(review_path)},
        "units": "source LAS coordinate units",
        "method": "Use each unique endpoint of an edge_exported row (row to next row) "
        "once. For each station, take the closest point on any finite 3D output "
        "leg. The subset comprises the linked IDs listed in reviewed_stations. "
        "The earlier axis is a comparison hypothesis.",
        "all_linked_stations": summary(distances),
        "section_reviewed_linked_stations": summary(reviewed_distances),
        "furthest_reference_station": furthest,
    }


def check_segments(xyz, first_record, segments, radius, records):
    tree = cKDTree(xyz)
    a, b = segments[:, 0], segments[:, 1]
    directions = b - a
    lengths = np.linalg.norm(directions, axis=1)
    candidates = tree.query_ball_point((a + b) / 2, lengths / 2 + radius)
    for start, delta, length, ids, row in zip(
        a, directions, lengths, candidates, records, strict=True
    ):
        if not ids:
            continue
        rel = xyz[ids] - start
        t = np.clip(np.einsum("ij,j->i", rel, delta) / length**2, 0, 1)
        off = rel - t[:, None] * delta
        distance_sq = np.einsum("ij,ij->i", off, off)
        inside = distance_sq <= radius**2
        row["returns_within_radius"] += int(inside.sum())
        if not inside.any():
            continue
        minimum = distance_sq.min()
        nearest = np.sqrt(minimum)
        if row["distance"] is None or nearest < row["distance"]:
            point = int(np.min(np.asarray(ids)[distance_sq == minimum]))
            row.update(
                distance=float(nearest),
                source_record_0based=first_record + point,
                source_xyz=xyz[point].tolist(),
            )


def draw(path, source_name, legs, images, origin, pixel, radius, records):
    fig = plt.figure(figsize=(13.5, 7), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=(1.12, 1), hspace=0.14)
    axes = [
        fig.add_subplot(grid[:, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[1, 1]),
    ]
    cmap = LinearSegmentedColormap.from_list("returns", ["white", "#b6bfc1", "#425362"])
    contacts = [
        row["leg_row_0based"] for row in records if row["returns_within_radius"]
    ]
    for ax, (x, y), image in zip(axes, PROJECTIONS, images, strict=True):
        rows, cols = image.shape
        extent = (
            origin[x],
            origin[x] + cols * pixel,
            origin[y],
            origin[y] + rows * pixel,
        )
        ax.imshow(
            np.log1p(image),
            extent=extent,
            origin="lower",
            interpolation="nearest",
            cmap=cmap,
            aspect="equal",
        )
        ax.add_collection(
            LineCollection(legs[:, :, [x, y]], colors="#006b87", linewidths=0.95)
        )
        if contacts:
            ax.add_collection(
                LineCollection(
                    legs[contacts][:, :, [x, y]], colors="#d64312", linewidths=1.25
                )
            )
        ax.set(
            xlabel=f"{'XYZ'[x]} [j. LAS]",
            ylabel=f"{'XYZ'[y]} [j. LAS]",
            title=f"Rzut {'XYZ'[x]}{'XYZ'[y]}",
        )
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(
        handles=[
            Line2D(
                [],
                [],
                marker="s",
                ls="none",
                color="#768691",
                label="odbicia LiDAR (gęstość)",
            ),
            Line2D([], [], color="#006b87", label="wygenerowany ciąg"),
            Line2D(
                [],
                [],
                color="#d64312",
                label=f"odcinek z odbiciem w promieniu {radius:g} j. LAS",
            ),
        ],
        loc="lower left",
        framealpha=0.92,
        fontsize=9,
    )
    fig.suptitle(f"{source_name}   •   sieć Survex na chmurze źródłowej", fontsize=15)
    fig.supxlabel(
        f"Rzuty w osiach LAS. Piksel obrazu = {pixel:.3g} j. LAS. "
        f"Sprawdzono odległość odbić od {len(records)} z {len(legs)} odcinków w 3D.",
        fontsize=9,
    )
    fig.savefig(path, dpi=155)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("network_directory", type=Path)
    parser.add_argument(
        "--radius",
        required=True,
        type=float,
        help="3D distance threshold in LAS coordinate units",
    )
    parser.add_argument(
        "--check-legs",
        type=int,
        default=64,
        help="number of evenly spaced rows to test from legs.csv (default 64)",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        nargs=2,
        metavar=("AXIS_CSV", "REVIEW_JSON"),
        help="compare against the project's earlier section-refined axis "
        "(station,x,y,z,edge_exported CSV and its kontrola.json)",
    )
    args = parser.parse_args()
    if not np.isfinite(args.radius) or args.radius <= 0 or args.check_legs < 1:
        parser.error("radius and check-legs must be positive")
    directory = args.network_directory
    manifest = json.loads((directory / "manifest.json").read_text())
    source_digest = digest(args.source)
    if source_digest != manifest["source"]["sha256"]:
        parser.error("LAS SHA-256 differs from the network manifest")
    names, legs = read_legs(directory)
    if args.reference:
        reference = compare_reference(args.reference, source_digest, names, legs)
        reference.update(
            legs_sha256=digest(directory / "legs.csv"),
            stations_sha256=digest(directory / "stations.csv"),
        )
        (directory / "reference-check.json").write_text(
            json.dumps(reference, indent=2) + "\n"
        )
    sample = np.linspace(0, len(legs) - 1, min(len(legs), args.check_legs), dtype=int)
    check = legs[sample]
    if np.any(np.all(check[:, 0] == check[:, 1], axis=1)):
        parser.error("coincident stations in selected leg")
    rows = [
        dict(
            zip(("from", "to"), names[i], strict=True),
            leg_row_0based=int(i),
            returns_within_radius=0,
            distance=None,
            source_record_0based=None,
            source_xyz=None,
        )
        for i in sample
    ]
    first_record = 0
    with laspy.open(args.source) as reader:
        origin = np.minimum(reader.header.mins, legs.min(axis=(0, 1)))
        hi = np.maximum(reader.header.maxs, legs.max(axis=(0, 1)))
        pixel = max(hi - origin) / 1100
        shape = np.floor((hi - origin) / pixel).astype(int) + 1
        images = [
            np.zeros((shape[y], shape[x]), dtype=np.int64) for x, y in PROJECTIONS
        ]
        for chunk in reader.chunk_iterator(500_000):
            xyz = np.column_stack((chunk.x, chunk.y, chunk.z))
            coords = np.floor((xyz - origin) / pixel).astype(np.intp)
            for (x, y), image in zip(PROJECTIONS, images, strict=True):
                cell = coords[:, y] * image.shape[1] + coords[:, x]
                image += np.bincount(cell, minlength=image.size).reshape(image.shape)
            check_segments(xyz, first_record, check, args.radius, rows)
            first_record += len(xyz)
    draw(
        directory / "preview.png",
        args.source.name,
        legs,
        images,
        origin,
        pixel,
        args.radius,
        rows,
    )
    report = {
        "command": [str(Path(__file__).name), *sys.argv[1:]],
        "script_sha256": digest(Path(__file__)),
        "source_name": args.source.name,
        "source_sha256": source_digest,
        "stations_sha256": digest(directory / "stations.csv"),
        "legs_sha256": digest(directory / "legs.csv"),
        "source_points_checked": first_record,
        "total_legs": len(legs),
        "sampled_legs": len(rows),
        "radius_las_units": args.radius,
        "legs_with_a_return_in_radius": sum(
            row["returns_within_radius"] > 0 for row in rows
        ),
        "method": "Evenly spaced CSV row indices. All LAS records tested against the finite "
        "3D segments including endpoints. Null distance means no return within the radius.",
        "preview": {
            "projections": ["XY", "XZ", "YZ"],
            "pixel_las_units": float(pixel),
            "origin_xyz": origin.tolist(),
            "counts": "log1p, each panel independently scaled",
        },
        "versions": {
            "python": sys.version.split()[0],
            **{
                pkg: version(pkg)
                for pkg in ("numpy", "scipy", "laspy", "lazrs", "matplotlib")
            },
        },
        "legs": rows,
    }
    (directory / "clearance.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"{first_record} points; {report['legs_with_a_return_in_radius']}/{len(rows)} legs "
        f"with a return within {args.radius:g} LAS units"
    )


if __name__ == "__main__":
    main()
