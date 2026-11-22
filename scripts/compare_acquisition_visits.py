"""Test local alignment of two disjoint visits from a section check; output a report only."""

import argparse
import csv
import hashlib
import json
from importlib.metadata import version
from pathlib import Path

import laspy
import numpy as np
from scipy.spatial import cKDTree


def _thin(points, origin, width):
    cells = np.floor((points - origin) / width).astype(np.int64)
    _, first = np.unique(cells, axis=0, return_index=True)
    return points[first]


def fit_translation(source, target, probe, *, extent, spacing):
    """Robust translation-only point-to-plane ICP with distinct evaluation samples."""
    if not np.isfinite([extent, spacing]).all() or min(extent, spacing) <= 0:
        raise ValueError("extent and spacing must be finite positive source distances")
    src = _thin(source[(abs(source - probe) <= extent).all(axis=1)], probe, spacing)
    # Split before target thinning and normal estimation. No raw target point
    # used in fitting may be used as its own validation match.
    target = target[(abs(target - probe) <= extent + 5 * spacing).all(axis=1)]
    fit_target = _thin(target[::2], probe, 0.625 * spacing)
    hold_target = target[1::2]
    patch = np.floor((src - probe) / (8 * spacing)).astype(np.int64).sum(axis=1) % 2
    train, hold = src[patch == 0], src[patch == 1]
    if min(len(fit_target), len(hold_target), len(train), len(hold)) < 32:
        raise ValueError("insufficient points to fit and independently check the shift")
    tree = cKDTree(fit_target)
    local = fit_target[tree.query(fit_target, k=16)[1]]
    centered = local - local.mean(axis=1, keepdims=True)
    eig, vectors = np.linalg.eigh(centered.swapaxes(1, 2) @ centered)
    normals = vectors[:, :, 0]
    planar = eig[:, 0] < 0.25 * eig[:, 1]
    shift = np.zeros(3)
    for iteration in range(30):  # noqa: B007 (reported after the loop)
        distance, j = tree.query(train + shift)
        n = normals[j]
        residual = np.einsum("ij,ij->i", n, train + shift - fit_target[j])
        use = planar[j] & (distance <= 5 * spacing)
        weight = use * np.minimum(1, 0.5 * spacing / np.maximum(abs(residual), 1e-15))
        lhs = (n.T * weight) @ n
        if use.sum() < 20 or np.linalg.cond(lhs) > 1e5:
            raise ValueError(
                "local correspondences do not constrain a translation in XYZ"
            )
        step = np.linalg.solve(lhs, (n.T * weight) @ (-residual))
        shift += step
        if np.linalg.norm(step) <= 0.005 * spacing:
            break
    else:
        raise ValueError("translation fit did not converge")
    held_tree = cKDTree(hold_target)
    return {
        "fit_shift_xyz": shift.tolist(),
        "iterations": iteration + 1,
        "source_fit_cells": len(train),
        "source_holdout_cells": len(hold),
        "target_fit_cells": len(fit_target),
        "target_holdout_returns": len(hold_target),
        "holdout_distance_before_p50_p90": np.quantile(
            held_tree.query(hold)[0], [0.5, 0.9]
        ).tolist(),
        "holdout_distance_after_p50_p90": np.quantile(
            held_tree.query(hold + shift)[0], [0.5, 0.9]
        ).tolist(),
    }


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _network_crossings(directory, report, axis):
    with (directory / "stations.csv").open() as f:
        stations = {r["station"]: r for r in csv.DictReader(f)}
    probe = np.asarray(report["probe_xyz"])
    if (
        json.loads((directory / "manifest.json").read_text())["source"]["sha256"]
        != report["source"]["sha256"]
    ):
        raise ValueError("survey and section have different source scans")
    crossings = []
    with (directory / "legs.csv").open() as f:
        for leg in csv.DictReader(f):
            ends = [stations[leg[k]] for k in ("from", "to")]
            a, b = (np.array([float(r[k]) for k in "xyz"]) for r in ends)
            if a[axis] == b[axis]:
                continue  # no unique intersection
            t = (probe[axis] - a[axis]) / (b[axis] - a[axis])
            pt = a + t * (b - a)
            if (
                0 <= t <= 1
                and np.linalg.norm(pt - probe) <= report["radius_source_units"]
            ):
                crossings.append(
                    {
                        "from": leg["from"],
                        "to": leg["to"],
                        "component": int(ends[0]["component"]),
                        "xyz": pt.tolist(),
                    }
                )
    return crossings


def _plot(path, report, visits, groups, intersections):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

    name = report["plane_axis"]
    axis = "XYZ".index(name)
    h, v = [i for i in range(3) if i != axis]
    probe = np.asarray(report["probe_xyz"])
    centers = np.array([visit["xyz"] for visit in visits])
    middle = centers.mean(axis=0)
    radius = report["radius_source_units"]
    fig, panels = plt.subplots(1, 2, figsize=(10.4, 5.8))
    fig.subplots_adjust(left=0.08, right=0.97, top=0.82, bottom=0.24, wspace=0.25)
    fig.suptitle(
        f"Przekrój przy {name} = {probe[axis]:g}   ·   {report['source']['name']}",
        x=0.075,
        ha="left",
        fontsize=15,
    )
    fig.text(
        0.075,
        0.895,
        f"Źródłowe XYZ. Warstwa ±{report['slab_width_source_units'] / 2:g} i promień wyboru {radius:g} w jednostkach LAS.",
        fontsize=10,
    )
    colors = ["#11658a", "#bb632d"]
    legend = []
    witnesses = []
    for group, visit, color in zip(groups, visits, colors, strict=True):
        centered = group - probe
        selected = np.flatnonzero(
            (abs(centered[:, axis]) <= report["slab_width_source_units"] / 2)
            & (np.linalg.norm(centered, axis=1) <= radius)
        )
        j = (
            selected[np.linalg.norm(group[selected] - visit["xyz"], axis=1).argmin()]
            if len(selected)
            else None
        )
        witnesses.append(
            None
            if j is None
            else {
                "record_0based": int(visit["source_record_range_half_open"][0] + j),
                "xyz": group[j].tolist(),
            }
        )
        ids = ", ".join(str(k) for k in sorted(set(visit["source_ids"])))
        legend.append(
            Line2D(
                [],
                [],
                linestyle="",
                marker=".",
                color=color,
                markersize=9,
                label=f"Odbicia ID {ids}",
            )
        )
        for ax in panels:
            ax.scatter(
                group[selected, h],
                group[selected, v],
                s=5,
                color=color,
                alpha=0.6,
                linewidths=0,
            )
            c = np.asarray(visit["xyz"])
            ax.scatter(
                c[h],
                c[v],
                s=95,
                marker="X",
                color=color,
                edgecolor="white",
                linewidths=1,
                zorder=4,
            )
    legend.append(
        Line2D(
            [],
            [],
            linestyle="",
            marker="X",
            color="#333",
            label="Środek rozkładu odbić",
        )
    )
    if intersections:
        legend.append(
            Line2D(
                [],
                [],
                linestyle="",
                marker="s",
                markerfacecolor="none",
                markeredgecolor="#181818",
                label="Przecięcie odcinka SVX",
            )
        )
    for ax, center, span, title in zip(
        panels,
        [probe, middle],
        [radius, radius / 3],
        ["Oba przejścia", "Zbliżenie przewężenia"],
        strict=True,
    ):
        for row in intersections:
            xyz = row["xyz"]
            ax.scatter(
                xyz[h],
                xyz[v],
                s=90,
                marker="s",
                edgecolor="#191919",
                facecolor="none",
                linewidth=1.5,
                zorder=5,
            )
        ax.set(
            xlabel=f"{'XYZ'[h]} [j. LAS]",
            ylabel=f"{'XYZ'[v]} [j. LAS]",
            xlim=center[h] + np.array([-span, span]),
            ylim=center[v] + np.array([-span, span]),
            aspect="equal",
        )
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.set_title(title, fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
    fig.legend(
        handles=legend,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.095),
        ncol=4,
        frameon=False,
        fontsize=9,
    )
    samples = " i ".join(str(w["record_0based"]) for w in witnesses if w)
    fig.text(
        0.075,
        0.045,
        f"Najbliższe rekordy w narysowanej warstwie mają numery {samples} (od zera).",
        fontsize=9,
    )
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return witnesses


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("sections", type=Path, help="JSON produced by lidar_svx.acquisition")
    p.add_argument("output", type=Path, help="new JSON comparison")
    p.add_argument("--visits", type=int, nargs=2, required=True, metavar=("I", "J"))
    p.add_argument(
        "--extent", type=float, required=True, help="cube half-side in source units"
    )
    p.add_argument(
        "--spacing", type=float, required=True, help="source thinning voxel side"
    )
    p.add_argument(
        "--plot", type=Path, help="new section PNG in original source coordinates"
    )
    p.add_argument(
        "--network",
        type=Path,
        help="optional network directory for the section overlay",
    )
    args = p.parse_args()
    if args.output.exists() or (args.plot and args.plot.exists()):
        p.error("output already exists")
    report = json.loads(args.sections.read_text())
    if digest(args.source) != report["source"]["sha256"]:
        p.error("section report and input scan have different SHA-256")
    if any(i < 0 or i >= len(report["crossings"]) for i in args.visits):
        p.error("invalid crossing index")
    ranges = [
        report["crossings"][i]["source_record_range_half_open"] for i in args.visits
    ]
    if max(r[0] for r in ranges) < min(r[1] for r in ranges):
        p.error("the selected visits share source records")
    with laspy.open(args.source) as reader:
        xyz = []
        for lo, hi in ranges:
            # A diagnostic bracket covers at most three source blocks of 250k.
            if not 0 <= lo < hi <= reader.header.point_count or hi - lo > 750_000:
                p.error("invalid source record bracket")
            reader.seek(lo)
            pts = reader.read_points(hi - lo)
            xyz.append(np.column_stack((pts.x, pts.y, pts.z)))
    first, second = xyz
    probe = np.array(report["probe_xyz"])
    try:
        result = {
            "source_sha256": report["source"]["sha256"],
            "section_report_sha256": digest(args.sections),
            "script_sha256": digest(Path(__file__)),
            "versions": {n: version(n) for n in ("numpy", "scipy", "laspy", "lazrs")},
            "probe_xyz": probe.tolist(),
            "source_ranges_half_open": ranges,
            "visit_indices_0based": args.visits,
            "cube_half_extent": args.extent,
            "spacing": args.spacing,
            "forward": fit_translation(
                first, second, probe, extent=args.extent, spacing=args.spacing
            ),
            "reverse": fit_translation(
                second, first, probe, extent=args.extent, spacing=args.spacing
            ),
            "method": "Translation-only point-to-plane ICP. Source cells split by alternating spatial patches of side 8*spacing; target split into disjoint alternating rows after clipping. Validation compares held-out source cells to raw held-out target returns. All vectors use original LAS axes and units.",
            "interpretation": "A reversible local offset that reduces held-out distances is consistent with misregistration. Differences in visible surfaces or moving objects can also contribute. Source coordinates and the survey graph are untouched.",
        }
        crossings = (
            _network_crossings(args.network, report, "XYZ".index(report["plane_axis"]))
            if args.network
            else []
        )
        if args.network:
            result["network_intersections"] = crossings
            result["network_legs_sha256"] = digest(args.network / "legs.csv")
            result["network_stations_sha256"] = digest(args.network / "stations.csv")
        if args.plot:
            visits = [report["crossings"][i] for i in args.visits]
            result["plotted_nearest_slab_witnesses"] = _plot(
                args.plot, report, visits, xyz, crossings
            )
            result["versions"]["matplotlib"] = version("matplotlib")
    except ValueError as error:
        p.exit(1, f"{error}\n")
    with args.output.open("x", encoding="utf-8") as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write("\n")


if __name__ == "__main__":
    main()
