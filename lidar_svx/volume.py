"""Estimate passage interiors with opposing wall returns on a sparse XYZ lattice."""

from itertools import product
from math import prod

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from .grid import SparseMask

# One of each pair of antipodal directions, cardinal axes first.
DIRECTIONS = sorted(
    (d for d in product((-1, 0, 1), repeat=3) if next((v for v in d if v), 0) == 1),
    key=lambda d: sum(v != 0 for v in d),
)


def _intervals(
    walls, direction, strides, local_span=None
) -> tuple[NDArray, NDArray, int]:
    """Flat-index gaps bounded by consecutive wall runs on the same ray."""
    d = np.asarray(direction)
    axis = np.flatnonzero(d)[0]
    a, b = [i for i in range(3) if i != axis]
    q = walls[:, axis]
    u, v = walls[:, a] - q * d[a], walls[:, b] - q * d[b]
    order = np.lexsort((q, v, u))
    q, u, v = q[order], u[order], v[order]
    new_ray = (u[1:] != u[:-1]) | (v[1:] != v[:-1])
    starts = np.r_[0, np.flatnonzero(new_ray | (np.diff(q) > 1)) + 1]
    ends = np.r_[starts[1:] - 1, len(q) - 1]
    boundary = new_ray[starts[1:] - 1]
    ray_first = np.r_[0, np.flatnonzero(boundary) + 1]
    run_count = np.diff(np.r_[ray_first, len(starts)])
    ordinal = np.arange(len(starts)) - np.repeat(ray_first, run_count)
    odd_left = (ordinal + 1) % 2 == 1
    odd_right = (np.repeat(run_count, run_count) - ordinal - 1) % 2 == 1
    if local_span is None:
        paired = np.flatnonzero(odd_left & odd_right)
    else:
        # A missing or extra wall run can spoil one end's parity test.
        # Requiring a bounded pair and an odd count from at least one end
        # recovers such intervals without filling the even/even rock gap
        # between two complete, closely spaced passage shells.
        paired = np.flatnonzero(~boundary & (odd_left | odd_right)[:-1])
        separation = q[starts[paired + 1]] - q[ends[paired]]
        paired = paired[separation <= local_span / np.linalg.norm(d)]
    last = ends[paired]
    lengths = q[starts[paired + 1]] - q[last] - 1
    step = int(d @ strides)
    return walls[order[last]] @ strides + step, lengths, step


def _span_batches(starts, lengths, step, max_cells):
    # Limit intermediate memory, including a malformed scan with two distant returns.
    if np.any(lengths > max_cells):
        raise ValueError(
            "one ray span exceeds --max-cells limit; crop or increase --voxel"
        )
    cumulative = np.r_[0, np.cumsum(lengths)]
    i = 0
    while i < len(lengths):
        hi = max(
            i + 1,
            np.searchsorted(cumulative, cumulative[i] + 250_000, side="right") - 1,
        )
        offsets = cumulative[i:hi] - cumulative[i]
        flat = np.repeat(starts[i:hi] - offsets * step, lengths[i:hi])
        flat += np.arange(len(flat)) * step
        yield flat
        i = hi


def _lookup(sorted_ids, ids):
    """Locate exact matches, including values outside the sorted table."""
    where = np.searchsorted(sorted_ids, ids)
    where = np.minimum(where, sorted_ids.size - 1)
    return where, sorted_ids[where] == ids


def _face_labels(ids, strides):
    """Label face-connected components of sorted padded voxel IDs."""
    starts, ends = [], []
    for step in strides:
        where, match = _lookup(ids, ids + step)
        starts.append(np.flatnonzero(match))
        ends.append(where[match])
    row, col = np.concatenate(starts), np.concatenate(ends)
    graph = coo_matrix(
        (np.ones(row.size, bool), (row, col)), shape=(len(ids), len(ids))
    )
    return connected_components(graph, directed=False)[1]


def _recover_gaps(core, residual, strides):
    """Keep a weaker region only where it borders two distinct confident cores."""
    if not len(core) or not len(residual):
        return residual[:0]
    core_labels = _face_labels(core, strides)
    sizes = np.bincount(core_labels)
    if np.count_nonzero(sizes >= 8) < 2:
        return residual[:0]
    regions = _face_labels(residual, strides)
    lowest = np.full(regions.max() + 1, len(sizes))
    highest = np.full(len(lowest), -1)
    for delta in (*strides, *-strides):
        ix, match = _lookup(core, residual + delta)
        adjacent = core_labels[ix]
        supported = match & (sizes[adjacent] >= 8)
        np.minimum.at(lowest, regions[supported], adjacent[supported])
        np.maximum.at(highest, regions[supported], adjacent[supported])
    return residual[(lowest < highest)[regions]]


def _component_boxes(ids, shape, strides, max_cells):
    """Keep each face-connected volume in its own XYZ box."""
    labels = _face_labels(ids, strides)
    order = np.argsort(labels, kind="stable")
    cuts = np.flatnonzero(np.diff(labels[order])) + 1
    for group in np.split(order, cuts):
        coords = np.column_stack(np.unravel_index(ids[group], shape))
        low = coords.min(axis=0)
        extent = tuple((coords.max(axis=0) - low + 1).tolist())
        if prod(extent) > min(1_000_000, max_cells):
            yield low - 1, SparseMask(coords - low, extent)
            continue
        mask = np.zeros(extent, bool)
        mask[tuple((coords - low).T)] = True
        yield low - 1, mask  # undo the exterior padding added by interior_components()


def interior_components(
    cells, shape, max_cells, *, local_span=None, barrier_cells=None
):
    """Return a lazy sequence of compact volume boxes and JSON-safe diagnostics."""
    if local_span is not None and (not np.isfinite(local_span) or local_span <= 0):
        raise ValueError("local wall span in voxels must be finite and positive")
    shape = tuple(int(n) + 2 for n in shape)
    if prod(shape) >= np.iinfo(np.int64).max:
        raise ValueError("voxel bounds exceed 64-bit addressing; increase --voxel")
    strides = np.array([shape[1] * shape[2], shape[2], 1])
    walls = np.asarray(cells, dtype=np.int64) + 1
    cardinal_batches = []
    for d in DIRECTIONS[:3]:
        starts, lengths, step = _intervals(walls, d, strides, local_span)
        if lengths.sum() > max_cells:
            raise ValueError(
                "cardinal ray cells exceed --max-cells limit; increase --voxel"
            )
        cardinal_batches.extend(_span_batches(starts, lengths, step, max_cells))
    info = {
        "wall_pair_method": "global_parity" if local_span is None else "local_brackets",
        "maximum_wall_span_voxels": local_span,
        "grid_shape": [n - 2 for n in shape],
        "occupied_cells": len(walls),
        "antipodal_directions": 13,
        "minimum_direction_votes": 6,
        "minimum_cardinal_votes": 1,
        "majority_neighbourhood": 27,
        "minimum_majority": 14,
        "minimum_votes_after_majority": 3,
        "interior_cells": 0,
    }
    if not cardinal_batches:
        return (), info
    # All admissible cells have at least one cardinal vote. Faraway outliers
    # therefore cost ray sorting work, but no dense bounding-box allocation.
    ids, votes = np.unique(np.concatenate(cardinal_batches), return_counts=True)
    del cardinal_batches
    if ids.size > max_cells:
        raise ValueError("candidate cells exceed --max-cells limit; increase --voxel")
    votes = votes.astype("u1")
    for d in DIRECTIONS[3:]:
        for flat in _span_batches(
            *_intervals(walls, d, strides, local_span), max_cells
        ):
            where, match = _lookup(ids, flat)
            votes[where[match]] += 1  # spans are disjoint within each direction
    where, match = _lookup(ids, walls @ strides)
    votes[where[match]] = 0
    if not np.any(votes >= 6):
        return (), info
    candidates = ids[votes >= 3]
    neighbors = np.zeros(len(candidates), dtype="u1")
    if local_span is not None:
        relaxed_neighbors = np.zeros_like(neighbors)
    for offset in product((-1, 0, 1), repeat=3):
        nearby, match = _lookup(ids, candidates + offset @ strides)
        scores = votes[nearby]
        neighbors += match & (scores >= 6)
        if local_span is not None:
            relaxed_neighbors += match & (scores >= 4)
    interior = candidates[neighbors >= 14]
    if local_span is not None:
        residual = candidates[(neighbors < 14) & (relaxed_neighbors >= 10)]
        before_veto = len(residual)
        if barrier_cells is not None and len(barrier_cells):
            # A low-count return cannot vote for a wall, but a weak link
            # may not pass through its observed voxel.
            barriers = np.unique((barrier_cells + 1) @ strides)
            residual = residual[~_lookup(barriers, residual)[1]]
        recovered = _recover_gaps(interior, residual, strides)
        info.update(
            core_cells=len(interior),
            relaxed_minimum_votes=4,
            relaxed_minimum_neighbors=10,
            minimum_anchor_component_cells=8,
            weak_cells_vetoed_by_returns=before_veto - len(residual),
            recovered_cells=len(recovered),
        )
        interior = np.union1d(interior, recovered)
    info["interior_cells"] = len(interior)
    return (
        _component_boxes(interior, shape, strides, max_cells) if len(interior) else ()
    ), info
