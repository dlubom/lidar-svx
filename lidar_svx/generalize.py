"""Reduce degree-two station chains, preserving graph topology and inferred interior."""

from collections import Counter
from itertools import product

import numpy as np


def _chains(edges, n):
    adjacency = [[] for _ in range(n)]
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    visited = set()
    chains = []
    # Walking from junctions first leaves only detached degree-two rings at the end.
    order = sorted(range(n), key=lambda i: len(adjacency[i]) == 2)
    for start in order:
        for end in sorted(adjacency[start]):
            if (start, end) in visited:
                continue
            chain = [start]
            a, b = start, end
            while True:
                chain.append(b)
                visited.update(((a, b), (b, a)))
                if len(adjacency[b]) != 2 or b == start:
                    break
                p, q = adjacency[b]
                a, b = b, p if q == a else q
            chains.append(np.asarray(chain))
    return chains


def _inside(start, delta, mask):
    """Check every traversed cell, including tangencies to a grid edge or corner."""
    crossings = [(np.arange(abs(d)) + 0.5) / abs(d) for d in delta if d]
    cuts = np.unique(np.concatenate(([0, 1], *crossings)))
    mid = (cuts[:-1] + cuts[1:]) / 2
    cells = np.floor(start + 0.5 + mid[:, None] * delta).astype(np.intp)
    if not mask[tuple(cells.T)].all():
        return False
    boundary = start + 0.5 + cuts[1:-1, None] * delta
    rounded = np.rint(boundary)
    tolerance = 16 * np.finfo(float).eps * max(mask.shape)
    on_plane = abs(boundary - rounded) <= tolerance
    for k in np.flatnonzero(on_plane.sum(axis=1) >= 2):
        options = [
            (int(r) - 1, int(r)) if tied else (int(np.floor(x)),)
            for x, r, tied in zip(boundary[k], rounded[k], on_plane[k])
        ]
        if not all(mask[index] for index in product(*options)):
            return False
    return True


def _reduce(ids, points, mask, tolerance):
    segments = []
    pending = [(0, len(ids) - 1)]
    while pending:
        lo, hi = pending.pop()
        a, b = ids[lo], ids[hi]
        if hi == lo + 1:
            segments.append((a, b))
            continue
        delta = points[b] - points[a]
        offsets = points[ids[lo + 1 : hi]] - points[a]
        t = np.clip((offsets @ delta) / float(delta @ delta), 0, 1)
        distances = np.linalg.norm(offsets - t[:, None] * delta, axis=1)
        worst = int(np.argmax(distances))
        if distances[worst] <= tolerance and _inside(points[a], delta, mask):
            segments.append((a, b))
        else:
            split = lo + 1 + worst
            pending.extend(((split, hi), (lo, split)))
    return segments


def generalize(points, edges, mask, tolerance_voxels):
    """Ramer–Douglas–Peucker in 3D with chord checks on the inferred volume."""
    chains = _chains(edges, len(points))
    pairs = Counter(tuple(sorted((c[0], c[-1]))) for c in chains)
    result = []
    for chain in chains:
        last = len(chain) - 1
        pair = tuple(sorted((chain[0], chain[-1])))
        cuts = [0, last]
        if pair[0] == pair[1]:
            cuts = [0, last // 3, (2 * last) // 3, last]
        elif pairs[pair] > 1 and last > 1:
            cuts = [0, last // 2, last]
        for lo, hi in zip(cuts, cuts[1:]):
            result.extend(_reduce(chain[lo : hi + 1], points, mask, tolerance_voxels))
    used, indices = np.unique(np.asarray(result), return_inverse=True)
    return points[used], indices.reshape(-1, 2)
