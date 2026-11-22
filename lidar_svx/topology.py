"""Extract a passage graph from connected bands of geodesic distance in 3D."""

import numpy as np
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra, minimum_spanning_tree

from .generalize import generalize
from .grid import SparseMask
from .volume import interior_components


def _face_edges(points, shape):
    """One undirected adjacency per pair of face-sharing interior cells."""
    padded = np.asarray(shape) + 2
    ids = np.ravel_multi_index((points + 1).T, padded)
    starts, ends = [], []
    for stride in (padded[1] * padded[2], padded[2], 1):
        targets = ids + stride
        neighbors = np.minimum(np.searchsorted(ids, targets), len(ids) - 1)
        matched = ids[neighbors] == targets
        starts.append(np.flatnonzero(matched))
        ends.append(neighbors[matched])
    return np.concatenate(starts), np.concatenate(ends)


def _interface_connections(pairs, band, clearance):
    """Retain a maximum-clearance forest of contacts per geodesic band boundary.

    The incoming and outgoing sides of a region receive separate IDs here.
    This keeps the forests disjoint: removing a contact requires an alternative
    that crosses only this boundary. Split-and-rejoin paths across multiple
    boundaries can therefore still form loops in the complete graph.
    """
    a, b = pairs.T
    forward = band[a] < band[b]
    lower = np.where(forward, a, b) + len(band)
    upper = np.where(forward, b, a)
    # Unique ranks favour wide contacts, with the input order as a stable tie-break.
    order = np.argsort(-clearance, kind="stable")
    rank = np.empty(len(pairs), dtype=np.intp)
    rank[order] = np.arange(1, len(pairs) + 1)
    size = 2 * len(band)
    interfaces = coo_matrix((rank, (lower, upper)), shape=(size, size))
    forest = minimum_spanning_tree(interfaces)
    keep = np.zeros(len(pairs), dtype=bool)
    keep[order[forest.data.astype(np.intp) - 1]] = True
    return keep


def _persistent_connections(pairs, region_count):
    """Suppress terminal branches spanning at most two distance bands."""
    adjacency: list[set[int]] = [set() for _ in range(region_count)]
    for a, b in pairs:
        adjacency[a].add(b)
        adjacency[b].add(a)
    while True:
        short = set()
        for tip, neighbors in enumerate(adjacency):
            if len(neighbors) != 1:
                continue
            chain = [tip]
            previous, node = tip, next(iter(neighbors))
            while len(chain) <= 2 and len(adjacency[node]) == 2:
                chain.append(node)
                previous, node = node, next(v for v in adjacency[node] if v != previous)
            if len(chain) <= 2 and len(adjacency[node]) >= 3:
                short.update(chain)
        if not short:
            break
        for u in sorted(short):
            for v in adjacency[u]:
                adjacency[v].remove(u)
            adjacency[u].clear()
    return np.array([b in adjacency[a] for a, b in pairs], dtype=bool)


def _route_regions(points, inside_graph, centers, ports, pairs):
    """Build one tree inside each band region, connected through its selected ports."""
    _, parent, _ = dijkstra(
        inside_graph,
        directed=False,
        indices=centers,
        min_only=True,
        return_predecessors=True,
    )
    legs = {tuple(sorted(pair)) for pair in ports}
    ends = set(ports.ravel())
    for node in sorted(ends):
        while parent[node] >= 0:
            edge = tuple(sorted((node, parent[node])))
            if edge in legs:
                break
            legs.add(edge)
            node = parent[node]

    # A non-terminal region's center serves as a routing seed only. Trim
    # back its single-ended stem so that it cannot become a false branch.
    degree = np.bincount(pairs.ravel(), minlength=len(centers))
    keep = ends | set(centers[degree == 1])
    neighbors: dict[int, set[int]] = {}
    for a, b in legs:
        neighbors.setdefault(a, set()).add(b)
        neighbors.setdefault(b, set()).add(a)
    dead = [u for u, adj in neighbors.items() if len(adj) == 1 and u not in keep]
    while dead:
        u = dead.pop()
        if len(neighbors[u]) != 1:
            continue
        v = neighbors[u].pop()
        neighbors[v].remove(u)
        legs.remove(tuple(sorted((u, v))))
        if v not in keep and len(neighbors[v]) == 1:
            dead.append(v)
    used, remapped = np.unique(sorted(legs), return_inverse=True)
    return points[used], remapped.reshape(-1, 2)


def _band_network(mask, width, *, boundary_diagnostics=None):
    points = mask.points if isinstance(mask, SparseMask) else np.argwhere(mask)
    empty = points[:0], np.empty((0, 2), dtype=int)
    if len(points) < 2:
        return empty
    a, b = _face_edges(points, mask.shape)
    shape = (len(points), len(points))
    grid = coo_matrix((np.ones(len(a), bool), (a, b)), shape=shape).tocsr()
    # Two sweeps give a reproducible distant seed and geodesic coordinate
    # without choosing entrances or declaring the farthest pair a survey path.
    first = dijkstra(grid, directed=False, indices=0, unweighted=True)
    distance = dijkstra(
        grid, directed=False, indices=int(np.argmax(first)), unweighted=True
    )
    band = np.floor(distance / width).astype(np.int64)
    same = band[a] == band[b]
    if same.all():
        return empty
    region_graph = coo_matrix(
        (np.ones(same.sum(), bool), (a[same], b[same])), shape=shape
    )
    nregions, labels = connected_components(region_graph, directed=False)
    # Distance is to the boundary of the inferred interior (includes occlusion).
    if isinstance(mask, SparseMask):
        clearance = mask.clearance()
    else:
        clearance = ndimage.distance_transform_edt(np.pad(mask, 1))[
            tuple((points + 1).T)
        ]
    middle = -np.abs((distance % width) - (width - 1) / 2)
    by_region = np.lexsort((middle, clearance, labels))
    centers = by_region[np.r_[np.diff(labels[by_region]) != 0, True]]

    links = np.column_stack((a[~same], b[~same]))
    pairs = np.sort(labels[links], axis=1)
    pair_key = pairs[:, 0].astype(np.int64) * nregions + pairs[:, 1]
    # For adjacent regions choose their widest face crossing; stable cell IDs
    # break ties. Distinct split and merge regions can form a cycle.
    order = np.lexsort((-links[:, 0], np.min(clearance[links], axis=1), pair_key))
    chosen = order[np.r_[np.diff(pair_key[order]) != 0, True]]
    if boundary_diagnostics is not None:
        keep = _interface_connections(
            pairs[chosen], band[centers], np.min(clearance[links[chosen]], axis=1)
        )
        boundary_diagnostics["single_boundary_cycles_removed"] += int((~keep).sum())
        chosen = chosen[keep]
    chosen = chosen[_persistent_connections(pairs[chosen], nregions)]
    if not len(chosen):
        return empty
    # Dijkstra favors short paths with clearance. Edges stay within a band.
    cost = 1 + 0.5 * width / np.minimum(clearance[a[same]], clearance[b[same]]) ** 2
    inside_graph = coo_matrix((cost, (a[same], b[same])), shape=shape).tocsr()
    return _route_regions(points, inside_graph, centers, links[chosen], pairs[chosen])


def _prune_branches(points, adj, threshold):
    """Peel terminal chains shorter than threshold, retaining cyclic cores."""
    while threshold > 0:
        discarded: set[tuple[int, int]] = set()
        for leaf, neighbors in enumerate(adj):
            if len(neighbors) != 1:
                continue
            prev, curr = leaf, next(iter(neighbors))
            chain = [(prev, curr)]
            length = float(np.linalg.norm(points[prev] - points[curr]))
            while length < threshold and len(adj[curr]) == 2:
                nxt = next(n for n in adj[curr] if n != prev)
                prev, curr = curr, nxt
                length += float(np.linalg.norm(points[prev] - points[curr]))
                chain.append((prev, curr))
            if length < threshold:
                discarded.update((min(a, b), max(a, b)) for a, b in chain)
        if not discarded:
            return
        for a, b in discarded:
            adj[a].remove(b)
            adj[b].remove(a)


def extract_network(
    occupied,
    voxel,
    *,
    min_branch_length=0.0,
    max_cells=48_000_000,
    band_width=None,
    tolerance=0.0,
    wall_span=None,
    global_parity=False,
    barrier_cells=None,
):
    """Return cell-center XYZ in source units, undirected station pairs and diagnostics.

    Add the physical origin of the input grid to the returned coordinates.
    Width, pruning length and generalization tolerance use source units.
    Low-count returns in barrier_cells veto weak links on the same voxel grid.
    """
    occupied = np.asarray(occupied)
    if occupied.ndim != 2 or occupied.shape[1] != 3 or not len(occupied):
        raise ValueError("Expected a non-empty N by 3 array of occupied cells")
    if not np.issubdtype(occupied.dtype, np.integer):
        raise ValueError("Expected integer voxel indices")
    if not np.isfinite(voxel) or voxel <= 0:
        raise ValueError("voxel must be finite and positive")
    if not np.isfinite(min_branch_length) or min_branch_length < 0:
        raise ValueError("min_branch_length must be finite and non-negative")
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("tolerance must be finite and non-negative")
    if not isinstance(max_cells, (int, np.integer)) or max_cells < 1:
        raise ValueError("max_cells must be a positive integer")
    if band_width is None:
        band_width = 5 * voxel
    if not np.isfinite(band_width) or band_width < voxel:
        raise ValueError("band_width must be finite and at least one voxel")
    low = occupied.min(axis=0)
    high = occupied.max(axis=0)
    shape = tuple(int(b) - int(a) + 1 for a, b in zip(low, high, strict=True))
    cells = occupied.astype(np.int64) - low.astype(np.int64)
    if global_parity and wall_span is not None:
        raise ValueError("choose either a local wall span or global parity")
    if not global_parity and wall_span is None:
        wall_span = 50 * voxel
    if wall_span is not None and (not np.isfinite(wall_span) or wall_span < voxel):
        raise ValueError("wall span must be finite and at least one voxel")
    local_span = None if wall_span is None else wall_span / voxel
    barriers = None
    if local_span is not None and barrier_cells is not None:
        barrier_cells = np.asarray(barrier_cells)
        if (
            barrier_cells.ndim != 2
            or barrier_cells.shape[1] != 3
            or not np.issubdtype(barrier_cells.dtype, np.integer)
        ):
            raise ValueError("Expected N by 3 integer voxel indices for barrier_cells")
        inside = ((barrier_cells >= low) & (barrier_cells <= high)).all(axis=1)
        barriers = barrier_cells[inside].astype(np.int64, copy=False) - low.astype(
            np.int64
        )
    components, info = interior_components(
        cells, shape, max_cells, local_span=local_span, barrier_cells=barriers
    )
    if not info["interior_cells"]:
        raise ValueError(
            "Insufficient enclosing wall evidence to infer a passage interior"
        )
    info.update(
        band_width=band_width,
        band_width_voxels=band_width / voxel,
        max_spur_bands=2,
        interior_components=0,
        interior_components_without_legs=0,
        edges_before_pruning=0,
        edges_before_generalization=0,
    )
    if local_span is not None:
        info["single_boundary_cycles_removed"] = 0
    all_points, all_edges = [], []
    station_count = 0
    for offset, mask in components:
        info["interior_components"] += 1
        points, edges = _band_network(
            mask,
            band_width / voxel,
            boundary_diagnostics=info if local_span is not None else None,
        )
        info["edges_before_pruning"] += len(edges)
        if len(edges) and min_branch_length > 0:
            adjacency: list[set[int]] = [set() for _ in points]
            for a, b in edges:
                adjacency[a].add(b)
                adjacency[b].add(a)
            _prune_branches(points, adjacency, min_branch_length / voxel)
            pairs = [
                (i, j)
                for i, neighbors in enumerate(adjacency)
                for j in sorted(neighbors)
                if i < j
            ]
            used, remapped = np.unique(
                np.asarray(pairs, dtype=int), return_inverse=True
            )
            points, edges = points[used], remapped.reshape(-1, 2)
        if not len(edges):
            info["interior_components_without_legs"] += 1
            continue
        info["edges_before_generalization"] += len(edges)
        if tolerance:
            points, edges = generalize(points, edges, mask, tolerance / voxel)
        all_points.append((points + offset + low + 0.5) * voxel)
        all_edges.append(edges + station_count)
        station_count += len(points)
    if not all_edges:
        raise ValueError(
            "No supported network edges survived band extraction and branch pruning"
        )
    return np.concatenate(all_points), np.concatenate(all_edges), info
