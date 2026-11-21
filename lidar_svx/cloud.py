"""Streaming, axis-preserving voxelisation. Each occupied cell keeps a raw witness."""

from dataclasses import dataclass
from itertools import islice
from math import prod
from pathlib import Path

import laspy
import numpy as np


@dataclass
class Cloud:
    origin: np.ndarray
    indices: np.ndarray
    counts: np.ndarray
    records: np.ndarray
    witness_xyz: np.ndarray
    metadata: dict
    barrier_cells: np.ndarray | None = None


def _text_chunks(path, size):
    with path.open(encoding="utf-8-sig") as stream:
        rows = (line.partition("#")[0].strip() for line in stream)
        rows = (line for line in rows if line)
        while batch := list(islice(rows, size)):
            xyz = np.loadtxt(batch, ndmin=2)
            if xyz.shape[1] != 3 or not np.isfinite(xyz).all():
                raise ValueError(
                    "XYZ requires three finite columns X Y Z on every data row"
                )
            yield xyz


def _las_chunks(path, size):
    with laspy.open(path) as reader:
        for chunk in reader.chunk_iterator(size):
            yield np.column_stack((chunk.x, chunk.y, chunk.z))


def _source_info(path, size):
    if path.suffix.lower() in (".las", ".laz"):
        with laspy.open(path) as reader:
            h = reader.header
            crs = h.parse_crs()
            if crs and (
                crs.is_geographic or any(c.is_geographic for c in crs.sub_crs_list)
            ):
                raise ValueError(
                    "reproject geographic LAS coordinates into a linear Cartesian frame first"
                )
            axes = crs.axis_info if crs else []
            factors = [a.unit_conversion_factor for a in axes]
            if factors and not np.allclose(factors, factors[0], rtol=1e-10, atol=0):
                raise ValueError("convert LAS axes to a common linear unit first")
            metadata = {
                "format": f"LAS {h.version}, point format {h.point_format.id}",
                "las_scales": h.scales.tolist(),
                "las_offsets": h.offsets.tolist(),
                "crs_wkt": crs.to_wkt() if crs else None,
                "crs_axis_units": [a.unit_name for a in axes],
                "crs_unit_factor": factors[0] if factors else None,
                "header_record_count": h.point_count,
            }
            return h.mins.copy(), h.maxs.copy(), metadata, _las_chunks
    if path.suffix.lower() not in (".xyz", ".txt"):
        raise ValueError(
            "supported inputs: .las, .laz, .xyz or .txt (three XYZ columns)"
        )
    low = np.full(3, np.inf)
    high = -low
    for points in _text_chunks(path, size):
        low = np.minimum(low, points.min(axis=0))
        high = np.maximum(high, points.max(axis=0))
    return low, high, {"format": "XYZ", "crs_wkt": None}, _text_chunks


_CELL = np.dtype(
    [("key", "i8"), ("count", "i8"), ("record", "i8"), ("xyz", "f8", (3,))]
)


def _merge(older, newer):
    """Sorted sparse histogram. Preserve the earliest record across batches."""
    keys = np.union1d(older["key"], newer["key"])
    result = np.empty(len(keys), dtype=_CELL)
    new_index = np.searchsorted(keys, newer["key"])
    old_index = np.searchsorted(keys, older["key"])
    result[new_index] = newer
    result[old_index] = older
    result["count"].fill(0)
    result["count"][old_index] = older["count"]
    result["count"][new_index] += newer["count"]
    return result


def read_cloud(
    path: Path,
    voxel: float,
    min_returns: int,
    chunk_size: int = 500_000,
    *,
    keep_barriers: bool = True,
) -> Cloud:
    if not np.isfinite(voxel) or voxel <= 0 or min_returns < 1 or chunk_size < 1:
        raise ValueError("voxel size, return threshold and batch size must be positive")
    low, high, metadata, chunks = _source_info(path, chunk_size)
    if not np.isfinite([low, high]).all() or np.any(high < low):
        raise ValueError("empty input or invalid XYZ bounds")
    shape = tuple(int(np.floor((b - a) / voxel)) + 1 for a, b in zip(low, high))
    if prod(shape) > np.iinfo(np.int64).max:
        raise ValueError(
            "voxel grid exceeds 64-bit indexing; increase --voxel or crop the scan"
        )
    # Merge pairs of equally sized batches to bound memory without rescanning old cells each time.
    levels = []
    count = 0
    for points in chunks(path, chunk_size):
        indices = np.floor((points - low) / voxel).astype(np.int64)
        try:
            flat = np.ravel_multi_index(indices.T, shape)
        except ValueError as error:
            raise ValueError(
                "point outside the source bounds; check LAS header extents"
            ) from error
        ids, first, counts = np.unique(flat, return_index=True, return_counts=True)
        cells = np.empty(len(ids), dtype=_CELL)
        cells["key"], cells["count"] = ids, counts
        cells["record"], cells["xyz"] = count + first, points[first]
        count += len(points)
        level = 0
        while level < len(levels) and levels[level] is not None:
            cells = _merge(levels[level], cells)
            levels[level] = None
            level += 1
        if level == len(levels):
            levels.append(None)
        levels[level] = cells
    if not count or count != metadata.get("header_record_count", count):
        raise ValueError("empty scan or point count inconsistent with LAS header")
    cells = None
    for batch in reversed(levels):
        if batch is not None:
            cells = batch if cells is None else _merge(cells, batch)
    metadata.update(
        record_count=count,
        occupied_before_filter=len(cells),
        bounds_xyz=[low.tolist(), high.tolist()],
    )
    retained = cells["count"] >= min_returns
    barrier_cells = None
    if keep_barriers and min_returns > 1:
        barrier_cells = np.column_stack(
            np.unravel_index(cells["key"][~retained], shape)
        )
    cells = cells[retained]
    if not len(cells):
        raise ValueError("no occupied voxels meet --min-returns")
    indices = np.column_stack(np.unravel_index(cells["key"], shape))
    metadata["retained_voxels"] = len(cells)
    metadata["returns_in_retained_voxels"] = int(cells["count"].sum())
    return Cloud(
        low,
        indices,
        cells["count"],
        cells["record"],
        cells["xyz"],
        metadata,
        barrier_cells,
    )
