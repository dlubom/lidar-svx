"""A connected volume can span a large box while filling a small part of it."""

import numpy as np
from scipy.spatial import cKDTree


class SparseMask:
    """Interior cells in XYZ lexicographic order, with membership and Euclidean clearance."""

    def __init__(self, points, shape):
        self.points = points
        self.shape = shape
        self.padded = np.asarray(shape, dtype=np.int64) + 2
        self.ids = np.ravel_multi_index((points + 1).T, self.padded)

    def _contains(self, ids):
        i = np.minimum(np.searchsorted(self.ids, ids), len(self.ids) - 1)
        return self.ids[i] == ids

    def __getitem__(self, xyz):
        # The caller checks chords between interior stations, hence every
        # queried coordinate is in the common bounding box.
        return self._contains(np.ravel_multi_index(np.asarray(xyz) + 1, self.padded))

    def clearance(self):
        # For an interior cell a nearest exterior cell must share a face
        # with some interior cell. A step from any other exterior cell
        # towards the query would give a strictly closer exterior cell.
        boundary = []
        for stride in (self.padded[1] * self.padded[2], self.padded[2], 1):
            for step in (-stride, stride):
                neighbor = self.ids + step
                boundary.append(neighbor[~self._contains(neighbor)])
        keys = np.unique(np.concatenate(boundary))
        outside = np.column_stack(np.unravel_index(keys, self.padded)) - 1
        distances, _ = cKDTree(outside).query(self.points)
        return distances
