import numpy as np

from lidar_svx.generalize import generalize


def test_bend_cannot_be_cut_through_an_unknown_cell():
    points = np.array([(x, 0, 0) for x in range(6)] + [(5, y, 0) for y in range(1, 6)])
    mask = np.zeros((6, 6, 1), bool)
    mask[tuple(points.T)] = True
    edges = np.column_stack((np.arange(10), np.arange(1, 11)))
    xyz, legs = generalize(points, edges, mask, 20.0)
    assert len(xyz) == 3 and len(legs) == 2
    assert (5, 0, 0) in map(tuple, xyz)


def test_pure_cycle_and_two_routes_between_junctions_survive():
    # Diamond with two tails attached to opposite corners, plus isolated square.
    points = np.array(
        [
            (2, 2, 0),
            (3, 3, 0),
            (4, 2, 0),
            (3, 1, 0),
            (1, 2, 0),
            (5, 2, 0),
            (7, 2, 0),
            (8, 2, 0),
            (8, 3, 0),
            (7, 3, 0),
        ]
    )
    edges = np.array(
        [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (2, 5), (6, 7), (7, 8), (8, 9), (9, 6)]
    )
    xyz, legs = generalize(points, edges, np.ones((10, 5, 1), bool), 8.0)
    assert len(set(map(tuple, np.sort(legs, axis=1)))) == len(legs)
    assert all(a != b for a, b in legs)
    assert len(legs) - len(xyz) + 2 == 2
    assert sorted(np.bincount(legs.ravel())).count(3) == 2


def test_chord_grazing_unknown_corner_is_rejected():
    points = np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0)])
    mask = np.zeros((2, 2, 1), bool)
    mask[tuple(points.T)] = True
    xyz, legs = generalize(points, [(0, 1), (1, 2)], mask, 3.0)
    assert len(xyz) == 3 and len(legs) == 2
