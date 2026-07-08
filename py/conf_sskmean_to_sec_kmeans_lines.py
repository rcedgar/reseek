#!/usr/bin/env python3
"""
Emit sec_kmeans tabbed lines for the fixed K-means centroids in src/conf.cpp
(SSKMEAN / InitLetter). Distances are in Angstroms; mean entries use sid_t
encoding: int(d*d*(100/16) + 0.5) per flat_dist_types.h::dist2sid.

Offset order matches conf.cpp Init() (ivalues / jvalues).
"""

# SSKMEAN rows from conf.cpp: (cluster, count ignored, x0..x8 Angstroms)
_MEANS_ANGSTROM = [
    [5.466, 5.218, 6.295, 5.472, 5.231, 5.478, 9.957, 5.326, 5.203],
    [6.561, 9.619, 12.5, 6.584, 9.54, 6.503, 17.84, 9.217, 9.588],
    [6.872, 10.26, 13.44, 6.863, 10.24, 6.839, 19.66, 10.12, 10.16],
    [6.003, 8.082, 10.41, 6.348, 8.912, 6.544, 15.41, 9.469, 8.833],
    [5.795, 8.276, 10.66, 6.402, 9.297, 6.554, 12.19, 9.211, 5.914],
    [6.581, 9.627, 12.72, 6.623, 9.836, 6.776, 15.52, 9.617, 7.687],
    [6.506, 9.369, 11.06, 6.583, 8.264, 5.669, 11.32, 5.711, 8.937],
    [5.573, 5.537, 6.667, 5.473, 5.418, 5.498, 11.02, 5.621, 8.287],
    [5.679, 7.569, 9.335, 6.127, 8.475, 5.949, 8.031, 7.051, 5.713],
    [6.423, 8.457, 7.479, 5.695, 5.693, 5.486, 10.18, 5.738, 9.094],
    [6.636, 9.453, 11.96, 6.413, 9.031, 5.945, 14.73, 6.715, 9.583],
    [5.474, 5.71, 7.441, 5.583, 6.551, 6.021, 12.0, 8.539, 5.595],
    [5.574, 5.368, 6.191, 5.682, 7.545, 6.115, 6.775, 8.909, 6.277],
    [6.387, 8.708, 8.938, 5.809, 6.408, 6.041, 13.26, 8.584, 9.058],
    [5.867, 6.271, 8.535, 5.848, 8.319, 6.471, 10.82, 9.066, 8.534],
    [6.515, 7.566, 6.319, 5.675, 5.571, 5.653, 6.643, 7.279, 8.784],
]

_OFFS1 = [-2, -2, -2, -1, -1, 0, -3, 0, -3]
_OFFS2 = [0, 1, 2, 1, 2, 2, 3, 3, 0]


def dist2sid(d: float) -> int:
    return int(d * d * (100.0 / 16.0) + 0.5)


def emit_lines() -> list[str]:
    k = len(_MEANS_ANGSTROM)
    d = len(_OFFS1)
    assert all(len(row) == d for row in _MEANS_ANGSTROM)
    sids = [dist2sid(x) for row in _MEANS_ANGSTROM for x in row]
    lines = [
        f"sec\t{k}",
        f"dim\t{d}",
        "offs1\t" + str(d) + "\t" + "\t".join(str(x) for x in _OFFS1),
        "offs2\t" + str(d) + "\t" + "\t".join(str(x) for x in _OFFS2),
        "mean\t" + str(k * d) + "\t" + "\t".join(str(x) for x in sids),
    ]
    return lines


def main() -> None:
    for line in emit_lines():
        print(line)


if __name__ == "__main__":
    main()
