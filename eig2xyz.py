#!/usr/bin/env python3

import numpy as np
from argparse import ArgumentParser

atomic_masses = [0.00, 1.008, 4.002, 6.94, 9.012, 10.81, 12.011,
                 14.007, 15.999, 18.998, 20.1797, 22.989, 24.305,
                 26.981, 28.085, 30.973, 32.06, 35.45, 39.948,
                 39.0983, 40.078, 44.955, 47.867, 50.9415, 51.9961,
                 54.938, 55.845, 58.933, 58.6934, 63.546, 65.38,
                 69.723, 72.63, 74.921, 78.971, 79.904, 83.798,
                 85.4678, 87.62, 88.905, 91.224, 92.906, 95.95, 97.0,
                 101.07, 102.905, 106.42, 107.8682, 112.414, 114.818,
                 118.71, 121.76, 127.6, 126.904, 131.293, 132.905,
                 137.327, 138.905, 140.116, 140.907, 144.242, 145.0,
                 150.36, 151.964, 157.25, 158.925, 162.5, 164.93,
                 167.259, 168.934, 173.045, 174.9668, 178.486,
                 180.947, 183.84, 186.207, 190.23, 192.217, 195.084,
                 196.966, 200.592, 204.38, 207.2, 208.98, 209.0,
                 210.0, 222.0, 223.0, 226.0, 227.0, 232.0377, 231.035,
                 238.028, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0,
                 252.0, 257.0, 258.0, 259.0, 262.0, 267.0, 270.0,
                 269.0, 270.0, 270.0, 278.0, 281.0, 281.0, 285.0,
                 286.0, 289.0, 289.0, 293.0, 293.0, 294.0]


def get_args():
    parser = ArgumentParser(description="Convert a GULP .eig file "
                            "into an .xyz file for visualisation with, "
                            "e.g., Jmol.")
    parser.add_argument('filename', type=str, help="The .eig file to read.")
    parser.add_argument('kindex', metavar='k', type=int,
                        help="The index of the k-point to display.")
    parser.add_argument('eigindex', metavar='N', type=int,
                        help="The index of the eigenvector to display.")
    return parser.parse_args()


class EigenvectorModel():
    def __init__(self, filename):
        with open(filename, 'r') as inf:
            natoms = int(inf.readline())

            Z = np.zeros(natoms, dtype=int)
            A = np.zeros(natoms, dtype=float)
            pos = np.zeros((natoms, 3), dtype=float)

            for i in range(natoms):
                atom = inf.readline().split()
                Z[i] = int(atom[0])
                A[i] = atomic_masses[Z[i]]
                pos[i, :] = [float(f) for f in atom[1:4]]

            nkpoints = int(inf.readline())
            nmodes = int(inf.readline())

            eigenvectors = np.zeros((nkpoints, nmodes, natoms, 3), dtype=float)

            for ikpoint in range(nkpoints):
                kpoint_header = inf.readline()
                for imode in range(nmodes):
                    mode_header = inf.readline()
                    frequency = float(inf.readline())
                    for i in range(natoms):
                        trans = inf.readline().split()
                        eigenvectors[ikpoint, imode, i, :] = \
                            [float(f) for f in trans]

        self.Z = Z
        self.A = A
        self.pos = pos
        self.eigenvectors = eigenvectors

        self.natoms = natoms
        self.nkpoints = nkpoints
        self.nmodes = nmodes

    def writexyz(self, filename, ikpoint, imode):
        with open(filename, 'w') as outf:
            outf.write(f"{self.natoms}\n\n")
            for i in range(self.natoms):
                displacement = self.eigenvectors[ikpoint, imode, i, :] / \
                    np.sqrt(self.A[i])
                outf.write(("{:12.6f}"*6 + '\n').format(
                    *self.pos[i, :], *displacement))


if __name__ == '__main__':
    args = get_args()

    model = EigenvectorModel(args.filename)

    if args.kindex > model.nkpoints:
        raise IndexError(f"k-point {args.kindex} requested, but only "
                         "{model.nkpoints} k-points in .eig file!")
    if args.eigindex > model.nmodes:
        raise IndexError(f"mode {args.eigindex} requested, but only "
                         "{model.nmodes} modes in .eig file!")

    outfile = '.'.join(args.filename.split('.'))[:-1] + '.xyz'

    model.writexyz(outfile, args.kindex - 1, args.eigindex - 1)
