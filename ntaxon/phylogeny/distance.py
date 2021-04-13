from scipy.spatial.distance import pdist as spdist, squareform
from Bio.Phylo.TreeConstruction import DistanceMatrix
from ntaxon.fingerprinting import BinaryMatrix
import numpy as np

def simple_matching(data: BinaryMatrix):
    print("Simple Matching")


def nei(x):
    # x is numpy ndarray rows = allele, cols = species
    # Returns a lower triangular matrix
    if type(x) != np.ndarray:
        x = np.array(x)

    totalalleles, spp = x.shape # get total allele and total species

    d = np.zeros((spp, spp))

    for i in range(1, spp + 1):
        for j in range(0, i - 1):
            s1 = 0.0
            s2 = 0.0
            s3 = 0.0
            for k in range(0, totalalleles):
                s1 += x[k][i - 1] * x[k][j]
                temp = x[k][i - 1]
                s2 += temp * temp

                temp = x[k][j]
                s3 += temp * temp

            if s1 <= 1.0e-20:
                d[j][i - 1] = -1.0
                print("\nWARNING: INFINITE DISTANCE BETWEEN SPECIES ")
                print(f"{i} AND {j}; -1.0 WAS WRITTEN\n")
            else:
                d[i - 1][j] = np.fabs(-1 * np.log(s1 / np.sqrt(s2 * s3)))
    return d


def cavalli(x):
    # x is numpy ndarray rows = allele, cols = species
    # Returns a lower triangular matrix
    if type(x) != np.ndarray:
        x = np.array(x)

    totalalleles, spp = x.shape  # get total allele and total species

    d = np.zeros((spp, spp))

    for i in range(1, spp + 1):
        for j in range(0, i - 1):
            s = 0.0
            for k in range(0, totalalleles):
                f = x[k][i - 1] * x[k][j]
                if f > 0:
                    s += np.sqrt(f)
            # TODO: get loci and df
            # d[i - 1][j] = 4 * (loci - s) / df;
    return d


def reynolds(x):
    # x is numpy ndarray rows = allele, cols = species
    # Returns a lower triangular matrix
    if type(x) != np.ndarray:
        x = np.array(x)

    totalalleles, spp = x.shape  # get total allele and total species

    d = np.zeros((spp, spp))

    for i in range(1, spp + 1):
        for j in range(0, i - 1):
            s1 = 0.0
            s2 = 0.0
            for k in range(0, totalalleles):
                temp = x[k][i - 1] - x[k][j]
                s1 += temp * temp
                s2 += x[k][i - 1] * x[k][j]
            # TODO: get loci
            # d[i - 1][j] = s1 / (loci * 2 - 2 * s2);
    return d


def squareform_to_dist(X, names = None):
    """
    Converts a matrix X of scipy squareform or vector-form distance vecto
    to Biopython DistanceMatrix

    Parameters
    ----------
    X : array_like
        Either a condensed or redundant distance matrix.

    names: list
        List of sample names

    Returns
    -------
    Y : DistanceMatrix
        Biopython Distance Matrix
    """

    X = np.ascontiguousarray(X)
    s = X.shape

    def _convert(m, names):
        lwtm = [] # Convert to lower triangular form
        for i in range(0, len(m)):
            j = i + 1
            lwtm.append(m[i,:j].tolist())
        if names is None:
            n = [f"S{n}" for n in range(1, len(m) + 1)]
        return DistanceMatrix(names=n, matrix=lwtm)

    if len(s) == 1:
        # Distance Vector
        return _convert(squareform(X), names)
    elif len(s) == 2:
        if s[0] != s[1]:
            raise ValueError('The matrix argument must be square.')
        # mxm Matrix
        return _convert(X, names)
