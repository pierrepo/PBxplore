#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

# Standard module
import subprocess
import sys

# Third-party module
import numpy

# Local module
from . import utils


class RError(RuntimeError):
    """
    Exception raised when something fails with a R script.
    """
    pass


def distance_matrix(sequences, substitution_mat):
    """
    Compute distances of all sequences against all the others

    The substitution matrix is expected to be expressed as similarity scores.
    """
    distance_mat = numpy.empty((len(sequences), len(sequences)), dtype='float')

    print("Building distance matrix")
    # Get similarity score
    for i, seqA in enumerate(sequences):
        sys.stdout.write("\r%.f%%" % (float(i+1)/len(sequences)*100))
        sys.stdout.flush()
        for j, seqB in enumerate(sequences[i:], start=i):
            score = utils.substitution_score(substitution_mat, seqA, seqB)
            distance_mat[i, j] = score
            distance_mat[j, i] = score
    print("")
    # Set the diagonal equal to its maximum value
    diag_maxi = numpy.max(distance_mat.diagonal())
    numpy.fill_diagonal(distance_mat, diag_maxi)
    # Convert similarity score into a distance
    mini = numpy.min(distance_mat)
    maxi = numpy.max(distance_mat)
    # Compute distance
    distance_mat = 1 - (distance_mat - mini)/(maxi - mini)
    # Check distance values are in expected range
    assert(numpy.min(distance_mat) >= 0.0)
    assert(numpy.max(distance_mat) <= 1.0)
    assert(numpy.sum(distance_mat.diagonal()) == 0.0)
    return distance_mat


def _matrix_to_str(distance_mat):
    numpy.set_printoptions(threshold=numpy.inf, precision=3, linewidth=100000)
    output_mat_str = numpy.array_str(distance_mat).replace('[', '').replace(']', '')
    return output_mat_str


def hclust(distance_mat, nclusters, method='ward'):
    """
    Hierachical clustering using R

    Parameters
    ----------
    distance_mat : 2D numpy array
        Distance matrix
    nclusters : int
        Number of cluster to build
    method : str
        Aggregation method for the clustering algorithm; must be a value
        valid for R hclust function

    Returns
    -------
    cluster_id : list of int
        Cluster ID for each item; starts at 1
    medoid_id : list of int
        Index of the medoid for each cluster

    Exceptions
    ----------
    RError : something wrong happened with R
    """
    # Convert the distance matrix into a string readable by R
    output_mat_str = _matrix_to_str(distance_mat)
    # Build the R script
    R_script = """
    connector = textConnection("{matrix}")
    distances = read.table(connector, header=FALSE)
    rownames(distances) = colnames(distances)

    clusters = cutree(hclust(as.dist(distances), method="{method}"), k={clusters})
    distances = as.matrix(distances)

    # function to find medoid in cluster i
    # http://www.biostars.org/p/11987/
    clust.medoid = function(i, distmat, clusters) {{
        ind = (clusters == i)

        if(length(distmat[ind, ind]) == 1){{
            names(clusters[ind])
        }} else {{
            names(which.min(rowSums( distmat[ind, ind] )))
            # c(min(rowMeans( distmat[ind, ind] )))
        }}
    }}

    medoids = sapply(unique(clusters), clust.medoid, distances, clusters)

    cat("cluster_id", clusters, "\n")
    cat("medoid_id", medoids)
    """.format(matrix=output_mat_str, clusters=nclusters, method=method)

    # Execute the R script
    command = "R --vanilla --slave"
    proc = subprocess.Popen(command, shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            stdin=subprocess.PIPE)
    (out, err) = proc.communicate(R_script.encode('utf-8'))
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    if err:
        raise RError("{0}".format(err))
    code = proc.wait()
    if code:
        raise RError("R returned with code {}".format(code))

    # Extract the output of the R script
    # only 2 lines are expected
    if len(out.split("\n")) != 2:
        raise RError("unexpected R ouput")
    cluster_id, medoid_id = out.split("\n")
    # As the input table is provided without headers, the sequences are named
    # V1, V2... with indices starting at 1. To get a integer index starting at
    # 0 from a sequence name, we need to remove the V prefix and to substract 1
    # from the remaining number. This applies to medoid_id that relies on the
    # sequence name, but not to cluster_id that is already a list of integers.
    cluster_id = [int(x) for x in cluster_id.split()[1:]]
    medoid_id = [int(x[1:]) - 1 for x in medoid_id.split()[1:]]

    return cluster_id, medoid_id
