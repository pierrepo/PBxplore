#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Cluster protein structures based on their PB sequences.

2013 - P. Poulain, A. G. de Brevern 
"""

#===============================================================================
# Modules
#===============================================================================
## Use print as a function for python 3 compatibility
from __future__ import print_function

## standard modules
import sys
import os
import subprocess
import argparse

## third-party module
import numpy

## local module
import PBlib as PB

## Python2/Python3 compatibility
# The range function in python 3 behaves as the range function in python 2
# and returns a generator rather than a list. To produce a list in python 3,
# one should use list(range). Here we change range to behave the same in
# python 2 and in python 3. In both cases, range will return a generator.
try:
    range = xrange
except NameError:
    pass

def compute_score_by_position(score_mat, seq1, seq2):
    """
    Computes similarity score between two sequences
    """
    assert len(seq1) == len(seq2), "sequences have different sizes:\n{}\nvs\n{}".format(seq1, seq2)
    score = []
    for pb1, pb2 in zip(seq1, seq2):
        # score is 0 for Z (dummy PB)
        if "z" in [pb1.lower(), pb2.lower()]:
            score.append(0)
        else:
            score.append( score_mat[PB.NAMES.index(pb1)][PB.NAMES.index(pb2)] )
    return score


def normalize_matrix(matrix):
    mini = numpy.min(matrix)
    maxi = numpy.max(matrix)
    # normalize substitution matrix between 0 and 9
    # 0 -> similar PBs
    # 9 -> different PBs
    mat_modified = (matrix + abs(mini))/(maxi - mini)
    mat_modified = 9 * (1 - mat_modified)
    mat_modified = mat_modified.astype(int)
    # set diagonal to 0
    for idx in range(len(mat_modified)):
        mat_modified[idx,idx] = 0
    return mat_modified


def user_input():
    #-------------------------------------------------------------------------------
    # get arguments
    #-------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description="Cluster protein structures based on their PB sequences.")

    # mandatory arguments
    parser.add_argument("-f", action="append", required=True, 
        help="name(s) of the PBs file (in fasta format)")
    parser.add_argument("-o", action="store", required=True,
        help="name for results")
    parser.add_argument("-c", action="store", required=True, type=int,
        help="number of wanted clusters")  

    # optional arguments
    parser.add_argument("--compare", action="store_true", default=False,
        dest="compare", help="compare the first sequence versus all others")

    # get all parameters
    options = parser.parse_args()

    #-------------------------------------------------------------------------------
    # check options
    #-------------------------------------------------------------------------------
    if options.c <= 0:
        parser.error("number of clusters must be > 0.")

    #-------------------------------------------------------------------------------
    # check input files
    #-------------------------------------------------------------------------------
    for name in options.f:
        if not os.path.isfile(name):
            sys.exit( "{0}: not a valid file. Bye".format(name) )

    return options


def distance_matrix(pb_seq, substitution_mat, distance_func):
    """
    Compute distance of all sequences against all
    """
    distance_mat = numpy.empty((len(pb_seq), len(pb_seq)), dtype='float')

    print( "Building distance matrix" )
    # get similarity score
    for i in range(len(pb_seq)):
        sys.stdout.write("\r%.f%%" % (float(i+1)/len(pb_seq)*100))
        sys.stdout.flush()
        for j in range(i, len(pb_seq)):
            score = distance_func(substitution_mat, pb_seq[i, 1], pb_seq[j, 1])
            distance_mat[i, j] = score
            distance_mat[j, i] = score 
    print( "" )

    # set equal the diagonal
    diag_mini =  numpy.min([distance_mat[i, i] for i in range(len(pb_seq))])
    for i in range(len(pb_seq)):
        distance_mat[i, i] = diag_mini
    return distance_mat


def substitution_score(substitution_matrix, seqA, seqB):
    return sum(compute_score_by_position(substitution_matrix, seqA, seqB))


def similarity_to_distance(similarity_mat):
    # convert similarity score to normalized distance between 0 and 1
    # dist = 1 means sequences are very different
    # dist = 0 means sequences are identical
    # dist = 1 - (score + abs(min)/(max - min)

    mini = numpy.min(similarity_mat)
    maxi = numpy.max(similarity_mat)
    return 1 - (similarity_mat + abs(mini))/(maxi - mini)


def matrix_to_str(distance_mat, pb_seq):
    numpy.set_printoptions(threshold=numpy.inf, precision = 3, linewidth = 100000)
    output_mat_str = numpy.array_str(distance_mat).replace('[', '').replace(']', '')
    # add sequence labels
    output_mat_str = " ".join(pb_seq[:,0])+"\n"+output_mat_str
    return output_mat_str


def write_distance_matrix(distance_matrix, pb_seq, fname):
    output_mat_str = matrix_to_str(distance_matrix, pb_seq)
    # write distance matrix
    f = open(fname, "w")
    f.write(output_mat_str)
    f.close()
    print("wrote {0}".format(fname))


def hclust(distance_mat, nclusters, pb_seq):
    output_mat_str = matrix_to_str(distance_mat, pb_seq)
    # build R script
    #-------------------------------------------------------------------------------
    # https://github.com/alevchuk/hclust-fasta/blob/master/003-hclust
    # and 
    # http://www.biostars.org/p/11987/
    # data
    R_script="""
    connector = textConnection("{matrix}")

    distances = read.table(connector, header = TRUE)
    rownames(distances) = colnames(distances)

    clusters = cutree(hclust(as.dist(distances), method = "ward"), k = {clusters})
    distances = as.matrix(distances)

    # function to find medoid in cluster i
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

    cat("seq_id", names(clusters), "\n")
    cat("cluster_id", clusters, "\n")
    cat("medoid_id", medoids)
    """.format(matrix=output_mat_str, clusters=nclusters)
    R_script = R_script.encode('utf-8')


    # execute R script
    #-------------------------------------------------------------------------------
    command="R --vanilla --slave"
    proc = subprocess.Popen(command, shell = True, 
    stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin = subprocess.PIPE)
    (out, err) = proc.communicate(R_script)
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    if err:
        print( "ERROR: {0}".format(err) )
    code = proc.wait()
    if code:
        print( "ERROR: exit code != 0" )
        print( "exit code: {0}".format(code) )
    else:
        print( "R clustering: OK" )

    # only 3 lines of output are expected
    if len(out.split("\n")) != 3:
        sys.exit("ERROR: wrong R ouput")

    seq_id, cluster_id, medoid_id = out.split("\n")
    seq_id = [int(x[3:]) for x in seq_id.split()[1:]]
    cluster_id = [int(x) for x in cluster_id.split()[1:]]
    medoid_id = [int(x[3:]) for x in medoid_id.split()[1:]]

    return seq_id, cluster_id, medoid_id


def display_clust_report(cluster_id):
    # count number of sequences in clusters
    cluster_count = {}
    for idx in cluster_id:
        cluster_count[idx] = cluster_count.get(idx, 0) + 1
    for idx in sorted(cluster_count):
        print("cluster {}: {} sequences ({:>2.0f}%)"
              .format(idx, cluster_count[idx],
                      1.0*cluster_count[idx]/len(cluster_id)*100))


def write_clusters(fname, seq_id, cluster_id, medoid_id, seq_names):
    with open(fname, "w") as outfile:
        for seq, cluster in zip(seq_id, cluster_id):
            outfile.write('SEQ_CLU  "{}"  {} \n'.format(seq_names[seq], cluster))
        for idx, med in enumerate(medoid_id):
            outfile.write('MED_CLU  "{}"  {} \n'.format(seq_names[med], idx+1))


def compare_to_first_sequence(headers, sequences, substitution_mat):
    """
    Compare all sequence to the first one
    """
    iheaders = iter(headers)
    isequences = iter(sequences)
    ref_name = next(iheaders)
    ref_seq = next(isequences)
    for target_header, target_seq in zip(iheaders, isequences):
        header = "%s vs %s" % (ref_name, target_header)
        score_lst = compute_score_by_position(substitution_mat, ref_seq, target_seq)
        yield header, score_lst


def compare(header_lst, seq_lst, substitution_mat, fname):
    PB.clean_file(fname)
    ref_name = header_lst[0]
    ref_seq = seq_lst[0]
    substitution_mat_modified = normalize_matrix(substitution_mat)
    print("Normalized substitution matrix (between 0 and 9)")
    print(substitution_mat_modified)
    print("Compare first sequence ({0}) with others".format(ref_name))
    for header, score_lst in compare_to_first_sequence(header_lst, seq_lst,
                                                       substitution_mat_modified):
        seq = "".join([str(s) for s in score_lst])
        PB.write_fasta(fname, seq, header)
    print("wrote {0}".format(fname))


def pbclust_cli():
    options = user_input()
    #-------------------------------------------------------------------------------
    # read PBs files
    #-------------------------------------------------------------------------------
    header_lst, seq_lst = PB.read_several_fasta(options.f)
    pb_seq = numpy.array(list(zip(header_lst, seq_lst)))

    #-------------------------------------------------------------------------------
    # load subtitution matrix
    #-------------------------------------------------------------------------------
    substitution_mat = PB.load_substitution_matrix(PB.SUBSTITUTION_MATRIX_NAME)

    #-------------------------------------------------------------------------------
    # --compare option
    # compare the first sequence (in the fasta file) versus all others
    #-------------------------------------------------------------------------------
    if options.compare:
        compare_file_name = options.o + ".PB.compare.fasta"
        compare(header_lst, seq_lst, substitution_mat, compare_file_name)
        sys.exit(0)

    # change sequence name for a better input in R
    seq_names = {}
    for i in range(pb_seq.shape[0]):
        new_name = "seq%d" % i
        seq_names[new_name] = pb_seq[i, 0]
        pb_seq[i, 0] = new_name

    similarity_mat = distance_matrix(pb_seq, substitution_mat, substitution_score)
    distance_mat = similarity_to_distance(similarity_mat)
    seq_id, cluster_id, medoid_id = hclust(distance_mat, nclusters=options.c, pb_seq=pb_seq)
    distance_fname = options.o + ".PB.dist"
    write_distance_matrix(distance_mat, pb_seq, distance_fname)
    display_clust_report(cluster_id)

    output_fname = options.o + ".PB.clust"
    write_clusters(output_fname, seq_id, cluster_id, medoid_id, header_lst)
    print("wrote {0}".format(output_fname))


if __name__ == '__main__':
    pbclust_cli()
