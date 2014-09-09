#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
PBclust.py 

2013 - P. Poulain, A. G. de Brevern 
"""

#===============================================================================
# load modules
#===============================================================================
import sys
import os
import subprocess
import optparse
# optparse in deprecated since Python 2.7 and has been replaced by argparse
# however many Python installations are still using Python < 2.7

import numpy

import PBlib as PB

#===============================================================================
# functions
#===============================================================================
def compute_score_by_position(score_mat, seq1, seq2):
    """computes similarity score between two sequences"""
    assert len(seq1) == len(seq2), "sequences have different sizes\n%s\n%s" %(seq1, seq2)
    score = []
    for pb1, pb2 in zip(seq1, seq2):
        # score is 0 for Z (dummy PB)
        if "z" in [pb1.lower(), pb2.lower()]:
            score.append(0)
        else:
            score.append( score_mat[PB.NAMES.index(pb1)][PB.NAMES.index(pb2)] )
    return score


#===============================================================================
# main - program starts here
#===============================================================================

#-------------------------------------------------------------------------------
# get options
#-------------------------------------------------------------------------------
parser = optparse.OptionParser(
    usage="%prog -f file.PB.fasta [options] -o output_root_name",
    version="1.0")
# mandatory arguments
mandatory_opts = optparse.OptionGroup(
    parser,
    'Mandatory arguments')
mandatory_opts.add_option("-f", action="append", type="string", 
help="name(s) of the PBs file (in fasta format)")
mandatory_opts.add_option("-o", action="store", type="string", 
help="root name for results")
parser.add_option_group(mandatory_opts)
# optional arguments
optional_opts = optparse.OptionGroup(
    parser,
    'Optional arguments')
optional_opts.add_option("--residue-shift", action="store", type="int",
    dest = "residue_shift", help="shift to adjust residue number")
optional_opts.add_option("--clusters", action="store", type="int",
    dest = "clusters_nb", help="number of clusters wanted")  
optional_opts.add_option("--compare", action="store_true", default=False,
    dest = "compare", help="compare the first sequence versus all others")
parser.add_option_group(optional_opts)
# get all parameters
(options, args) = parser.parse_args()

#-------------------------------------------------------------------------------
# check options
#-------------------------------------------------------------------------------
if not options.f:
    parser.print_help()
    parser.error("option -f is mandatory")

if not options.o:
    parser.print_help()
    parser.error("option -o is mandatory")

if options.residue_shift and options.residue_shift < 0:
	parser.error("residue shift must be positive")
if options.residue_shift:
    residue_shift = options.residue_shift
else:
    residue_shift = 0

if options.clusters_nb and options.clusters_nb <= 0:
    parser.error("number of clusters must be strictly positive")

if not options.clusters_nb:
    # set default clusters number to 5
    options.clusters_nb = 5
    
#-------------------------------------------------------------------------------
# check input files
#-------------------------------------------------------------------------------
for name in options.f:
    if not os.path.isfile(name):
        sys.exit("%s does not appear to be a valid file.\nBye" % name)

#-------------------------------------------------------------------------------
# read PBs files
#-------------------------------------------------------------------------------
header_lst = []
seq_lst = []
for name in options.f:
    header, seq =  PB.read_fasta(name)
    header_lst += header
    seq_lst += seq

pb_seq = numpy.array( zip(header_lst, seq_lst) )

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
    ref_name = pb_seq[0,0]
    ref_seq = pb_seq[0,1]
    mini = numpy.min(substitution_mat)
    maxi = numpy.max(substitution_mat)
    # normalize substitution matrix between 0 and 9
    # 0 -> similar PBs
    # 9 -> different PBs
    substitution_mat_modified = (substitution_mat + abs(mini))/(maxi - mini)
    substitution_mat_modified = 9 * (1 - substitution_mat_modified)
    substitution_mat_modified = substitution_mat_modified.astype(int)
    # set diagonal to 0
    for idx in xrange(len(substitution_mat_modified)):
        substitution_mat_modified[idx,idx] = 0
    print "Normalized substitution matrix (between 0 and 9)"
    print substitution_mat_modified
    print "Compare first sequence (%s) with others" % ref_name
    for target in pb_seq[1:,]:
        header = "%s vs %s" % (ref_name, target[0])
        score_lst = compute_score_by_position(substitution_mat_modified, ref_seq, target[1] )
        seq = "".join([str(s) for s in score_lst])
        PB.write_fasta(compare_file_name, seq, header)
    print "wrote %s" % compare_file_name
    # name = options.o + ".PB.compare.data"
    # f = open(name, "w")
    # for idx, score in enumerate(score_lst):
    #     f.write( "%4d  %d\n" % (idx + 1 + residue_shift, score) )
    # f.close()
    # print "wrote %s" % (name)
    sys.exit(0)

# change sequence name for a better input in R
seq_names = {}
for i in xrange(len(pb_seq)):
    new_name = "seq%d" % i
    seq_names[new_name] = pb_seq[i, 0]
    pb_seq[i, 0] = new_name

#-------------------------------------------------------------------------------
# compute distance of all sequences against all
#-------------------------------------------------------------------------------
distance_mat = numpy.empty((len(pb_seq), len(pb_seq)), dtype='float')

print "Building distance matrix"
# get similarity score
for i in xrange(len(pb_seq)):
    sys.stdout.write("\r%.f%%" % (float(i+1)/len(pb_seq)*100))
    sys.stdout.flush()
    for j in xrange(i, len(pb_seq)):
        score = sum( compute_score_by_position(substitution_mat, pb_seq[i, 1], pb_seq[j, 1]) )
        distance_mat[i, j] = score
        distance_mat[j, i] = score 
print ""

# set equal the diagonal
diag_mini =  numpy.min([distance_mat[i, i] for i in xrange(len(pb_seq))])
for i in xrange(len(pb_seq)):
    distance_mat[i, i] = diag_mini

# convert similarity score to normalized distance between 0 and 1
# dist = 1 means sequences are very different
# dist = 0 means sequences are identical
# dist = 1 - (score + abs(min)/(max - min)

mini = numpy.min(distance_mat)
maxi = numpy.max(distance_mat)
distance_mat = 1 - (distance_mat + abs(mini))/(maxi - mini)

numpy.set_printoptions(threshold=numpy.inf, precision = 3, linewidth = 100000)
output_mat_str = numpy.array_str(distance_mat).translate(None, '[]')

# add sequence labels
output_mat_str = " ".join(pb_seq[:,0])+"\n"+output_mat_str

# write distance matrix
name = options.o + ".PB.dist"
f = open(name, "w")
f.write(output_mat_str)
f.close()
print "wrote", name

# build R script
#-------------------------------------------------------------------------------
# https://github.com/alevchuk/hclust-fasta/blob/master/003-hclust
# and 
# http://www.biostars.org/p/11987/
# data
R_script="""
connector = textConnection("%s")

distances = read.table(connector, header = TRUE)
rownames(distances) = colnames(distances)

clusters = cutree(hclust(as.dist(distances)), k = %d)
distances = as.matrix(distances)

# function to find medoid in cluster i
clust.medoid = function(i, distmat, clusters) {
    ind = (clusters == i)

    if(length(distmat[ind, ind]) == 1){
        names(clusters[ind])
    } else {
        names(which.min(rowSums( distmat[ind, ind] )))
        # c(min(rowMeans( distmat[ind, ind] )))
    }
}

medoids = sapply(unique(clusters), clust.medoid, distances, clusters)

cat("seq_id", names(clusters), "\n")
cat("cluster_id", clusters, "\n")
cat("medoid_id", medoids)
""" % (output_mat_str, options.clusters_nb)


# execute R script
#-------------------------------------------------------------------------------
command="R --vanilla --slave"
proc = subprocess.Popen(command, shell = True, 
stdout = subprocess.PIPE, stderr = subprocess.PIPE, stdin = subprocess.PIPE)
(out, err) = proc.communicate(R_script)
if err:
    print "ERROR:", err
code = proc.wait()
if code:
    print "ERROR: exit code != 0"
    print "exit code:", code
else:
    print "R clustering: OK"

# only 3 lines of output are expected
if len(out.split("\n")) != 3:
    sys.exit("ERROR: wrong R ouput")

seq_id, cluster_id, medoid_id = out.split("\n")
seq_id = seq_id.split()[1:]
cluster_id = cluster_id.split()[1:]
medoid_id = medoid_id.split()[1:]

# count number of sequences in clusters
cluster_count = {}
for idx in cluster_id:
    cluster_count[idx] = cluster_count.get(idx, 0) + 1
for idx in sorted(cluster_count):
    print "cluster %3s: %5d sequences (%3d%%)" %(idx, cluster_count[idx], 1.0*cluster_count[idx]/len(seq_lst)*100)


name = options.o + ".PB.clust"
f = open(name, "w")
for seq, cluster in zip(seq_id, cluster_id):
    f.write('SEQ_CLU  "%s"  %s \n' %(seq_names[seq], cluster))
for idx, med in enumerate(medoid_id):
    f.write('MED_CLU  "%s"  %d \n' %(seq_names[med], idx+1))
f.close()
print "wrote", name
