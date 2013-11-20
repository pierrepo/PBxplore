#! /usr/bin/env python

"""
PBclust.py 

2013 - P. Poulain, A. G. de Brevern 
"""

#===============================================================================
# load modules
#===============================================================================
import PBlib as PB
import optparse
# optparse in deprecated since Python 2.7 and has been replaced by argparse
# however many Python installations are steal using Python < 2.7
import sys
import os
import subprocess
import numpy

#===============================================================================
# functions
#===============================================================================
def check_symetry(mat):
    """check if matrix is symetric"""
    for i in xrange(len(mat)):
        for j in xrange(len(mat[0])):
            if mat[i][j] != mat[j][i]:
                print i, j
                print mat[i][j], mat[j][i]
                sys.exit("ERROR: matrix is not symetric - idx %i and %i" % (i, j))
    print "matrix is symetric"

def compute_score(seq1, seq2):
    """computes similarity score between two sequences"""
    assert len(seq1) == len(seq2), "sequences have different sizes\n%s\n%s" %(seq1, seq2)
    score = 0
    # remove Z (dummy PB)
    seq1 = seq1.lower().replace("z", "")
    seq2 = seq2.lower().replace("z", "")
    for pb1, pb2 in zip(seq1, seq2):
        score += substitution_matrix[PB.NAMES.index(pb1)][PB.NAMES.index(pb2)]
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
pb_seq = []
header = []
seq = []

for name in options.f:
    header, seq = PB.read_fasta(name)
    pb_seq += zip(header, seq)

pb_seq = numpy.array(pb_seq)

# change sequence names for test
for i in xrange(len(pb_seq)):
	pb_seq[i, 0] =  "S" + str(i)

#-------------------------------------------------------------------------------
# load subtitution matrix
#-------------------------------------------------------------------------------
try:
    substitution_matrix = numpy.loadtxt(PB.SUBSTITUTION_MATRIX_NAME, 
                                        dtype=int, skiprows=2)
except:
    sys.exit("ERROR: cannot read %s" % PB.SUBSTITUTION_MATRIX_NAME)

assert len(substitution_matrix) == 16, 'wrong substitution matrix size'
assert len(substitution_matrix[0]) == 16, 'wrong substitution matrix size'
check_symetry(substitution_matrix)


#-------------------------------------------------------------------------------
# compute distance
#-------------------------------------------------------------------------------
distance_mat = numpy.empty((len(pb_seq), len(pb_seq)), dtype='float')

# get similarity score
for i in xrange(len(pb_seq)):
	for j in xrange(i, len(pb_seq)):
		score = compute_score(pb_seq[i, 1], pb_seq[j, 1])
		distance_mat[i, j] = score
		distance_mat[j, i] = score 

# set equal the diagonal
diag_mini =  numpy.min([distance_mat[i, i] for i in xrange(len(pb_seq))])
for i in xrange(len(pb_seq)):
    distance_mat[i, i] = diag_mini

# convert similarity score to distance between 0 and 1
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

print "<debug>"
print out
print "</debug>"
if len(out.split("\n")) != 3:
    sys.exit("ERROR: wrong R ouput")
    
seq_id, cluster_id, medoid_id = out.split("\n")
seq_id = seq_id.split()[1:]
cluster_id = cluster_id.split()[1:]
medoid_id = medoid_id.split()[1:]

print "%d clusters" % (len(medoid_id))

name = options.o + ".PB.clust"
f = open(name, "w")
for seq, cluster in zip(seq_id, cluster_id):
    f.write("SEQ_CLU  %s  %s \n" %(seq, cluster))
for idx, med in enumerate(medoid_id):
    f.write("MED_CLU  %s  %s \n" %(idx+1, med))
f.close()
print "wrote", name




