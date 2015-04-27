#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python library to handle Protein Blocks

2013 - P. Poulain, A. G. de Brevern 
"""
#===============================================================================
# Modules
#===============================================================================
## Use print as a function for python 3 compatibility
from __future__ import print_function

## standard modules
import os
import subprocess
import sys

## third-party modules
import numpy

#===============================================================================
# Python2/Python3 compatibility
#===============================================================================

# The range function in python 3 behaves as the range function in python 2
# and returns a generator rather than a list. To produce a list in python 3,
# one should use list(range). Here we change range to behave the same in
# python 2 and in python 3. In both cases, range will return a generator.
try:
    range = xrange
except NameError:
    pass

#===============================================================================
# Data
#===============================================================================
# real module directory
__location__ = os.path.realpath(
               os.path.join(os.getcwd(), os.path.dirname(sys.argv[0])) )

# Protein Blocks reference angles
# taken from A. G. de Brevern, C. Etchebest and S. Hazout. 
# "Bayesian probabilistic approach for predicting backbone structures 
# in terms of protein blocks"
# Proteins, 41: 271-288 (2000)
REFERENCES = {
'a': [ 41.14,    75.53,   13.92,   -99.80,   131.88,   -96.27,  122.08,   -99.68],
'b': [108.24,   -90.12,  119.54,   -92.21,   -18.06,  -128.93,  147.04,   -99.90],
'c': [-11.61,  -105.66,   94.81,  -106.09,   133.56,  -106.93,  135.97,  -100.63],
'd': [141.98,  -112.79,  132.20,  -114.79,   140.11,  -111.05,  139.54,  -103.16],
'e': [133.25,  -112.37,  137.64,  -108.13,   133.00,   -87.30,  120.54,    77.40], 
'f': [116.40,  -105.53,  129.32,   -96.68,   140.72,   -74.19,  -26.65,   -94.51],
'g': [  0.40,   -81.83,    4.91,  -100.59,    85.50,   -71.65,  130.78,    84.98], 
'h': [119.14,  -102.58,  130.83,   -67.91,   121.55,    76.25,   -2.95,   -90.88],
'i': [130.68,   -56.92,  119.26,    77.85,    10.42,   -99.43,  141.40,   -98.01],
'j': [114.32,  -121.47,  118.14,    82.88,  -150.05,   -83.81,   23.35,   -85.82],
'k': [117.16,   -95.41,  140.40,   -59.35,   -29.23,   -72.39,  -25.08,   -76.16],
'l': [139.20,   -55.96,  -32.70,   -68.51,   -26.09,   -74.44,  -22.60,   -71.74],
'm': [-39.62,   -64.73,  -39.52,   -65.54,   -38.88,   -66.89,  -37.76,   -70.19], 
'n': [-35.34,   -65.03,  -38.12,   -66.34,   -29.51,   -89.10,   -2.91,    77.90],   
'o': [-45.29,   -67.44,  -27.72,   -87.27,     5.13,    77.49,   30.71,   -93.23],  
'p': [-27.09,   -86.14,    0.30,    59.85,    21.51,   -96.30,  132.67,   -92.91]
}
# PB  psi(n-2) phi(n-1)  psi(n-1)   phi(n)   psi(n)   phi(n+1)  psi(n+1)  phi(n+2) 

# names of the 16 PBs
NAMES = ["a", "b", "c", "d", "e", "f", "g", "h",
           "i", "j", "k", "l", "m", "n", "o", "p"]
NUMBER = len(NAMES)

SUBSTITUTION_MATRIX_NAME = os.path.join(__location__, "PBs_substitution_matrix.dat")


# line width for fasta format
FASTA_WIDTH = 60


#===============================================================================
# Exceptions
#===============================================================================

class SizeError(AssertionError):
    """
    Exception raised when a sequence does not have the expected length.
    """
    pass


class InvalidBlockError(ValueError):
    """
    Exception raised when encounter an invalid protein block.
    """
    def __init__(self, block=None):
        self.block = block

    def __repr__(self):
        if block is None:
            return "Invald block"
        else:
            return "Ivalid block '{}'".format(self.block)

#===============================================================================
# Functions
#===============================================================================


#-------------------------------------------------------------------------------    
def read_fasta(name):
    """
    Read fasta file and output sequences in a list.
    
    Parameters
    ----------
    name : str
        Name of file containing sequences in fasta format.
    
    Returns
    -------
    header_lst : list
        List of headers (str)
    sequence_lst : list
        List of sequences (str)
    
    """
    assert os.path.exists(name), name + ' does not exist'
    sequence_lst = []
    header_lst = []
    header = ""
    sequence = ""
    f_in = open(name, "rt")
    for line in f_in:
        data = line.strip()
        # jump empty lines
        if not data:
            continue
        # store header and sequence when a new header (i.e. sequence) is found
        if sequence and header and data.startswith(">"):
            header_lst.append(header)
            sequence_lst.append(sequence)
            # reset header and sequence
            header = ""
            sequence = ""
        # save header of sequence
        if data.startswith(">"):
            header = data[1:]
        # save sequence
        if ">" not in data:
            sequence += data
    f_in.close()
    # save last sequence
    if header and sequence:
        header_lst.append(header)
        sequence_lst.append(sequence)
    # outputs
    assert len(header_lst) == len(sequence_lst), \
           "cannot read same number of headers and sequences"
    print("read %d sequences in %s" % (len(sequence_lst), name))
    if len(sequence_lst) == 0:
        print("WARNING: %s seems empty of sequence" %(name))
    return header_lst, sequence_lst


def read_several_fasta(input_files):
    """
    Read several fasta files

    Note that each fasta file may contain several sequences.

    Parameters
    ----------
    input_files: a list of fasta file paths.

    Returns
    -------
    pb_name: a list of the headers
    pb_seq: a list of the sequences
    """
    pb_seq = []
    pb_name = []
    for name in input_files:
        header, seq = read_fasta(name)
        pb_name += header
        pb_seq += seq
    return pb_name, pb_seq


def assert_same_size(sequences):
    """
    Raise an exception is all sequence are not the same length.

    Parameters
    ----------
    sequences: a list of sequences

    Exceptions
    ----------
    SizeError: not all the sequences are the same length.
    """
    seq_size = len(sequences[0])
    for seq in sequences:
        if len(seq) != seq_size:
            raise SizeError

#-------------------------------------------------------------------------------
def load_substitution_matrix(name):
    """
    Load PB substitution matrix.
    
    Parameters
    ----------
    name : str
        Name of the file containing the PBs susbtitution matrix.
    
    Returns
    -------
    mat : numpy array
        Array of floats.
    """
    try:
        mat = numpy.loadtxt(name, dtype=float, skiprows=2)
    except:
        sys.exit("ERROR: cannot read %s" % name)
    assert len(mat) == 16, 'wrong substitution matrix size'
    assert len(mat[0]) == 16, 'wrong substitution matrix size'
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            if mat[i][j] != mat[j][i]:
                print(i, j)
                print(mat[i][j], mat[j][i])
                sys.exit("ERROR: matrix is not symetric - idx %i and %i" % (i, j))
    print("read substitution matrix")
    return mat

#-------------------------------------------------------------------------------
def clean_file(name):
    """
    Clean existing file.
    
    Parameters
    ----------
    name : str
        Name of file to remove.
    """
    if os.path.exists(name):
        os.remove(name)

#-------------------------------------------------------------------------------
def write_fasta(name, seq, comment):
    """
    Format seq and comment to fasta format and write file.
    
    Parameters
    ----------
    name : str
        Name of file to write.
    seq : str
        Sequence to format.
    comment : str
        Comment to make header of sequence.
    
    """
    fasta_content  = ">"+comment+"\n"
    fasta_content += "\n".join( [seq[i:i+FASTA_WIDTH] for i in range(0, len(seq), FASTA_WIDTH)] )
    fasta_content += "\n"
    f_out = open(name, "a")
    f_out.write(fasta_content)
    f_out.close()


def count_to_transfac(identifier, count_content):
    """
    Convert a table of PB frequencies into transfac format
    
    http://meme.sdsc.edu/meme/doc/transfac-format.html

    Parameters
    ----------
    identifier : str
        Chain used for the ID property in the output.
    count_content : 
        Content of the count file outputed by PBcount as a list of lines.

    Return
    ------
    The frequency matrix as a string in the transfac format.
    """
    residue_lst = []
    transfac_content  = "ID %s\n" % identifier
    transfac_content += "BF unknown\n"
    transfac_content += "P0" + count_content[0][2:]
    for line in count_content[1:]:
        item = line.split()
        residue = int(item[0])
        residue_lst.append(residue)
        transfac_content += "%05d " % residue + line[5:-1] +  "    X" + "\n"
    transfac_content += "XX\n"
    transfac_content += "//"
    return transfac_content


def assign(dihedrals, pb_ref):
    """
    Assign Protein Blocks.

    Dihedral angles are provided as a dictionnary with one item per residue. The
    key is the residue number, and the value is a dictionnary with phi and psi
    as keys, and the dihedral angles as values.

    The protein block definitions are provided as a dictionnary. Each key is a
    block name, the values are lists of dihedral angles.

    Parameters
    ----------
    dihedrals : dict
        Phi and psi dihedral angles for each residue.
    pb_ref : dict
        The definition of the protein blocks.
    """
    pb_seq = ""
    # iterate over all residues
    for res in sorted(dihedrals):
        angles = []
        # try to get all eight angles required for PB assignement
        try:
            angles.append(dihedrals[res-2]["psi"])
            angles.append(dihedrals[res-1]["phi"])
            angles.append(dihedrals[res-1]["psi"])
            angles.append(dihedrals[res  ]["phi"])
            angles.append(dihedrals[res  ]["psi"])
            angles.append(dihedrals[res+1]["phi"])
            angles.append(dihedrals[res+1]["psi"])
            angles.append(dihedrals[res+2]["phi"])
            # check for bad angles 
            # (error while calculating torsion: missing atoms)
            if None in angles:
                pb_seq += "Z"
                continue 
           
        # cannot get required angles (Nter, Cter or missign residues)
        # -> cannot assign PB
        # jump to next residue
        except KeyError:
            pb_seq += "Z"
            continue
        
        # convert to array
        angles = numpy.array(angles)

        # compare to reference PB angles
        rmsda_lst = {}
        for block in pb_ref:
            diff = pb_ref[block] - angles
            diff2 = angle_modulo_360_vect(diff)
            rmsda = numpy.sum(diff2**2)
            rmsda_lst[rmsda] = block
        pb_seq += rmsda_lst[min(rmsda_lst)]
    return pb_seq


def angle_modulo_360(angle):
    """keep angle in the range -180 / +180 [degrees]
    """
    if angle > 180.0:
        return angle - 360.0
    elif angle < -180.0:
        return angle + 360.0
    else:
        return angle
    

def write_phipsi(name, torsion, com):
    """save phi and psi angles
    """
    f_out = open(name, "a")
    for res in sorted(torsion):
        try:
            phi_str = "%8.2f" % torsion[res]["phi"]
        except TypeError:
            phi_str = "    None"
        try:
            psi_str = "%8.2f" % torsion[res]["psi"]
        except TypeError:
            psi_str = "    None"
        f_out.write("%s %6d %s %s \n" % (com, res, phi_str, psi_str))
    f_out.close()


def write_flat(name, seq):
    """write flat sequence to file 
    """
    f_out = open(name, "a")
    f_out.write(seq + "\n")
    f_out.close()


def count_matrix(pb_seq):
    """
    Count the occurences of each block at each position.

    The occurence matrix has one row per sequence, and one column per block.
    The columns are ordered in as PB.NAMES.

    Parameters
    ----------
    pb_seq: a list of PB sequences.

    Returns
    -------
    The occurence matrix.
    """
    assert_same_size(pb_seq)
    pb_count = numpy.zeros((len(pb_seq[0]),  len(NAMES)))
    for seq in pb_seq:
        for idx, block in enumerate(seq):
            if block in NAMES:
                pb_count[idx, NAMES.index(block)] += 1.0
            elif block not in ["Z", "z"]:
                raise InvalidBlockError(block=block)
    return pb_count


def write_count_matrix(pb_count, outfile, first=1):
    """
    Write a PB occurence matrix in a file.

    Parameters
    ----------
    pb_count: an occurence matrix as a 2D numpy array.
    outfile: an open file where to write the matrix.
    first: the residue number of the first position.
    """
    # write the header (PB names)
    print("    " + "".join(["%6s" % name for name in NAMES]), file=outfile)
    # write the data table
    for residue_idx, residue_pb in enumerate(pb_count):
        print("%-5d" % (residue_idx + first) +
              " ".join("%5d" % i for i in residue_pb), file=outfile)


def compute_score_by_position(score_mat, seq1, seq2):
    """
    Computes substitution score between two sequences position per position

    The substitution score can represent a similarity or a distance depending
    on the score matrix provided. The score matrix should be provided as a 2D
    numpy array with score[i, j] the score to swich the PB at the i-th position
    in PB.NAMES to the PB at the j-th position in PB.NAMES.

    The function returns the result as a list of substitution scores to go from
    `seq1` to `seq2` for each position. Both sequences must have the same
    length.

    ..note:

        The score to move from or to a Z block (dummy block) is always 0.
    """
    assert len(seq1) == len(seq2), "sequences have different sizes:\n{}\nvs\n{}".format(seq1, seq2)
    score = []
    for pb1, pb2 in zip(seq1, seq2):
        # score is 0 for Z (dummy PB)
        if "z" in [pb1.lower(), pb2.lower()]:
            score.append(0)
        else:
            score.append( score_mat[NAMES.index(pb1)][NAMES.index(pb2)] )
    return score


def matrix_to_single_digit(matrix):
    """
    Convert a similarity substitution matrix to a single digit distance
    substitution matrix

    This function takes a substitution matrix expressed as similarity scores,
    and expresses it as integer distances in the range [0; 9]. Where 0 means
    similar PBs, and 9 means different PBs.

    Such single digit substitution matrix is convenient to display sequences
    of substitution scores.

    ..note:

            If a substitution matrix expressed with distances is provided
            rather than a matrix expressed with similarity scores, then the
            returned matrix is expressed with similarity scores.
    """
    mini = numpy.min(matrix)
    maxi = numpy.max(matrix)
    # Normalize between 0 and 1
    mat_modified = (matrix + abs(mini))/(maxi - mini)
    # Convert similarity scores to distances and change the range to [0; 9]
    mat_modified = 9 * (1 - mat_modified)
    # Convert to integers
    mat_modified = mat_modified.astype(int)
    # Set diagonal to 0
    for idx in range(len(mat_modified)):
        mat_modified[idx,idx] = 0
    return mat_modified


def distance_matrix(pb_seq, substitution_mat, distance_func):
    """
    Compute distances of all sequences against all the others
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
    print("")

    # set equal the diagonal
    diag_mini =  numpy.min([distance_mat[i, i] for i in range(len(pb_seq))])
    for i in range(len(pb_seq)):
        distance_mat[i, i] = diag_mini
    return distance_mat


def substitution_score(substitution_matrix, seqA, seqB):
    """
    Compute the substitution score to go from `seqA` to `seqB`

    Both sequences must have the same length.

    The score is either expressed as a similarity or a distance depending on
    the substitution matrix.
    """
    return sum(compute_score_by_position(substitution_matrix, seqA, seqB))


def similarity_to_distance(similarity_mat):
    """
    Convert a substitution matrix from similarity scores to distances

    The resulting substitution matrix is expressed as normalized distances
    between 0 and 1, where 0 means the blocks are identical, and 1 means the
    block are the most different.
    """
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


#-------------------------------------------------------------------------------
# vertorize function
#-------------------------------------------------------------------------------
angle_modulo_360_vect = numpy.vectorize(angle_modulo_360)


