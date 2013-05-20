#! /usr/bin/env python

"""
PBassign.py read Protein Data Bank (PDB) structures
and assigns protein blocs (PBs).

2012 - P. Poulain, A. G. de Brevern 
"""


#===============================================================================
# load modules
#===============================================================================
import optparse 
# optparse in deprecated since Python 2.7 and has been replaced by argparse
# however many Python installations are steal using Python < 2.7
import os
import sys
import numpy 
import math
import glob

#===============================================================================
# data
#===============================================================================
# Protein Blocks angle definition
PB_DATA = """
#PB psi(n-2) phi(n-1)  psi(n-1)   phi(n)   psi(n)   phi(n+1)  psi(n+1)  phi(n+2) 
a    41.14      75.53     13.92   -99.80   131.88     -96.27    122.08    -99.68  
b   108.24     -90.12    119.54   -92.21   -18.06    -128.93    147.04    -99.90  
c   -11.61    -105.66     94.81  -106.09   133.56    -106.93    135.97   -100.63 
d   141.98    -112.79    132.20  -114.79   140.11    -111.05    139.54   -103.16 
e   133.25    -112.37    137.64  -108.13   133.00     -87.30    120.54     77.40   
f   116.40    -105.53    129.32   -96.68   140.72     -74.19    -26.65    -94.51  
g     0.40     -81.83      4.91  -100.59    85.50     -71.65    130.78     84.98   
h   119.14    -102.58    130.83   -67.91   121.55      76.25     -2.95    -90.88  
i   130.68     -56.92    119.26    77.85    10.42     -99.43    141.40    -98.01  
j   114.32    -121.47    118.14    82.88  -150.05     -83.81     23.35    -85.82  
k   117.16     -95.41    140.40   -59.35   -29.23     -72.39    -25.08    -76.16  
l   139.20     -55.96    -32.70   -68.51   -26.09     -74.44    -22.60    -71.74  
m   -39.62     -64.73    -39.52   -65.54   -38.88     -66.89    -37.76    -70.19  
n   -35.34     -65.03    -38.12   -66.34   -29.51     -89.10     -2.91     77.90   
o   -45.29     -67.44    -27.72   -87.27     5.13      77.49     30.71    -93.23  
p   -27.09     -86.14      0.30    59.85    21.51     -96.30    132.67    -92.91 
"""

# line width for fasta format
FASTA_WIDTH = 60

#===============================================================================
# classes
#===============================================================================
class PdbAtomCls:
    """class for atoms in PDB format"""
    def __init__(self):
        """default constructor"""
        self.id = 0
        self.name = None
        self.resalt = ' '
        self.resname = None
        self.chain = None
        self.resid = 0
        self.resins = ''
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.occupancy = 0.0
        self.tempfactor = 0.0
        self.element = '  '
        self.charge =  '  '
        
    def read(self, line):
        """read ATOM data from a PDB file line"""
        self.id = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21:22].strip()
        self.resid = int(line[22:26].strip())
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip()) 
        
    def __repr__(self):
        """representation for atom"""
        return 'atom %4d %4s in %4d %3s' \
        % (self.id, self.name, self.resid, self.resname)
        
    def format(self):
        """atom data formated as PDB"""
        return '%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s' \
        % ('ATOM  ', self.id, self.name, self.resalt,self.resname, self.chain, 
        self.resid, self.resins, self.x, self.y, self.z, 
        self.occupancy, self.tempfactor, self.element, self.charge)
        
    def coord(self):
        return [self.x, self.y, self.z]
    

class PdbStructureCls:
    """class for a Protein Data Bank Structure
    - computes dihedral angles
    """
    def __init__(self):
        """default constructor for PDB structure"""
        self.chain = ""
        self.atoms = []
        self.comment = ""
        self.PB = ""

    def add_atom(self, atom):
        """add atom to the structure"""
        # update chain name
        if not self.atoms:
            self.chain = atom.chain
        elif self.chain != atom.chain:
            print "WARNING: several chains in the same structure"
        # add atom to structure
        self.atoms.append(atom)

    def clean(self):
        """clean up chain, atoms and comment"""
        self.chain = ""
        self.atoms = []
        self.comment = ""
        
    def size(self):
        """get number of atoms from a structure"""
        return len(self.atoms)
               
    def get_all_dihedral(self):
        """compute phi and psi angles for a protein"""
        # extract backbone atoms
        backbone = {}
        for atom in self.atoms:
            if atom.name in ["CA", "C", "O", "N"]:
                resid = atom.resid
                if resid in backbone:
                    backbone[resid][atom.name] = atom
                else:
                    backbone[resid] = {atom.name: atom}
        
        # get dihedrals 
        phi_psi_angles = {}
        for res in sorted(backbone.iterkeys()):
            # phi : C(i-1) - N(i) - CA(i) - C(i)
            try:
                phi = get_dihedral(backbone[res-1]["C"], backbone[res]["N"], backbone[res]["CA"], backbone[res]["C"])
            except:
                phi = None
            # psi : N(i) - CA(i) - C(i) - N(i+1)
            try:
                psi = get_dihedral(backbone[res]["N"], backbone[res]["CA"], backbone[res]["C"], backbone[res+1]["N"])
            except:
                psi = None
            #print res, phi, psi
            phi_psi_angles[res] = {"phi":phi, "psi":psi}
        return phi_psi_angles
                
       
#===============================================================================
# functions
#===============================================================================
def get_dihedral(atomA, atomB, atomC, atomD) :
    """compute dihedral angle between 4 atoms (A, B, C, D)
    output is in degree in the range -180 / +180
    """
    
    # get numpy objects
    A = numpy.array(atomA.coord())
    B = numpy.array(atomB.coord())
    C = numpy.array(atomC.coord())
    D = numpy.array(atomD.coord())
    
       
    # vectors
    AB = B - A 
    BC = C - B 
    CD = D - C 

    # normal vectors
    n1 = numpy.cross(AB, BC)
    n2 = numpy.cross(BC, CD)

    # normalize normal vectors
    n1 /= numpy.linalg.norm(n1)
    n2 /= numpy.linalg.norm(n2)
    
    # angle between normals
    cosine = numpy.sum(n1*n2) / (numpy.linalg.norm(n1) * numpy.linalg.norm(n2))
    torsion = math.acos(cosine)

    # convert radion to degree
    torsion = torsion * 180.0 / math.pi 

    # find if the torsion is clockwise or counterclockwise
    #if numpy.sum(n1 * CD) < 0.0:
    if numpy.dot(n1, CD) < 0.0:
        torsion = 360 - torsion
    if torsion == 360.0:
        torsion = 0.0
    
    # get range -180 / +180
    if torsion > 180.0:
        torsion = torsion - 360
    if torsion < -180.0:
        torsion = torsion + 360
   
    return torsion


#-------------------------------------------------------------------------------
def angle_modulo_360(angle):
    """keep angle in the range -180 / +180 [degrees]
    """
    if angle > 180.0:
        return angle - 360.0
    elif angle < -180.0:
        return angle + 360.0
    else:
        return angle
    
#-------------------------------------------------------------------------------
def write_fasta(name, seq, comment):
    """format seq and comment to fasta format
    and write file
    """
    fasta_content  = ">"+comment+"\n"
    fasta_content += "\n".join( [seq[i:i+FASTA_WIDTH] for i in xrange(0, len(seq), FASTA_WIDTH)] )
    fasta_content += "\n"
    f_out = open(name, "a")
    f_out.write(fasta_content)
    f_out.close()

#-------------------------------------------------------------------------------
def write_phipsi(name, torsion, com):
    """save phi and psi angles
    """
    f_out = open(name, "a")
    for res in sorted(torsion.iterkeys()):
        try:
            phi_str = "%8.2f" % torsion[res]["phi"]
        except:
            phi_str = "    None"
        try:
            psi_str = "%8.2f" % torsion[res]["psi"]
        except:
            psi_str = "    None"
        f_out.write("%s %6d %s %s \n" % (com, res, phi_str, psi_str))
    f_out.close()

#-------------------------------------------------------------------------------
def write_flat(name, seq):
    """write flat sequence to file 
    """
    f_out = open(name, "a")
    f_out.write(seq + "\n")
    f_out.close()

#-------------------------------------------------------------------------------
def clean_file(name):
    """clean existing file
    """
    if os.path.exists(name):
        os.remove(name)

#-------------------------------------------------------------------------------
def PB_assign(pb_ref, structure, comment):
    """assign Protein Blocks (PB) from phi and psi angles
    """
    # get phi and psi angles from structure
    dihedrals = structure.get_all_dihedral()
    # write phi and psi angles
    if options.phipsi:
        write_phipsi(phipsi_name, dihedrals, comment)

    pb_seq = ""
    # iterate over all residues
    for res in sorted(dihedrals.iterkeys()):
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
        except:
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

    # write PBs in fasta file
    write_fasta(fasta_name, pb_seq, comment)
    
    # write PBs in flat file
    if options.flat:
        write_flat(flat_name, pb_seq)
 
    print "PBs assigned for", comment 
             
#-------------------------------------------------------------------------------
# vertorize function
#-------------------------------------------------------------------------------
angle_modulo_360_vect = numpy.vectorize(angle_modulo_360)

#===============================================================================
# MAIN - program starts here
#===============================================================================

#-------------------------------------------------------------------------------
# manage parameters
#-------------------------------------------------------------------------------
parser = optparse.OptionParser(
    usage="%prog [options] -p file.pdb|dir [-p file2.pdb] -o output_root_name",
    version="1.0")
# mandatory arguments
mandatory_opts = optparse.OptionGroup(
    parser,
    'Mandatory arguments')
mandatory_opts.add_option("-p", action="append", type="string", 
    help="name of pdb file or directory containing pdb files")
mandatory_opts.add_option("-o", action="store", type="string", 
    help="root name for results")
parser.add_option_group(mandatory_opts)
# optional arguments
optional_opts = optparse.OptionGroup(
    parser,
    'Optional arguments')
optional_opts.add_option("--phipsi", action="store_true", default=False,
    help="print phi and psi angle")
optional_opts.add_option("--flat", action="store_true", default=False,
    help="output with one sequence per line")
parser.add_option_group(optional_opts)
# get all parameters
(options, args) = parser.parse_args()

# check options
if not options.p:
    parser.print_help()
    parser.error("options -p is mandatory")

if not options.o:
    parser.print_help()
    parser.error("option -o is mandatory")

#-------------------------------------------------------------------------------
# check files
#-------------------------------------------------------------------------------
pdb_name_lst = []

for name in options.p:
    if os.path.isfile(name):
        pdb_name_lst.append(name)
    elif os.path.isdir(name):
        pdb_name_lst += glob.glob(name + "/*.pdb")
    elif (not os.path.isfile(name) or not os.path.isdir(name)):
        print "%s does not appear to be a valid file or directory" % name

print "%d PDB file(s) to process" % (len(pdb_name_lst))
if not pdb_name_lst:
    sys.exit("Nothing to do. Bye.")


#-------------------------------------------------------------------------------
# read PB definitions
#-------------------------------------------------------------------------------
pb_def = {}
for line in PB_DATA.split("\n"):
    if line and "#" not in line:
        items = line.split()
        pb_def[items[0]] = numpy.array([float(items[i]) for i in xrange(1, len(items))])
print "read PB definitions: %d PBs x %d angles " % (len(pb_def), len(pb_def["a"]))

#-------------------------------------------------------------------------------
# prepare fasta file for output
#-------------------------------------------------------------------------------
fasta_name = options.o + ".PB.fasta"
clean_file(fasta_name)

#-------------------------------------------------------------------------------
# prepare phi psi file for output
#-------------------------------------------------------------------------------
if options.phipsi:
    phipsi_name = options.o + ".PB.phipsi"
    clean_file(phipsi_name)
 
#-------------------------------------------------------------------------------
# prepare flat file for output
#-------------------------------------------------------------------------------
if options.flat:
    flat_name = options.o + ".PB.flat"
    clean_file(flat_name)

#-------------------------------------------------------------------------------
# read PDB files
#-------------------------------------------------------------------------------
structure = PdbStructureCls()
model = ""
chain = " "
comment = ""

for pdb_name in pdb_name_lst:
    print pdb_name 
    f_in = open(pdb_name, 'r')
    for line in f_in:
        flag = line[0:6].strip()
        if flag == "MODEL":
            model = line.split()[1]
        if flag == "ATOM":
            atom = PdbAtomCls()
            atom.read(line)
            # assign structure upon new chain
            if structure.size() != 0 and structure.chain != atom.chain:
                PB_assign(pb_def, structure, comment)
                structure.clean()
            # append structure with atom
            structure.add_atom(atom)
            # define structure comment
            # when the structure contains 1 atom
            if structure.size() == 1:
                comment = pdb_name 
                if model:
                    comment += " | model %s" % (model)
                    model = ""
                if atom.chain:
                    comment += " | chain %s" % (atom.chain)
                    atom.chain = ""
        # assign structure after end of model (or chain)
        if structure.size() != 0 and flag in ["TER", "ENDMDL"]:
            PB_assign(pb_def, structure, comment)
            structure.clean()
    # assign last structure
    if structure.size() != 0:
        PB_seq = PB_assign(pb_def, structure, comment)

    f_in.close()   
print "wrote %s" % (fasta_name)
if options.flat:
    print "wrote %s" % (flat_name)
if options.phipsi:
    print "wrote %s" % (phipsi_name)

