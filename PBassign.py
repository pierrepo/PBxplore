#! /usr/bin/env python

"""
PBassign.py read Protein Data Bank (PDB) structures
and assigns protein blocs (PBs).

2012 - P. Poulain, A. G. de Brevern 
"""


#===============================================================================
# load modules
#===============================================================================
from optparse import OptionParser
import os
import sys
import copy
import numpy 
import math
import glob

#===============================================================================
# data
#===============================================================================
# Protein Blocks angle definition
PBdata = """
#PB psi(n-2)  phi(n-1)  psi(n-1)  phi(n)   psi(n)   phi(n+1)  psi(n+1)  phi(n+2)
a   41.14     75.53     13.92     -99.80   131.88   -96.27    122.08    -99.68  
b   108.24    -90.12    119.54    -92.21   -18.06   -128.93   147.04    -99.90  
c   -11.61    -105.66   94.81     -106.09  133.56   -106.93   135.97    -100.63 
d   141.98    -112.79   132.20    -114.79  140.11   -111.05   139.54    -103.16 
e   133.25    -112.37   137.64    -108.13  133.00   -87.30    120.54    77.40   
f   116.40    -105.53   129.32    -96.68   140.72   -74.19    -26.65    -94.51  
g   0.40      -81.83    4.91      -100.59  85.50    -71.65    130.78    84.98   
h   119.14    -102.58   130.83    -67.91   121.55   76.25     -2.95     -90.88  
i   130.68    -56.92    119.26    77.85    10.42    -99.43    141.40    -98.01  
j   114.32    -121.47   118.14    82.88    -150.05  -83.81    23.35     -85.82  
k   117.16    -95.41    140.40    -59.35   -29.23   -72.39    -25.08    -76.16  
l   139.20    -55.96    -32.70    -68.51   -26.09   -74.44    -22.60    -71.74  
m   -39.62    -64.73    -39.52    -65.54   -38.88   -66.89    -37.76    -70.19  
n   -35.34    -65.03    -38.12    -66.34   -29.51   -89.10    -2.91     77.90   
o   -45.29    -67.44    -27.72    -87.27   5.13     77.49     30.71     -93.23  
p   -27.09    -86.14    0.30      59.85    21.51    -96.30    132.67    -92.91
"""

# line width for fasta format
FASTA_WIDTH = 60

#===============================================================================
# classes
#===============================================================================
class PdbAtomCls:
    """class for atoms in PDB format"""
    def __init__(self):
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
        self.id = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21:22].strip()
        self.resid = int(line[22:26].strip())
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip()) 
        
    def __repr__(self):
        return 'atom %4d %4s in %4d %3s' % (self.id, self.name, self.resid, self.resname)
        
    def format(self):
        return '%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s' \
        % ('ATOM  ', self.id, self.name, self.resalt,self.resname, self.chain, self.resid,\
        self.resins, self.x, self.y, self.z, self.occupancy, self.tempfactor, self.element, self.charge)
        
    def coord(self):
        return [self.x, self.y, self.z]
    

class PdbStructureCls:
    """class for a Protein Data Bank Structure
    - computes dihedral angles
    """
    def __init__(self):
        self.chain = ""
        self.atoms = []
        self.comment = ""
        self.PB = ""

    def add_atom(self, atom):
        # update chain name
        if not self.atoms:
            self.chain = atom.chain
        elif self.chain != atom.chain:
            print "WARNING: several chains in the same structure"
        # add atom to structure
        self.atoms.append(atom)

    def clean(self):
        self.chain = ""
        self.atoms = []
        self.comment = ""
        
    def size(self):
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
        phiPsi = {}
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
            phiPsi[res] = {"phi":phi, "psi":psi}
        return phiPsi
                
       
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
    cosine = numpy.sum(n1 * n2) / (numpy.linalg.norm(n1) * numpy.linalg.norm(n2))
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
def PBload(data):
    """load PB definition"""
    PB = {}
    for line in PBdata.split("\n"):
        if line and "#" not in line:
            items = line.split()
            PB[items[0]] = numpy.array([float(items[i]) for i in xrange(1, len(items))])
    return PB

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
    fasta_content += "\n".join( [seq[i:i+FASTA_WIDTH] for i in xrange(0,len(seq),FASTA_WIDTH)] )
    fasta_content += "\n"
    f_out = open(name, "a")
    f_out.write(fasta_content)
    f_out.close()

#-------------------------------------------------------------------------------
def write_phi_psi(name, torsion):
    """save phi and psi angles
    """
    f_out = open(name, "w")
    for res in sorted(torsion.iterkeys()):
        try:
            phi_str = "%8.2f" % torsion[res]["phi"]
        except:
            phi_str = "    None"
        try:
            psi_str = "%8.2f" % torsion[res]["psi"]
        except:
            psi_str = "    None"
        f_out.write("%-d %s %s \n" % (res, phi_str, psi_str))
    f_out.close()
    print "wrote", name

#-------------------------------------------------------------------------------
def write_flat(name, seq):
    """write flat sequence to file 
    """
    f_out = open(name, "a")
    f_out.write(seq + "\n")
    f_out.close()

#-------------------------------------------------------------------------------
def PB_assign(PB, structure, comment):
    """assign Protein Blocks (PB) from phi and psi angles
    """
    # get phi and psi angles from structure
    dihedrals = structure.get_all_dihedral()
    # write phi and psi angles
    if options.phi_psi:
        write_phi_psi(phi_psi_name, dihedrals)

    PB_seq = ""
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
                PB_seq += "Z"
                continue 
           
        # cannot get required angles (Nter, Cter or missign residues)
        # -> cannot assign PB
        # jump to next residue
        except:
            PB_seq += "Z"
            continue
        
        # convert to array
        angles = numpy.array(angles)

        # compare to reference PB angles
        RMSDAlst = {}
        for block in PB:
            diff = PB[block] - angles
            diff2 = angle_modulo_360_vect(diff)
            RMSDA = numpy.sum(diff2**2)
            RMSDAlst[RMSDA] = block
        PB_seq += RMSDAlst[min(RMSDAlst)]

    # write PBs in fasta file
    write_fasta(fasta_name, PB_seq, comment)
    
    # write PBs in flat file
    if options.flat:
        write_flat(flat_name, PB_seq)
 
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
parser = OptionParser(usage="%prog -f file.pdb -d directory -o output_root_name")
parser.add_option("-f", action="store", type="string", dest="pdb_name",
help="name of the pdb file")
parser.add_option("-d", action="store", type="string", dest="dir_name",
help="name of directory that contains pdb files")
parser.add_option("-o", action="store", type="string", dest="out_name",
help="root name for results")
parser.add_option("-p", action="store_true", default=False, dest="phi_psi",
help="[optional] print phi and psi angle")

parser.add_option("--flat", action="store_true", default=False, dest="flat",
help="[optional] output with one sequence per line")

# get all parameters
(options, args) = parser.parse_args()

# check options
if (not options.pdb_name) and (not options.dir_name):
    parser.print_help()
    parser.error("options -f or -d are mandatory")

if not options.out_name:
    parser.print_help()
    parser.error("option -o is mandatory")

out_name = options.out_name

#-------------------------------------------------------------------------------
# check files
#-------------------------------------------------------------------------------
pdb_name_lst = []

if options.pdb_name and not os.path.isfile(options.pdb_name):
    sys.exit("%s does not appear to be a valid file" % (options.pdb_name))

if options.dir_name and not os.path.isdir(options.dir_name):
    sys.exit("%s does not appear to be a valid directory" % (options.dir_name))

if options.pdb_name:
    pdb_name_lst.append(options.pdb_name)
    
if options.dir_name:
    # search all pdb files in directory
    pdb_name_lst += glob.glob(options.dir_name + "/*.pdb")
    # sort filenames alphabetically 
    pdb_name_lst.sort()

print "%d pdb files to process" % (len(pdb_name_lst))

#-------------------------------------------------------------------------------
# read PB definitions
#-------------------------------------------------------------------------------
PB_def = PBload(PBdata)
print "read %d PB definition" % (len(PB_def))

#-------------------------------------------------------------------------------
# prepare fasta file for output
#-------------------------------------------------------------------------------
fasta_name = options.out_name + ".PB.fasta"
if os.path.exists(fasta_name):
    os.remove(fasta_name)

#-------------------------------------------------------------------------------
# prepare phi psi file for output
#-------------------------------------------------------------------------------
if options.phi_psi:
    phi_psi_name = options.out_name + ".phipsi"
    if os.path.exists(phi_psi_name):
        os.remove(phi_psi_name)

#-------------------------------------------------------------------------------
# prepare flat file for output
#-------------------------------------------------------------------------------
if options.flat:
    flat_name = options.out_name + ".PB.flat"
    if os.path.exists(flat_name):
        os.remove(flat_name)


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
                PB_assign(PB_def, structure, comment)
                structure.clean()
            # append structure with atom
            structure.add_atom(atom)
            # define structure comment
            # when the structure contains 1 atom
            if structure.size() == 1:
                comment = pdb_name 
                if model:
                    comment += " | model %s" % (model)
                if atom.chain:
                    comment += " | chain %s" % (atom.chain)
        # assign structure after end of model (or chain)
        if structure.size() != 0 and flag in ["TER", "ENDMDL"]:
            PB_assign(PB_def, structure, comment)
            structure.clean()
    # assign last structure
    if structure.size() != 0:
        PB_seq = PB_assign(PB_def, structure, comment)

    f_in.close()   
print "wrote %s" % (fasta_name)

