# Introduction to Protein Blocks

Protein Blocks are structural prototypes defined by [de Brevern](http://www.dsimb.inserm.fr/~debrevern/index.php) *et al* in 2000 [1]. Their main interest is to modelize the 3-dimensional local structure of the protein backbone into a 1-dimension sequence. In principle, any conformation of any amino acid could be represented by one of the sixteen available Protein Blocks (see Figure 1).

![PBs](img/PBs.jpg "PBs")

**Figure 1.** Schematic representation of the sixteen protein blocks, labeled from *a* to *p* (Creative commons CC BY).

For instance, the 3D structure of the barstar protein in Figure 2 is represented in term of Protein Blocks:

    ZZdddfklpcbfklmmmmmmmmnopafklgoiaklmmmmmmmmpacddddddehkllmmmmnnomm
    mmmmmmmmmmmmnopacddddZZ

![Barstar protein (PDB 1AY7)](img/1AY7_B.png "Barstar protein (PDB 1AY7)")

**Figure 2.** 3D representation of the barstar protein (PDB ID [1AY7](http://www.rcsb.org/pdb/explore/explore.do?pdbId=1AY7), chain B) (Creative commons CC BY)

The conformations of the 89 residues are translated into a sequence of 89 protein blocks. Note that "Z" is a default name given to amino acid for which a protein block cannot be assigned. Indeed, the assignation procedure for a given residue requires the conformation of the two residues placed before and the two residues placed after the residue under consideration. Therefore, a protein block cannot be assigned to the two first (N-termini) and two last (C-termini) residues of a polypeptide chain.


