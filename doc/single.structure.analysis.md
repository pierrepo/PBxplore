# Single structure analysis

**PBassign.py** assigns a PB sequence to a protein structure.

    ./PBassign.py -p test/static/3ICH.pdb -o 3ICH

Output is

    1 PDB file(s) to process
    read PB definitions: 16 PBs x 8 angles
    test/static/3ICH.pdb
    PBs assigned for test/static/3ICH.pdb | chain A
    wrote 3ICH.PB.fasta

Content of `3ICH.PB.fasta` is:

    >test/static/3ICH.pdb | chain A
    ZZccdfbdcdddddehjbdebjcdddddfklmmmlmmmmmmmmnopnopajeopacfbdc
    ehibacehiamnonopgocdfkbjbdcdfblmbccfbghiacdddebehiafkbccddfb
    dcfklgokaccfbdcfbhklmmmmmmmpccdfkopafbacddfbgcddddfbacddddZZ

Note that Protein Blocs assignment is only possible for proteins (as its name suggests). As a consequence, processed PDB files must contain protein structures **only** (please remove any other molecule). In addition, the PDB parser implemented here is pretty straightforward. Be sure your PDB files complies with the [ATOM field](http://www.wwpdb.org/documentation/format33/sect9.html#ATOM) of the [PDB format](http://www.wwpdb.org/documentation/format33/v3.3.html) and that the protein structure is coherent.

PBassign.py can take several options:

|:---------------------|--------------------------------------------------------------------------------------|
| `-h` or `--help`     | shows help message                                                                   |
| `--version`          | shows program version                                                                |
| `-p` **(mandatory)** | defines the name of the PDB file or the name of the directory containing PDB files   |
| `-o` **(mandatory)** | defines the root name for results (do not specify any extension)                     |
| `--phipsi`           | prints phi and psi angles                                                            |
| `--flat`             | outputs will be one sequence per line (not compatible with fasta format)             |

## `-p` option
can be used several times. For instance:

    ./PBassign.py -p test/static/3ICH.pdb -p test/static/1BTA.pdb -p test/static/1AY7.pdb -o test1

Output is:

    3 PDB file(s) to process
    read PB definitions: 16 PBs x 8 angles 
    test/static/3ICH.pdb
    PBs assigned for test/static/3ICH.pdb | chain A
    test/static/1BTA.pdb
    PBs assigned for test/static/1BTA.pdb | chain A
    test/static/1AY7.pdb
    PBs assigned for test/static/1AY7.pdb | chain A
    PBs assigned for test/static/1AY7.pdb | chain B
    wrote test1.PB.fasta

All PB assignments are written in the same output file. If a PDB file contains several chains and/or models, PBs assignments are also written in a single output file. From the previous example, the ouput of `test1.PB.fasta` is:

    >test/static/3ICH.pdb | chain A
    ZZccdfbdcdddddehjbdebjcdddddfklmmmlmmmmmmmmnopnopajeopacfbdc
    ehibacehiamnonopgocdfkbjbdcdfblmbccfbghiacdddebehiafkbccddfb
    dcfklgokaccfbdcfbhklmmmmmmmpccdfkopafbacddfbgcddddfbacddddZZ
    >test/static/1BTA.pdb | chain A
    ZZdddfklonbfklmmmmmmmmnopafklnoiaklmmmmmnoopacddddddehkllmmm
    mngoilmmmmmmmmmmmmnopacdcddZZ
    >test/static/1AY7.pdb | chain A
    ZZbjadfklmcfklmmmmmmmmnnpaafbfkgopacehlnomaccddehjaccdddddeh
    klpnbjadcdddfbehiacddfegolaccdddfkZZ
    >test/static/1AY7.pdb | chain B
    ZZcddfklpcbfklmmmmmmmmnopafklgoiaklmmmmmmmmpacddddddehkllmmm
    mnnommmmmmmmmmmmmmnopacddddZZ

One can also use the `-p` option to provide a directory containing PDB files as an input:

    ./PBassign.py -p test/static/ -o test2

`PBassign.py` will process all PDB files located in the `test/static` directory:

    4 PDB file(s) to process
    read PB definitions: 16 PBs x 8 angles 
    test/static/2LFU.pdb
    PBs assigned for test/static/2LFU.pdb | model 1 | chain A
    PBs assigned for test/static/2LFU.pdb | model 2 | chain A
    PBs assigned for test/static/2LFU.pdb | model 3 | chain A
    PBs assigned for test/static/2LFU.pdb | model 4 | chain A
    PBs assigned for test/static/2LFU.pdb | model 5 | chain A
    PBs assigned for test/static/2LFU.pdb | model 6 | chain A
    PBs assigned for test/static/2LFU.pdb | model 7 | chain A
    PBs assigned for test/static/2LFU.pdb | model 8 | chain A
    PBs assigned for test/static/2LFU.pdb | model 9 | chain A
    PBs assigned for test/static/2LFU.pdb | model 10 | chain A
    test/static/1AY7.pdb
    PBs assigned for test/static/1AY7.pdb | chain A
    PBs assigned for test/static/1AY7.pdb | chain B
    test/static/3ICH.pdb
    PBs assigned for test/static/3ICH.pdb | chain A
    test/static/1BTA.pdb
    PBs assigned for test/static/1BTA.pdb | chain A
    wrote test2.PB.fasta

## `--phipsi` option

generates an additionnal file with the [phi and psi angles](http://en.wikipedia.org/wiki/Dihedral_angle#Dihedral_angles_of_biological_molecules) for each residue.

    ./PBassign.py -p test/static/1BTA.pdb -o 1BTA --phipsi

Output is:

    1 PDB file(s) to process
    read PB definitions: 16 PBs x 8 angles 
    test/static/1BTA.pdb
    PBs assigned for test/static/1BTA.pdb | chain A
    wrote 1BTA.PB.fasta
    wrote 1BTA.PB.phipsi

Content of `1BTA.PB.phipsi` is:

    test/static/1BTA.pdb | chain A      1     None  -171.66 
    test/static/1BTA.pdb | chain A      2  -133.80   153.74 
    test/static/1BTA.pdb | chain A      3  -134.66   157.30 
    test/static/1BTA.pdb | chain A      4  -144.49   118.60 
    test/static/1BTA.pdb | chain A      5  -100.13    92.99 
    test/static/1BTA.pdb | chain A      6   -83.49   104.24 
    test/static/1BTA.pdb | chain A      7   -64.77   -43.25 
    test/static/1BTA.pdb | chain A      8   -44.48   -25.89 
    test/static/1BTA.pdb | chain A      9   -94.91   -47.18 
    test/static/1BTA.pdb | chain A     10   -41.31   133.74 
    [snip]

The first part of the line is the comment also found in the fasta file. The last thee columns are, from left to right, the residue number, the phi angle and the psi angle. The phi angle of the first residue and the psi angle of the last residue cannot be obtained.


## `--flat` option

formats the PBs assignment with one sequence per line. 

    ./PBassign.py -p test/static/1BTA.pdb -o 1BTA --flat

Output is:

    1 PDB file(s) to process
    read PB definitions: 16 PBs x 8 angles 
    test/static/1BTA.pdb
    PBs assigned for test/static/1BTA.pdb | chain A
    wrote 1BTA.PB.fasta
    wrote 1BTA.PB.flat

Content of `1BTA.PB.flat` is:

    ZZdddfklonbfklmmmmmmmmnopafklnoiaklmmmmmnoopacddddddehkllmmmmngoilmmmmmmmmmmmmnopacdcddZZ


