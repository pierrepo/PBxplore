PBassign
==========

``PBassign`` assigns a PB sequence to a protein structure.

.. note:: The following examples use ``PBdata`` and the demo files.
          See :ref:`Demo files <demo>` for more information.


Example
-------

.. code-block:: bash

    $ PBassign -p `PBdata`/3ICH.pdb -o 3ICH
    Read 1 chain(s) in demo/3ICH.pdb
    wrote 3ICH.PB.fasta

Content of `3ICH.PB.fasta`: ::

    >demo1/3ICH.pdb | chain A
    ZZccdfbdcdddddehjbdebjcdddddfklmmmlmmmmmmmmnopnopajeopacfbdc
    ehibacehiamnonopgocdfkbjbdcdfblmbccfbghiacdddebehiafkbccddfb
    dcfklgokaccfbdcfbhklmmmmmmmpccdfkopafbacddfbgcddddfbacddddZZ

Note that Protein Blocs assignment is only possible for proteins (as its name suggests).
As a consequence, processed PDB files must contain protein structures **only** (please remove any other molecule).
In addition, the PDB parser implemented here is pretty straightforward.
Be sure your PDB files complies with the `ATOM field <http://www.wwpdb.org/documentation/format33/sect9.html#ATOM>`_
of the `PDB format <http://www.wwpdb.org/documentation/format33/v3.3.html) and that the protein structure is coherent>`_.


Usage
-----

Here’s the ``PBassign`` help text. ::

    Usage: PBassign [options] -p file.pdb|dir [-p file2.pdb] -o output_root_name -g gro_file -x xtc_file

    Options:
      --version   show program's version number and exit
      -h, --help  show this help message and exit

      Mandatory arguments:
        -p P      name of pdb file or directory containing pdb files
        -o O      root name for results
        -x X      name of xtc file (Gromacs)
        -g G      name of gro file (Gromacs)

      Optional arguments:
        --phipsi  writes phi and psi angle
        --flat    writes one PBs sequence per line


`-p` option
```````````

can be used several times. For instance:

.. code-block:: bash

    $ PBassign -p `PBdata`/3ICH.pdb -p `PBdata`/1BTA.pdb -p `PBdata`/1AY7.pdb -o test1
    3 PDB file(s) to process
    Read 1 chain(s) in demo/3ICH.pdb
    Read 1 chain(s) in demo/1BTA.pdb
    Read 2 chain(s) in demo/1AY7.pdb
    wrote test1.PB.fasta


All PB assignments are written in the same output file. If a PDB file contains several chains
and/or models, PBs assignments are also written in a single output file.
From the previous example, the ouput of `test1.PB.fasta` is: ::

    >demo/3ICH.pdb | chain A
    ZZccdfbdcdddddehjbdebjcdddddfklmmmlmmmmmmmmnopnopajeopacfbdc
    ehibacehiamnonopgocdfkbjbdcdfblmbccfbghiacdddebehiafkbccddfb
    dcfklgokaccfbdcfbhklmmmmmmmpccdfkopafbacddfbgcddddfbacddddZZ
    >demo/1BTA.pdb | chain A
    ZZdddfklonbfklmmmmmmmmnopafklnoiaklmmmmmnoopacddddddehkllmmm
    mngoilmmmmmmmmmmmmnopacdcddZZ
    >demo/1AY7.pdb | chain A
    ZZbjadfklmcfklmmmmmmmmnnpaafbfkgopacehlnomaccddehjaccdddddeh
    klpnbjadcdddfbehiacddfegolaccdddfkZZ
    >demo/1AY7.pdb | chain B
    ZZcddfklpcbfklmmmmmmmmnopafklgoiaklmmmmmmmmpacddddddehkllmmm
    mnnommmmmmmmmmmmmmnopacddddZZ


One can also use the `-p` option to provide a directory containing PDB files as an input.
``PBassign`` will process all PDB files located in the `PBdata` directory:

.. code-block:: bash

    $ PBassign -p `PBdata`/ -o test2
    8 PDB file(s) to process
    Read 2 chain(s) in demo/1AY7.pdb
    Read 90 chain(s) in demo/psi_md_traj_1.pdb
    Read 10 chain(s) in demo/2LFU.pdb
    Read 90 chain(s) in demo/psi_md_traj_2.pdb
    Read 1 chain(s) in demo/3ICH.pdb
    Read 90 chain(s) in demo/psi_md_traj_3.pdb
    Read 190 chain(s) in demo/beta3_IEGF12.pdb
    Read 1 chain(s) in demo/1BTA.pdb
    wrote test2.PB.fasta


`-x` and `-g` options
`````````````````````

.. warning:: These options require the installation of python library `MDAnalysis <http://www.mdanalysis.org/>`_

Instead using the `-p` option, the protein structures could come
from a molecular dynamics simulation file from Gromacs.
For this, you have to specify a '.xtc' file with the `-x` option and a '.gro' file with the `-g` option.

.. code-block:: bash

    $ PBassign -x `PBdata`/md_traj_4.xtc -g `PBdata`/md_traj_4.gro -o md_traj_4
    PBs assigned for demo/md.xtc | frame 1
    PBs assigned for demo/md.xtc | frame 2
    PBs assigned for demo/md.xtc | frame 3
    PBs assigned for demo/md.xtc | frame 4
    ...
    PBs assigned for demo/md.xtc | frame 198
    PBs assigned for demo/md.xtc | frame 199
    PBs assigned for demo/md.xtc | frame 200
    PBs assigned for demo/md.xtc | frame 201
    wrote md_traj_4.PB.fasta


`--phipsi` option
`````````````````

generates an additionnal file with the
`phi and psi angles <http://en.wikipedia.org/wiki/Dihedral_angle#Dihedral_angles_of_biological_molecules>`_
for each residue.

.. code-block:: bash

    $ PBassign -p `PBdata`/1BTA.pdb -o 1BTA --phipsi
    1 PDB file(s) to process
    Read 1 chain(s) in demo/1BTA.pdb
    wrote 1BTA.PB.fasta
    wrote 1BTA.PB.phipsi

Content of `1BTA.PB.phipsi`: ::

    demo/1BTA.pdb | chain A      1     None  -171.66
    demo/1BTA.pdb | chain A      2  -133.80   153.74
    demo/1BTA.pdb | chain A      3  -134.66   157.30
    demo/1BTA.pdb | chain A      4  -144.49   118.60
    demo/1BTA.pdb | chain A      5  -100.13    92.99
    demo/1BTA.pdb | chain A      6   -83.49   104.24
    demo/1BTA.pdb | chain A      7   -64.77   -43.25
    demo/1BTA.pdb | chain A      8   -44.48   -25.89
    demo/1BTA.pdb | chain A      9   -94.91   -47.18
    demo/1BTA.pdb | chain A     10   -41.31   133.74
    [snip]

The first part of the line is the comment also found in the fasta file.
The last thee columns are, from left to right, the residue number, the phi angle and the psi angle.
The phi angle of the first residue and the psi angle of the last residue cannot be computed.


`--flat` option
```````````````

formats the PBs assignment with one sequence per line.

.. code-block:: bash

    $ PBassign -p `PBdata`/1BTA.pdb -o 1BTA --flat
    1 PDB file(s) to process
    Read 1 chain(s) in demo/1BTA.pdb
    wrote 1BTA.PB.fasta
    wrote 1BTA.PB.flat

Content of `1BTA.PB.flat`: ::

    ZZdddfklonbfklmmmmmmmmnopafklnoiaklmmmmmnoopacddddddehkllmmmmngoilmmmmmmmmmmmmnopacdcddZZ


