PBassign
==========

``PBassign`` assigns a PB sequence to a protein structure.


Example
-------


First download data:

.. code-block:: bash

    $ wget https://files.rcsb.org/view/3ICH.pdb

Then perform PBs assignment:

.. code-block:: bash

    $ PBassign -p 3ICH.pdb -o 3ICH
    1 PDB file(s) to process
    Read 1 chain(s) in 3ICH.pdb
    wrote 3ICH.PB.fasta

Content of `3ICH.PB.fasta` : ::

    >3ICH.pdb | chain A
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

    optional arguments:
      -h, --help     show this help message and exit
      -p P           name of a pdb file or name of a directory containing pdb
                     files
      -o O           name for results
      -v, --version  show program's version number and exit

    other options to handle molecular dynamics trajectories:
      -x X           name of the topology file
      -g G           name of the trajectory file


``-p`` option
`````````````

can be used several times. For instance:

.. code-block:: bash

    $ wget https://files.rcsb.org/view/3ICH.pdb
    $ wget https://files.rcsb.org/view/1BTA.pdb
    $ wget https://files.rcsb.org/view/1AY7.pdb
    $ PBassign -p 3ICH.pdb -p 1BTA.pdb -p 1AY7.pdb -o test1
    3 PDB file(s) to process
    Read 1 chain(s) in 3ICH.pdb
    Read 1 chain(s) in 1BTA.pdb
    Read 2 chain(s) in 1AY7.pdb
    wrote test1.PB.fasta


All PB assignments are written in the same output file. If a PDB file contains several chains
and/or models, PBs assignments are also written in a single output file.
From the previous example, the ouput of ``test1.PB.fasta`` is: ::

    >3ICH.pdb | chain A
    ZZccdfbdcdddddehjbdebjcdddddfklmmmlmmmmmmmmnopnopajeopacfbdc
    ehibacehiamnonopgocdfkbjbdcdfblmbccfbghiacdddebehiafkbccddfb
    dcfklgokaccfbdcfbhklmmmmmmmpccdfkopafbacddfbgcddddfbacddddZZ
    >1BTA.pdb | chain A
    ZZdddfklonbfklmmmmmmmmnopafklnoiaklmmmmmnoopacddddddehkllmmm
    mngoilmmmmmmmmmmmmnopacdcddZZ
    >1AY7.pdb | chain A
    ZZbjadfklmcfklmmmmmmmmnnpaafbfkgopacehlnomaccddehjaccdddddeh
    klpnbjadcdddfbehiacddfegolaccdddfkZZ
    >1AY7.pdb | chain B
    ZZcddfklpcbfklmmmmmmmmnopafklgoiaklmmmmmmmmpacddddddehkllmmm
    mnnommmmmmmmmmmmmmnopacddddZZ


One can also use the ``-p`` option to provide a directory containing PDB files as an input.
``PBassign`` will process all PDB files located in the `PBdata` directory:

.. code-block:: bash

    $ wget https://files.rcsb.org/view/1AY7.pdb -P demo
    $ wget https://files.rcsb.org/view/2LFU.pdb -P demo
    $ wget https://files.rcsb.org/view/3ICH.pdb -P demo
    $ wget https://files.rcsb.org/view/1BTA.pdb -P demo
    $ PBassign -p demo/ -o test2
    4 PDB file(s) to process
    Read 1 chain(s) in demo/3ICH.pdb
    Read 2 chain(s) in demo/1AY7.pdb
    Read 1 chain(s) in demo/1BTA.pdb
    Read 10 chain(s) in demo/2LFU.pdb
    wrote test2.PB.fasta


``-x`` and ``-g`` options
`````````````````````````

.. warning:: These options use the `MDAnalysis <http://www.mdanalysis.org/>`_ library which is installed by PBxplore.

Instead using the ``-p`` option, protein structures could come from a molecular dynamics simulation trajectory file.
For this, you have to specify a trajectory file with the ``-x`` option and a topology file with the ``-g`` option.
It will accept any trajectory file format handled by the MDAnalysis library. See their [table of supported formats](https://pythonhosted.org/MDAnalysis/documentation_pages/coordinates/init.html#id1) for the full list.
Here an example with GROMACS files.

.. code-block:: bash

    $ wget https://raw.githubusercontent.com/pierrepo/PBxplore/master/demo_doc/psi_md_traj.gro
    $ wget https://raw.githubusercontent.com/pierrepo/PBxplore/master/demo_doc/psi_md_traj.xtc
    $ PBassign -x psi_md_traj.xtc -g psi_md_traj.gro -o psi_md_traj
    Frame 1/225.
    Frame 100/225.
    Frame 200/225.
    Frame 225/225.
    wrote psi_md_traj.PB.fasta

If needed, you can download ``psi_md_traj.PB.fasta`` [here](https://raw.githubusercontent.com/pierrepo/PBxplore/master/demo_doc/psi_md_traj.PB.fasta).


Tips'n tricks
-------------

To flatten the PB sequences obtained in FASTA format, i.e. get PB sequences in a single line each, one solution could be:

.. code-block:: bash

    $ wget https://files.rcsb.org/view/1AY7.pdb
    $ PBassign -p 1AY7.pdb -o 1AY7
    $ cat 1AY7.PB.fasta | sed "s/^>.*/\t/" | tr -d "\n" | tr "\t" "\n" > 1AY7.PB.flat

Content of `1AY7.PB.flat` : ::

    ZZbjadfklmcfklmmmmmmmmnnpaafbfkgopacehlnomaccddehjaccdddddehklpnbjadcdddfbehiacddfegolaccdddfkZZ
    ZZcddfklpcbfklmmmmmmmmnopafklgoiaklmmmmmmmmpacddddddehkllmmmmnnommmmmmmmmmmmmmnopacddddZZ

