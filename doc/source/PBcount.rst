PBcount
=======

``PBcount`` computes the frequency of PBs at each position along the amino acid sequence.

.. note:: The following examples use ``PBdata`` and the demo files.
          See :ref:`Demo files <demo>` for more information.

Example
-------

.. code-block:: bash

    $ PBcount -f `PBdata`/psi_md_traj_1.PB.fasta -o psi_md_traj_1
    read 90 sequences in demo/psi_md_traj_1.PB.fasta
    wrote psi_md_traj_1.PB.count

Content of `psi_md_traj_1.PB.count`: ::

             a     b     c     d     e     f     g     h     i     j     k     l     m     n     o     p
    1        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
    2        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
    3        0     0     0     0     0    90     0     0     0     0     0     0     0     0     0     0
    4        0     0     0     0     0     1     0     0     0     0    89     0     0     0     0     0
    [snip]
    51       0     0     0     0     0    22     0    40     0     0    28     0     0     0     0     0
    52       0    23     0     0     0     0     0     0    38     1     1    27     0     0     0     0
    53      62     0    21     0     0     0     0     0     0     0     0     0     0     0     0     7
    54       0     0    90     0     0     0     0     0     0     0     0     0     0     0     0     0
    55       0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
    56       0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0

Note that residues 1, 2, 55 and 56 have a null count of all PBs.
These residues are the first and last residues of the structure and no PB can be assigned to them.

Usage
-----

Here’s the ``PBcount`` help text. ::

    usage: PBcount [-h] -f F -o O [--first-residue FIRST_RESIDUE]

    Compute PB frequency along protein sequence.

    optional arguments:
      -h, --help            show this help message and exit
      -f F                  name(s) of the PBs file (in fasta format)
      -o O                  name for results
      --first-residue FIRST_RESIDUE
                            define first residue number (1 by default)


`-f` option
```````````

can be used several times:

.. code-block:: bash

    $ PBcount -f `PBdata`/psi_md_traj_1.PB.fasta -f `PBdata`/psi_md_traj_2.PB.fasta -f `PBdata`/psi_md_traj_3.PB.fasta -o psi_md_traj_all
    read 90 sequences in demo/psi_md_traj_1.PB.fasta
    read 90 sequences in demo/psi_md_traj_2.PB.fasta
    read 90 sequences in demo/psi_md_traj_3.PB.fasta
    wrote psi_md_traj_all.PB.count


`--first-residue` option
````````````````````````

By default, the number of the first residue is 1, this option allows
to adjust the number associated to the first residue (and to the followings automaticaly).

.. code-block:: bash

    $ PBcount --first-residue 5 -f `PBdata`/psi_md_traj_1.PB.fasta -o psi_md_traj_1_shifted
    read 90 sequences in demo/psi_md_traj_1.PB.fasta
    wrote psi_md_traj_1_shifted.PB.count


Content of `psi_md_traj_1_shifted.PB.count`: ::

             a     b     c     d     e     f     g     h     i     j     k     l     m     n     o     p
    5        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
    6        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
    7        0     0     0     0     0    90     0     0     0     0     0     0     0     0     0     0
    8        0     0     0     0     0     1     0     0     0     0    89     0     0     0     0     0
    9        0    89     0     0     0     0     0     0     0     0     0     1     0     0     0     0
    10       0     0    86     0     0     3     0     0     0     0     0     0     1     0     0     0
    [snip]
