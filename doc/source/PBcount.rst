PBcount
=======

``PBcount`` computes the frequency of PBs at each position along the amino acid sequence.

.. note:: 

    The following examples require ``psi_md_traj.PB.fasta`` and 
          ``psi_md_traj2.PB.fasta`` obtained with ``PBassign``:
          
    .. code-block:: bash

        $ wget https://raw.githubusercontent.com/pierrepo/PBxplore/master/demo_doc/psi_md_traj_1.gro
        $ wget https://raw.githubusercontent.com/pierrepo/PBxplore/master/demo_doc/psi_md_traj_1.xtc
        $ PBassign -x psi_md_traj_1.xtc -g psi_md_traj_1.gro -o psi_md_traj_1
        Frame 1/225.
        Frame 100/225.
        Frame 200/225.
        Frame 225/225.
        wrote psi_md_traj_1.PB.fasta

        $ wget https://raw.githubusercontent.com/pierrepo/PBxplore/master/demo_doc/psi_md_traj_2.gro
        $ wget https://raw.githubusercontent.com/pierrepo/PBxplore/master/demo_doc/psi_md_traj_2.xtc
        $ PBassign -x psi_md_traj_2.xtc -g psi_md_traj_2.gro -o psi_md_traj_2
        Frame 1/225.
        Frame 100/225.
        Frame 200/225.
        Frame 225/225.
        wrote psi_md_traj_2.PB.fasta

If needed, you can download [psi_md_traj_1.PB.fasta](https://raw.githubusercontent.com/pierrepo/PBxplore/master/demo_doc/psi_md_traj_1.PB.fasta) and [psi_md_traj_2.PB.fasta](https://raw.githubusercontent.com/pierrepo/PBxplore/master/demo_doc/psi_md_traj_2.PB.fasta).

          
Example
-------

.. code-block:: bash

    $ PBcount -f psi_md_traj_1.PB.fasta -o psi_md_traj_1
    read 225 sequences in psi_md_traj_1.PB.fasta
    wrote psi_md_traj_1.PB.count

Content of `psi_md_traj_1.PB.count`: ::

             a     b     c     d     e     f     g     h     i     j     k     l     m     n     o     p
    1        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
    2        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
    3        0     0     0     4     0   221     0     0     0     0     0     0     0     0     0     0
    4        0     0     0     0     0     5     0     0     0     0   220     0     0     0     0     0
    [snip]
    51       0     0     0     0     0    56     0    98     0     0    71     0     0     0     0     0
    52       0    56     0     0     0     0     0     0    94     3     3    69     0     0     0     0
    53     144     0    60     2     0     0     0     0     0     0     0     0     1     0     0    18
    54       0     0   225     0     0     0     0     0     0     0     0     0     0     0     0     0
    55       0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
    56       0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0

Note that residues 1, 2, 55 and 56 have a null count of all PBs.
These residues are the first and last residues of the structure and no PB can be assigned to them.

If needed, you can download [psi_md_traj_1.PB.count](https://raw.githubusercontent.com/pierrepo/PBxplore/master/demo_doc/psi_md_traj_1.PB.count).


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
      -v, --version         show program's version number and exit


`-f` option
```````````

can be used several times:

.. code-block:: bash

    $ PBcount -f psi_md_traj_1.PB.fasta -f psi_md_traj_2.PB.fasta -o psi_md_traj_all
    read 225 sequences in psi_md_traj_1.PB.fasta
    read 225 sequences in psi_md_traj_2.PB.fasta
    wrote psi_md_traj_all.PB.count

If needed, you can download [psi_md_traj_all.PB.count](https://raw.githubusercontent.com/pierrepo/PBxplore/master/demo_doc/psi_md_traj_all.PB.count).

`--first-residue` option
````````````````````````

By default, the number of the first residue is 1, this option allows
to adjust the number associated to the first residue (and to the followings automaticaly).

.. code-block:: bash

    $ PBcount --first-residue 5 -f psi_md_traj_1.PB.fasta -o psi_md_traj_1_shifted
    read 225 sequences in psi_md_traj_1.PB.fasta
    wrote psi_md_traj_1_shifted.PB.count


Content of `psi_md_traj_1_shifted.PB.count`: ::

             a     b     c     d     e     f     g     h     i     j     k     l     m     n     o     p
    5        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
    6        0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0
    7        0     0     0     4     0   221     0     0     0     0     0     0     0     0     0     0
    8        0     0     0     0     0     5     0     0     0     0   220     0     0     0     0     0
    9        0   222     0     0     0     0     0     0     0     0     0     3     0     0     0     0
    10       6     0   201     0     0     5     0     0     0     0     0     0    12     0     1     0
    [snip]

If needed, you can download [psi_md_traj_1_shifted.PB.count](https://raw.githubusercontent.com/pierrepo/PBxplore/master/demo_doc/psi_md_traj_1_shifted.PB.count).

