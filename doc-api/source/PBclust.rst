PBclust
=======

.. warning:: NOT UPDATED.

Once converted to PB sequences, conformations of a same protein can be clustered
based on PB similarities.


Example
-------

.. code-block:: bash

    $ PBclust -f `PBdata`/psi_md_traj_all.PB.fasta -o psi_md_traj_all --clusters 5
    read 270 sequences in demo2/psi_md_traj_all.PB.fasta
    read substitution matrix
    Building distance matrix
    100%
    wrote psi_md_traj_all.PB.dist
    R clustering: OK
    cluster 1: 90 sequences (33%)
    cluster 2: 55 sequences (20%)
    cluster 3: 35 sequences (13%)
    cluster 4: 35 sequences (13%)
    cluster 5: 55 sequences (20%)
    wrote psi_md_traj_all.PB.clust


Cluster 1 is the biggest cluster with 33% of all conformations.
`psi_md_traj_all.PB.dist` contains the matrix distance between all PB sequences.

Content of `psi_md_traj_all.PB.clust` (clustering results): ::

    SEQ_CLU  "psi_md_traj_1.pdb | model 0"  1
    SEQ_CLU  "psi_md_traj_1.pdb | model 1"  1
    SEQ_CLU  "psi_md_traj_1.pdb | model 2"  1
    [snip]
    ...
    [snip]
    SEQ_CLU  "psi_md_traj_3.pdb | model 31"  4
    SEQ_CLU  "psi_md_traj_3.pdb | model 32"  4
    SEQ_CLU  "psi_md_traj_3.pdb | model 33"  5
    SEQ_CLU  "psi_md_traj_3.pdb | model 34"  5
    [snip]
    ...
    [snip]
    SEQ_CLU  "psi_md_traj_3.pdb | model 88"  5
    SEQ_CLU  "psi_md_traj_3.pdb | model 89"  5
    MED_CLU  "psi_md_traj_1.pdb | model 65"  1
    MED_CLU  "psi_md_traj_2.pdb | model 33"  2
    MED_CLU  "psi_md_traj_2.pdb | model 74"  3
    MED_CLU  "psi_md_traj_3.pdb | model 0"  4
    MED_CLU  "psi_md_traj_3.pdb | model 87"  5


Usage
-----

Here’s the ``PBclust`` help text. ::

    usage: PBclust [-h] -f F -o O (--clusters CLUSTERS | --compare)

    Cluster protein structures based on their PB sequences.

    optional arguments:
      -h, --help           show this help message and exit
      -f F                 name(s) of the PBs file (in fasta format)
      -o O                 name for results
      --clusters CLUSTERS  number of wanted clusters
      --compare            compare the first sequence versus all others


`--compare` option
``````````````````

compares, position by position, the first sequence found in the fasta file against all others.
The result of the comparison is a score between O (identical) and 9 (different).

.. code-block:: bash

    $ PBclust -f `PBdata`/psi_md_traj_all.PB.fasta -o psi_md_traj_all --compare
    read 270 sequences in demo2/psi_md_traj_all.PB.fasta
    read substitution matrix
    Normalized substitution matrix (between 0 and 9)
    [[0 3 2 3 4 3 3 4 2 3 5 3 5 4 3 3]
     [3 0 3 3 3 4 3 2 2 3 3 2 5 3 3 2]
     [2 3 0 3 4 3 2 4 3 4 5 5 5 4 3 2]
     [3 3 3 0 2 3 4 4 3 3 5 5 9 6 5 4]
     [4 3 4 2 0 2 2 2 4 3 3 4 7 4 5 5]
     [3 4 3 3 2 0 3 3 4 2 3 3 5 5 4 5]
     [3 3 2 4 2 3 0 3 3 3 4 3 3 2 2 1]
     [4 2 4 4 2 3 3 0 3 1 2 3 5 4 2 4]
     [2 2 3 3 4 4 3 3 0 2 2 2 5 3 3 2]
     [3 3 4 3 3 2 3 1 2 0 2 2 4 4 3 3]
     [5 3 5 5 3 3 4 2 2 2 0 3 3 3 4 4]
     [3 2 5 5 4 3 3 3 2 2 3 0 3 2 2 4]
     [5 5 5 9 7 5 3 5 5 4 3 3 0 2 3 3]
     [4 3 4 6 4 5 2 4 3 4 3 2 2 0 2 2]
     [3 3 3 5 5 4 2 2 3 3 4 2 3 2 0 2]
     [3 2 2 4 5 5 1 4 2 3 4 4 3 2 2 0]]
    Compare first sequence (psi_md_traj_1.pdb | model 0) with others
    wrote psi_md_traj_all.PB.compare.fasta

Content of `psi_md_traj_all.PB.compare.fasta`: ::

    >psi_md_traj_1.pdb | model 0 vs psi_md_traj_1.pdb | model 1
    00000002000000000020000000000002000200000000000230002000
    >psi_md_traj_1.pdb | model 0 vs psi_md_traj_1.pdb | model 2
    00000002000000000005000000000002000243000000055230000000
    >psi_md_traj_1.pdb | model 0 vs psi_md_traj_1.pdb | model 3
    00000002000000000020000000000002000200000000055230002000
    [snip]
    ...
    [snip]
    >psi_md_traj_1.pdb | model 0 vs psi_md_traj_3.pdb | model 87
    00302523340000000005000000035032000323300000335220000000
    >psi_md_traj_1.pdb | model 0 vs psi_md_traj_3.pdb | model 88
    00302523350500000005000000032232000323300000555225000000
    >psi_md_traj_1.pdb | model 0 vs psi_md_traj_3.pdb | model 89
    00333522250000000025000000035032000323300002035020002000
