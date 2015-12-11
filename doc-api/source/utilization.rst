Utilization
===========

2 ways to use `PBxplore`:
  - command-line scripts
  - API

Command-line scripts
--------------------

PBxplore is bundle with a suite of command-line tools.
Once installed, they should be available in your ``$PATH``.
Here the list:

- :doc:`PBassign <PBassign>` assign :doc:`Protein Blocks (PBs) <intro_PB>`
  from a set of protein structures.
- PBcount computes the frequency of PBs at each position along the amino acid sequence.
- PBstast generates frequency and logo plots, and estimates the (_Neq_).
- PBclust use clustering algorithm (k-means) to re-group similar PBs sequences.


API
---



.. _demo:

Demo files
----------

PBxplore provides scripts to demonstrate its functionalities. These scripts
guide the user through the different command line tools of PBxplore.

3 demonstration scripts are available:

* `run_demo1_assignation.sh` demonstrates the how to use ``PBassign`` to assign
  Protein Block sequences to protein structures;
* `run_demo2_statistics.sh` demonstrates how to analyse protein dynamics using
  PBxplore;
* finally, `run_demo2_clusters.sh` demonstrates how to use ``PBclust`` to cluster
  Protein Block sequences.

In addition to the scripts, PBxplore bundles input files to test and get
accustom with the software and the python library. Call the ``PBdata`` program to
display where these files are stored. With the `--list` and `--list-abs`
arguments ``PBdata`` list the available name and absolute path, respectively.
