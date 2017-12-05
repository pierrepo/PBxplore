.. PBxplore API documentation master file, created by
   sphinx-quickstart on Tue Oct  6 15:44:18 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PBxplore's documentation!
====================================

**PBxplore is a suite of tools dedicated to Protein Block analysis.**

Protein Blocks (PBs) are structural prototypes defined by
`de Brevern <http://www.dsimb.inserm.fr/~debrevern/index.php>`_ *et al* [#]_.
The 3-dimensional local structure of a protein backbone can be modelized as an 1-dimensional
sequence of PBs.
In principle, any conformation of any amino acid could be represented
by one of the sixteen available Protein Blocks.

.. figure:: img/PBs.jpg
    :align: center

    Schematic representation of the sixteen protein blocks, labeled from *a* to *p*
    (`Creative commons 4.0 CC-BY <https://creativecommons.org/licenses/by/4.0/>`_).




PBxplore provides both a Python library and command-line tools (:doc:`Utilization <utilization>`).
For some features (reading trajectory, plots), PBxplore requires
some optional dependencies (:doc:`Installation <installation>`).
Basically, PBxplore can:

* **assign PBs** from either a PDB or either a molecular dynamics trajectory (:doc:`PB assignation <./notebooks/Assignement>`).
* use analysis tools to perform **statistical analysis** on PBs (:doc:`Statistics <PBstat>`).
* use analysis tools to **study protein flexibility and deformability** (:doc:`Analyzis <./notebooks/Deformability>`).


.. raw:: html

   <div style="display:none">

.. toctree::
   :maxdepth: 1

   intro_PB
   installation
   utilization
   PBassign
   PBcount
   PBstat
   api_cookbook
   api_reference

.. raw:: html

   </div>


Citation
--------

If you use PBxplore, please cite this tool as:

| Barnoud J, Santuz H, Craveur P, Joseph AP, Jallu V, de Brevern AG, Poulain P,
| PBxplore: a tool to analyze local protein structure and deformability with Protein Blocks
| *PeerJ*  5:e4013 `<https://doi.org/10.7717/peerj.4013>`_ (2017).


Contact & Support
-----------------

PBxplore is a research software and has been developped by:

* Pierre Poulain, "Mitochondria, Metals and Oxidative Stress" group, Institut Jacques Monod, UMR 7592, Univ. Paris Diderot, CNRS, Sorbonne Paris Cité, France.
* Jonathan Barnoud, University of Groningen, Groningen, The Netherlands.
* Hubert Santuz, Laboratoire de Biochimie Théorique, CNRS UPR 9080, Institut de Biologie Physico-Chimique, Paris, France.
* Alexandre G. de Brevern, DSIMB, INSERM, U 1134, INTS, Paris, France.

If you want to report a bug or request a feature,
please use the `GitHub issue system <https://github.com/pierrepo/PBxplore/issues>`_.


Licence
-------

PBxplore is licensed under `The MIT License <https://github.com/pierrepo/PBxplore/blob/master/LICENSE>`_.


.. [#] A. G. de Brevern, C. Etchebest, and S. Hazout. Bayesian Probabilistic Approach for Predicting Backbone Structures in Terms of Protein Blocks. *Proteins* **41**:271-87 (2000).
