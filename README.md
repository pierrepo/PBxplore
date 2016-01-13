# PBxplore [![Build Status](https://travis-ci.org/pierrepo/PBxplore.svg?branch=master)](https://travis-ci.org/pierrepo/PBxplore)

---

**PBxplore** is a suite of tools dedicated to Protein Block analysis.
Protein Blocks are structural prototypes defined by
[de Brevern](http://www.dsimb.inserm.fr/~debrevern/index.php) *et al* [1]. The 3-dimensional local
structure of a protein backbone can be modelized as an 1-dimensional sequence of PBs.
In principle, any conformation of any amino acid could be represented by one of
the sixteen available Protein Blocks (see Figure 1).

![PBs](doc/source/img/PBs.jpg "PBs")

**Figure 1.** Schematic representation of the sixteen protein blocks,
              labeled from *a* to *p* ([Creative commons 4.0 CC-BY](https://creativecommons.org/licenses/by/4.0/)).


PBxplore provides both a python library and command-line tools. Basically, PBxplore can:

* **assign PBs** from either a PDB or either a molecular dynamics trajectory.
* use analysis tools to perform **statistical analysis** on PBs.
* use analysis tools to **study protein flexibility and deformability**.

For details, see the documentation at [https://pbxplore.readthedocs.org/en/latest/](https://pbxplore.readthedocs.org/en/latest/).

## Download

- [Get latest zip archive](https://github.com/pierrepo/PBxplore/archive/master.zip)
- Clone repository: `git clone git@github.com:pierrepo/PBxplore.git` or `git clone https://github.com/pierrepo/PBxplore.git`

## Requirements

PBxplore requires:

* Python 2.7 or Python 3.x (>= 3.3)
* the [NumPy](http://numpy.scipy.org/ "NumPy") Python library,

Optionally, PBxplore can use:

* the [MDAnalysis](https://code.google.com/p/mdanalysis/) Python library (version >= 0.11) to read MD trajectories generated by Gromacs (.xtc files),
* [Mapltolib](http://matplotlib.org/) to generate plots.
* [WebLogo 3](http://weblogo.threeplusone.com/) to create logo from PB sequences.

## Documentation

All documentation are hosted by Read The Docs and can be found [here](https://pbxplore.readthedocs.org/en/latest/).

## Contact & Support

PBxplore is a research software and has been developped by:

* Pierre Poulain, DSIMB, Ets Poulain, Pointe-Noire, Congo
* Jonathan Barnoud, University of Groningen, Groningen, The Netherlands
* Hubert Santuz, DSIMB, Paris, France
* Alexandre G. de Brevern, DSIMB, Paris, France

If you want to report a bug, request a feature,
use the [GitHub issue system](https://github.com/pierrepo/PBxplore/issues>).


## Licence

PBxplore is licensed under [The MIT License](https://github.com/pierrepo/PBxplore/blob/master/LICENSE).



## References
[1] A. G. de Brevern, C. Etchebest and S. Hazout. Bayesian probabilistic approach for predicting backbone structures in terms of protein blocks. *Proteins* **41**: 271-288 (2000).
