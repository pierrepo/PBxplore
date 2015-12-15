Installation
============

Although it's not possible to cover all possible ways to install PBxplore on any operating system,
we try in this document to provide few guidelines regarding PBxplore setup.

Supported Platforms
-------------------

Currently, `PBxplore` run with Python 2.7, 3.3 and 3.4 on Linux and Mac OS X.


Dependencies
------------

To use `PBxplore`, the following libraries have to be installed.

    `NumPy <http://numpy.scipy.org/>`_ >= 1.6.0
        Numpy is the base package for numerical computing in python.

Optionally, `PBxplore` can use the following packages:

    `MDAnalysis <http://www.mdanalysis.org/>`_ [#]_ >= 0.11
        We use MDAnalysis for loading trajectories.
        See the full supported list
        `here <https://pythonhosted.org/MDAnalysis/documentation_pages/coordinates/init.html#id1>`_.

    `Matplotlib <http://matplotlib.org/>`_ >= 1.4.0
        All ploting functions use `matplotlib` package.

    `Weblogo3 <http://weblogo.threeplusone.com/>`_ [#]_
        `Weblogo3` is required to create logo from PB sequences.


Installing PBxplore
-------------------

Once dependencies installed, the most straightforward way is to use `pip`:

.. code-block:: bash

    $ pip install pbxplore


PBxplore can also be installed for the current user only:

.. code-block:: bash

    $ pip install --user pbxplore


See the documentation of `pip <https://pip.pypa.io/en/stable/>`_ for more information.
You may also want to look at `virtualenv <https://virtualenv.readthedocs.org/en/latest/>`_.


.. [#] N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
       MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
       *J. Comput. Chem.* **32** (2011), 2319–2327. doi:10.1002/jcc.21787

.. [#] G. E. Crooks, G. Hon, J.-M. Chandonia, and S. E. Brenner.
       WebLogo: A Sequence Logo Generator.
       *Genome Research* **14**: 1188–90 (2004) doi:10.1101/gr.849004
