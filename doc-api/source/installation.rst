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

    `MDAnalysis <http://www.mdanalysis.org/>`_ >= 0.11
        We use MDAnalysis for loading trajectories.
        See the full supported list
        `here <https://pythonhosted.org/MDAnalysis/documentation_pages/coordinates/init.html#id1>`_.

    `Matplotlib <http://matplotlib.org/>`_ >= 1.4.0
        All ploting functions use `matplotlib` package.

    `Weblogo3 <http://weblogo.threeplusone.com/>`_
        `Weblogo3` is required to create logo from PB sequences.


Installing PBxplore
-------------------

Once dependencies installed, the most straightforward way is to use `pip`:

.. code-block:: bash

    $ pip install pbxplore
