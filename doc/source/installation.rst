Installation
============

Although it's not possible to cover all possible ways to install PBxplore on any operating system,
we try in this document to provide few guidelines regarding PBxplore setup.

Supported Platforms
-------------------

Currently, `PBxplore` run with Python 2.7, 3.4, 3.5 and 3.6 on Linux and Mac OS X.


Dependencies
------------

`PBxplore` requires the following libraries (which will be installed through ` pip`).

    `NumPy <http://numpy.scipy.org/>`_ >= 1.6.0
        NumPy is the base package for numerical computing in python.

    `matplotlib <http://matplotlib.org/>`_ [#]_ >= 1.4.0
        All ploting functions use the `matplotlib` package.

    `MDAnalysis <http://www.mdanalysis.org/>`_ [#]_ >= 0.11
        We use MDAnalysis to read Gromacs molecular dynamics trajectories.
        Many other trajectory files are also supported, see list
        `here <https://pythonhosted.org/MDAnalysis/documentation_pages/coordinates/init.html#id1>`_.


Optionally, `PBxplore` can use the following packages:

    `Weblogo3 <http://weblogo.threeplusone.com/>`_ [#]_
        `Weblogo3` is required to create logo from PB sequences. It has to be installed by the user.


Installing PBxplore
-------------------

The most straightforward way is to use `pip`. It will also install the required dependencies:

.. code-block:: bash

    $ pip install --user pbxplore


If none of the dependencies were previously installed, the installation process can take up to several minutes.

.. note::

    The former command will install the PBxplore command-line scripts in:

    - ``$HOME/.local/bin`` on Linux.
    - ``~/Library/Python/X.Y/bin`` on Mac OSX, where ``X.Y`` stands for the Python version (for instance ``2.7`` for Python 2.7 and ``3.5`` for Python 3.5).

    Do not forget to update your ``$PATH`` environment variable to make all PBxplore tools accessible. As an example, with the Bash shell, this gives:

    .. code-block:: bash

        $ # for Linux
        $ export PATH=$PATH:$HOME/.local/bin

        $ # for Mac OSX
        $ export PATH=$PATH:~/Library/Python/X.Y/bin


See the documentation of `pip <https://pip.pypa.io/en/stable/>`_ for more information.
You may also want to look at `virtualenv <https://virtualenv.readthedocs.org/en/latest/>`_.



Upgrading PBxplore
---------------------

Be sure to always have the latest version of PBxplore with:

.. code-block:: bash

    $ pip install --user --upgrade pbxplore


Testing PBxplore
----------------

`PBxplore` comes with unit tests and regression tests. It requires the package
`pytest <https://docs.pytest.org>`_. You can run the tests within Python:

.. code-block:: python

    import pbxplore
    pbxplore.test()


Uninstalling PBxplore
---------------------

Run the simple command:

.. code-block:: bash

    $ pip uninstall pbxplore



PBxplore for advanced users
---------------------------

You can clone PBxplore from GitHub:

.. code-block:: bash

    $ git clone --depth 1 https://github.com/pierrepo/PBxplore.git

Once in the ``PBxplore`` directory, we advise you to create a virtual environment:

.. code-block:: bash

    $ pip3 install --user virtualenv
    $ virtualenv -p python3 pbxplore-py3
    $ source venv/bin/activate

You can then install the latest version of PBxplore as a Python module:

.. code-block:: bash

    $ pip install -e .

You can also run unit tests and regression tests:

.. code-block:: bash

    $ pip install pytest
    $ pytest -v pbxplore/tests

or

.. code-block:: bash

    $ pip install pytest
    $ python setup.py test


.. [#] J. D. Hunter.
       Matplotlib: A 2D graphics environment.
       *Computing In Science and Engineering* **9** (2007), 90-95. doi:10.1109/MCSE.2007.55

.. [#] N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
       MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
       *J. Comput. Chem.* **32** (2011), 2319–2327. doi:10.1002/jcc.21787

.. [#] G. E. Crooks, G. Hon, J.-M. Chandonia, and S. E. Brenner.
       WebLogo: A Sequence Logo Generator.
       *Genome Research* **14**: 1188–90 (2004) doi:10.1101/gr.849004
