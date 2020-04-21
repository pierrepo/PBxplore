## How to release

Before any release, double-check all tests had run successfully.

Then install needed tools:

    pip install --user bumpversion twine


#### Update version number

We use the tool [bumpversion](https://github.com/peritus/bumpversion) to synchronize the version number
across different files:

    bumpversion --verbose --config-file devtools/bumpversion.cfg  patch
    git push origin
    git push origin --tags

#### Publish on PyPI

Initial setup:

- Create an account on PyPI and create a projet (here [pbxplore](https://pypi.python.org/pypi/pbxplore)).

- Create a .pypirc files into your HOME directory with the following lines:

        [distutils]
        index-servers = pypi

        [pypi]
        repository=https://upload.pypi.org/legacy/
        username=your_username

- Prepare package files (at minimum `setup.py`, `setup.cfg`, `README.md`, `MANIFEST.in` and `LICENSE.txt`).

- Register package to PyPI :

        python setup.py register -r pypi

Build package:

    python setup.py sdist

Upload package to PyPI:

    twine upload dist/*

Enter your password when required.

Clean local package:

    rm -f dist/*.tar.gz

Doc:

- https://packaging.python.org/distributing/#uploading-your-project-to-pypi


### Add new release on GitHub

On [GitHub release page](https://github.com/pierrepo/PBxplore/releases) :

- Click the *Draft a release* button.
- Select the latest version as *tag version*.
- Add release version as *Release title* (e.g.: v1.3.7).
- Copy and paste the content of the CHANGELOG in the *Describe this release* field.
- Hit the *Publish Release* button.


### Publish in Zenodo 

- Login on Zenodo: <https://zenodo.org/login>
- Check new PBxplore version coming from GitHub: <https://zenodo.org/deposit?page=1&size=20>


#### Publish Documentation

Documentation is automatically retrieved and build by ReadTheDocs.

Just check revision at the page footer corresponds to the last commit: 
<https://pbxplore.readthedocs.io/en/latest/index.html>


#### Publish on Conda


1. Build the package (inside `conda` directory)

    `conda-build --python 2.7 --python 3.4 --python 3.5 --python 3.6 -c bioconda -c mdanalysis pbxplore`


2. The path of the created archive will be given at the end of the process


3. Convert it for all plateforms (where `X.X.X` is the pbxplore version and `YY` is the Python version)

    `conda-convert --platform all /path/to/archive/pbxplore-X.X.X-pyYY_0.tar.bz2`


4. Upload it to anaconda

    `anaconda upload /path/to/archive/pbxplore-X.X.X-pyYY_0.tar.bz2`


This has to be done for each package created.
