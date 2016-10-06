## How to release

#### Update version number


#### Publish Documentation


#### Publish on Pypi


#### Publish on Conda

1. Build the package (inside `conda` directory)

    conda build --python 2.7 --python 3.3 --python 3.4 --python 3.5 -c bioconda -c mdanalysis pbxplore

2. The path of the created archive will be given at the end of the process

3. Convert it for all plateforms (where `X.X.X` is the pbxplore version and `YY` is the Python version)

    conda convert --platform all /path/to/archive/pbxplore-X.X.X-pyYY_0.tar.bz2

4. Upload it to anaconda

    anaconda upload /path/to/archive/pbxplore-X.X.X-pyYY_0.tar.bz2

This has to be done for each package created.
