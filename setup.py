#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import
from setuptools import setup, find_packages
import os

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'README.rst')) as f:
    readme = f.read()


# Extras requirements for optional dependencies
extras = {
    'analysis': ['weblogo', 'matplotlib'],
    'trajectories': ['MDAnalysis>=0.11'],
    'all': ['weblogo', 'matplotlib', 'MDAnalysis>=0.11']
}

# Version number must be in sync with the one in pbxplore/__init__.py
setup(
    name='pbxplore',
    version='1.3.1',

    description="PBxplore is a suite of tools dedicated to Protein Block analysis.",
    long_description=readme,

    url='https://github.com/pierrepo/PBxplore',

    # Author details
    author='Pierre Poulain',
    author_email='pierre.poulain@cupnet.net',

    license='MIT',

    classifiers=[
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],

    install_requires=['numpy'],
    tests_require=['nose', 'coverage'],
    # List additional groups of dependencies here
    # To install, use
    # $ pip install -e .[analysis]
    extras_require=extras,

    packages=find_packages(exclude=['test']),
    include_package_data=True,
    package_data={'pbxplore':['demo/*']},

    entry_points={
        'console_scripts': [
            'PBassign = pbxplore.scripts.PBassign:pbassign_cli',
            'PBclust  = pbxplore.scripts.PBclust:pbclust_cli',
            'PBcount  = pbxplore.scripts.PBcount:pbcount_cli',
            'PBstat   = pbxplore.scripts.PBstat:pbstat_cli',
            'PBdata   = pbxplore.scripts.PBdata:pbdata_cli',
        ],
    },

)
