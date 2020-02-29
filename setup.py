#! /usr/bin/env python3
'''Just a little setup file with the basics, in case we ever need that.'''

import sys
from setuptools import setup

if sys.version_info[0] < 3:
    print('STOP USING PYTHON 2!i')
    sys.exit(-1)

with open("README.md", "r") as fh:
    long_description = fh.read()

print("Sifter unittests haven't been built into setup.py scripts. "
      "Please manually run '$ pytest tests/test_*.py' to ensure everything "
      "works before production use.")

dependencies = ['numpy >= 1.15.0',
                'json >= 2.0.9',
                'argparse >= 1.1',
                'sqlite >= 2.6.0',
                'pytest >= 4.0.1',
                'astropy >= 3.0.3',
                'astropy_healpix >= 0.3.1'
                ]

classifiers = ["Programming Languagge :: Python :: 3",
               "License :: OSI Approved :: MIT License",
               "Operating System :: OS Independent",
               "Development Status :: 2 - Pre-Alpha",
               "Intended Audience :: Science/Research",
               "Natural Language :: English",
               "Topic :: Scientific/Engineering :: Astronomy"
               ]

setup(name='sifter',
      version='0.1',
      author='Matt Payne, Mike Alexandersen',
      author_email='mpayne@cfa.harvard.edu, mike.alexandersen@alumni.ubc.ca',
      description='Sift through ITF/observations to find matches.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/matthewjohnpayne/sifter',
      download_url='https://github.com/fraserw/trippy/archive/master.zip',
      packages=['sifter'],
      classifiers=classifiers,
      keywords=['Sift', 'ITF', 'obs80'],
      license='MIT',
      install_requires=dependencies,
      )
