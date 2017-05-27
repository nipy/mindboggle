#!/usr/bin/env python
"""
This file contains parameters for Mindboggle to fill settings in setup.py,
the Mindboggle top-level docstring, and for building the docs.
In setup.py we execute this file, so it cannot import mindboggle.
"""

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
from mindboggle.version import __version__ as __version__

CLASSIFIERS = ["Development Status :: Beta",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: Apache v2.0",
               "Operating System :: Linux",
               "Programming Language :: Python 3",
               "Topic :: Scientific/Engineering"]

description  = "Automated human brain image feature extraction, labeling, and shape analysis"

# Note: this long_description is actually a copy/paste from the top-level
# README.rst, so that it shows up nicely on PyPI.  So please remember to edit
# it only in one place and sync it correctly.
long_description = """
==========
Mindboggle
==========

Mindboggle is a software package for automated feature extraction, anatomical
labeling, and morphometry of human brain magnetic resonance images,
licensed under the terms of the Apache v2.0 license.
Current information can always be found at the Mindboggle website,
http://mindboggle.info, and on the Github main repository,
http://github.com/nipy/mindboggle

"""

# Main setup parameters
NAME                = 'Mindboggle'
MAINTAINER          = "Arno Klein"
MAINTAINER_EMAIL    = "arno@mindboggle.info"
DESCRIPTION         = description
LONG_DESCRIPTION    = long_description
URL                 = "http://mindboggle.info/"
DOWNLOAD_URL        = "http://mindboggle.info/"
LICENSE             = "Apache v2.0"
CLASSIFIERS         = CLASSIFIERS
AUTHOR              = "Arno Klein"
AUTHOR_EMAIL        = "arno@mindboggle.info"
PLATFORMS           = "Linux"
VERSION             = __version__
PROVIDES            = ["mindboggle"]
#REQUIRES            = ["numpy (>={0})".format(NUMPY_MIN_VERSION)]

