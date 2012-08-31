"""
This file contains parameters for Mindboggle to fill settings in setup.py,
the Mindboggle top-level docstring, and for building the docs.
In setup.py we execute this file, so it cannot import mindboggle.
"""

# Mindboggle version information.  An empty _version_extra corresponds to a
# full release.  '.dev' as a _version_extra string means a development version
_version_major = 0
_version_minor = 1
_version_micro = 0
_version_extra = 'dev'
#_version_extra = ''

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
__version__ = "%s.%s.%s%s" % (_version_major,
                              _version_minor,
                              _version_micro,
                              _version_extra)

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: Apache v2.0",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

description  = "Automated human brain image anatomical labeling and shape analysis"

# Note: this long_description is actually a copy/paste from the top-level
# README.rst, so that it shows up nicely on PyPI.  So please remember to edit
# it only in one place and sync it correctly.
long_description = """
==========
Mindboggle
==========

Mindboggle is a package for automated anatomical labeling and morphometry
of human brain images.

Website
=======

Current information can always be found at the Mindboggle website::

    http://mindboggle.info

Code
====

You can find our sources and single-click downloads:

* `Main repository`_ on Github.
* Documentation_ for all releases and current development tree.
* Download as a tar/zip file the `current trunk`_.
* Downloads of all `available releases`_.

.. _main repository: http://github.com/binarybottle/mindboggle
.. _Documentation: http://mindboggle.info
.. _available releases: http://github.com/binarybottle/mindboggle/downloads

License
=======

Mindboggle is licensed under the terms of the Apache v2.0 license.

"""

# versions for dependencies
NUMPY_MIN_VERSION='1.2'

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
PLATFORMS           = "OS Independent"
MAJOR               = _version_major
MINOR               = _version_minor
MICRO               = _version_micro
ISRELEASE           = _version_extra == ''
VERSION             = __version__
PROVIDES            = ["mindboggle"]
REQUIRES            = ["numpy (>=%s)" % NUMPY_MIN_VERSION]

