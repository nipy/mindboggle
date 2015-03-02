.. _howto_make_html:

==================
Mindboggle website
==================

Note::
 make html
 rsync -avz --sparse --exclude-from=/homedir/.rsync-exclude -e /usr/bin/ssh build/html/* binarybottle@binarybottle.com:/home/binarybottle/mindboggle.info

This is the top level build directory for the Mindboggle documentation.
The Mindboggle documentation is adapted from the documentation for nipy_,
and is written using Sphinx_, a python documentation
system built on top of reST_.  In order to build the documentation,
you must have Sphinx v0.5 or greater installed.

This directory contains:

* Makefile - the build script to build the HTML or PDF docs.
  Type ``make help`` for a list of options.

* links.txt - reST document with hyperlink targets for common
  links used throughout the documentation

* .rst files - some top-level documentation source files

* conf.py - the sphinx configuration

* _static - used by the sphinx build system.

* _templates - used by the sphinx build system


.. _Sphinx: http://sphinx.pocoo.org/
.. _reST: http://docutils.sourceforge.net/rst.html
.. _numpy: http://www.scipy.org/NumPy
