.. _install:

======================
 Download and install
======================

This page covers the necessary steps to install Mindboggle.

Download
--------

http://mindboggle.info/download

To check out the latest development version::

        git clone git://github.com/binarybottle/mindboggle.git

Install
-------

Step 1: Install Mindboggle Python software
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can skip this step if you checked out the development version.
If you downloaded the source distribution named something like
``mindboggle-x.y.tar.gz``, then unpack the tarball, change into the
``mindboggle-x.y`` directory and install mindboggle using::

    python setup.py install

**Note:** Depending on permissions you may need to use ``sudo``.

Step 2: Compile C++ code
~~~~~~~~~~~~~~~~~~~~~~~~

Move the mindboggle_cpp directory to wherever you install your software
(let's call it YOUR_SOFTWARE_PATH), then run cmake and make::

>> mv mindboggle_cpp YOUR_SOFTWARE_PATH
>> cd YOUR_SOFTWARE_PATH/mindboggle_cpp/bin/
>> cmake ..  [OR "ccmake .." and select appropriate architecture]
>> make

Step 3: Set environment variable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For your system to recognize the C++ code we just compiled,
include the following lines in your system's equivalent to a .bashrc file::

MINDBOGGLE_CPP=YOUR_SOFTWARE_PATH/mindboggle_cpp/bin
PATH=${MINDBOGGLE_CPP}:${PATH}
export MINDBOGGLE_CPP

Step 4: Move directory
~~~~~~~~~~~~~~~~~~~~~~

To take advantage of the FreeSurfer_ software, move the contents of the
data/ directory to FreeSurfer's subjects/ directory, and make sure that
the environment variable SUBJECTS_DIR is set in your system's equivalent
to a .bashrc file::

>> mv data/* $FREESURFER_HOME/subjects/

Dependencies
------------

Below is a list of required dependencies, along with additional software
recommendations.

Must Have
~~~~~~~~~

Basic compilation software::
  - easy_install
    ????sudo easy_install --always-unzip --upgrade setuptools
  - cmake (http://cmake.org)
  - C++ compiler such as gcc

`Visualization Toolkit <http://vtk.org>`_

Python_ 2.6 -2.7

Nipype_ 1.0
  Neuroimaging software pipeline framework

Nibabel_ 1.0
  Neuroimaging file i/o library

NumPy_ 1.3 - 1.6

SciPy_ 0.7 - 0.10
  Numpy and Scipy are high-level, optimized scientific computing libraries

NetworkX_ 1.0 - 1.4
  Python package for working with complex networks

.. note::

    Full distributions such as pythonxy_ or EPD_ provide the above packages,
    except Nibabel_ and Nipype_.

Strong Recommendations
~~~~~~~~~~~~~~~~~~~~~~

FreeSurfer_
  FreeSurfer version 4 and higher

IPython_ 0.10.2 - 0.12
  Interactive python environment. This is necessary for some parallel
  components of the pipeline engine.

Matplotlib_ 1.0
  Plotting library

Sphinx_
  Required for building the documentation

`Graphviz <http://www.graphviz.org/>`_
  Required for building the documentation
