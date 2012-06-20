.. _install:

======================
 Download and install
======================

This page covers the necessary steps to install Mindboggle.

Download
--------

Release 0.5.3: [`zip <http://github.com/binarybottle/mindboggle/zipball/0.5.3>`_ `tar
<http://github.com/binarybottle/mindboggle/tarball/0.5.3>`_]

Development: [`zip <http://github.com/binarybottle/mindboggle/zipball/master>`_ `tar
<http://github.com/binarybottle/mindboggle/tarball/master>`_]

`Prior downloads <http://github.com/binarybottle/mindboggle/tags>`_

To check out the latest development version::

        git clone git://github.com/binarybottle/mindboggle.git

Install
-------

The installation process is similar to other Python packages.

If you already have a Python environment setup that has the dependencies listed
below, you can do::

	easy_install mindboggle

or::

	pip install mindboggle

Debian and Ubuntu
~~~~~~~~~~~~~~~~~

Add the `NeuroDebian <http://neuro.debian.org>`_ repository and install
the ``python-mindboggle`` package using ``apt-get`` or your favourite package
manager.

Mac OS X
~~~~~~~~

The easiest way to get mindboggle running on MacOSX is to install EPD_ and then add
nibabel and mindboggle by executing::

	easy_install nibabel
	easy_install mindboggle

If you are running a 64 bit version of EPD, you will need to compile
ETS. Instructions for a 64-bit boot mode are available:
https://gist.github.com/845545


From source
~~~~~~~~~~~

If you downloaded the source distribution named something
like ``mindboggle-x.y.tar.gz``, then unpack the tarball, change into the
``mindboggle-x.y`` directory and install mindboggle using::

    python setup.py install

**Note:** Depending on permissions you may need to use ``sudo``.


Testing the install
-------------------

The best way to test the install is to run the test suite.  If you have
nose_ installed, then do the following::

    python -c "import mindboggle; mindboggle.test()"

you can also test with nosetests::

    nosetests --with-doctest /software/binarybottle-repo/mastermindboggle/mindboggle
    --exclude=external --exclude=testing

All tests should pass (unless you're missing a dependency). If SUBJECTS_DIR
variable is not set some FreeSurfer related tests will fail. If any tests
fail, please report them on our `bug tracker
<http://github.com/binarybottle/mindboggle/issues>`_.


Dependencies
------------

Below is a list of required dependencies, along with additional software
recommendations.

Must Have
~~~~~~~~~

Python_ 2.6 -2.7

Nibabel_ 1.0
  Neuroimaging file i/o library

NetworkX_ 1.0 - 1.4
  Python package for working with complex networks.

NumPy_ 1.3 - 1.6

SciPy_ 0.7 - 0.10
  Numpy and Scipy are high-level, optimized scientific computing libraries.

Enthought_ Traits_ 4.0.0

.. note::

    Full distributions such as pythonxy_ or EPD_ provide the above packages,
    except Nibabel_.

Strong Recommendations
~~~~~~~~~~~~~~~~~~~~~~

IPython_ 0.10.2 - 0.12
  Interactive python environment. This is necessary for some parallel
  components of the pipeline engine.

Matplotlib_ 1.0
  Plotting library

Sphinx_
  Required for building the documentation

`Graphviz <http://www.graphviz.org/>`_
  Required for building the documentation

Interface Dependencies
~~~~~~~~~~~~~~~~~~~~~~

These are the software packages that mindboggle.interfaces wraps:

FreeSurfer_
  FreeSurfer version 4 and higher
  
nipy_ 
  0.1.2+20110404 or later

