.. _INSTALL:

============
Installation  **(UNDER CONSTRUCTION)**
============
Mindboggle comes as a single installation program, `setup_mindboggle <https://github.com/binarybottle/mindboggle/blob/master/setup_mindboggle>`_,
that creates a virtual machine where everything is installed and
configured to run Mindboggle.  See below for directions on its use.
(For those who just *have* to install from source, see guidelines in
`setup_mindboggle.sh <https://github.com/binarybottle/mindboggle/blob/master/setup_mindboggle.sh>`_.)

Dependencies
------------
`Vagrant <http://www.vagrantup.com/downloads.html>`_ manages virtual machines.
    Vagrant provides reproducible and portable work environments
    that isolate dependencies and their configuration within a single
    disposable, consistent environment that can run on
    Linux, Mac OS X, or Windows.

`Virtualbox <https://www.virtualbox.org/wiki/Downloads>`_ provides virtual machines used by Vagrant.
    Alternative backend providers for Vagrant include VMware and Amazon Web Services.

Software to generate input data for Mindboggle
----------------------------------------------
Mindboggle takes in outputs from the recon-all command in `FreeSurfer <http://surfer.nmr.mgh.harvard.edu>`_
(`download <http://surfer.nmr.mgh.harvard.edu>`_)
and the ``antsCorticalThickness.sh`` command in `ANTs <http://stnava.github.io/ANTs/>`_
(`download <http://brianavants.wordpress.com/2012/04/13/updated-ants-compile-instructions-april-12-2012/>`_).
See `README <http://mindboggle.info/users/README.html>`_ for
instructions on running these commands.

Configure virtual machine with mounted directories
-----------------------------------------------------------------------------
Mindboggle's single installation program,
`setup_mindboggle <https://github.com/binarybottle/mindboggle/blob/master/setup_mindboggle>`_,
generates a configuration file called ``Vagrantfile`` for the mindboggle virtual machine.
Copy ``setup_mindboggle`` to wherever you want the startup directory to be,
and make sure it is executable by typing in a terminal::

    chmod +x setup_mindboggle

Configure the virtual machine to access your local FreeSurfer subjects
directory, your local ANTs subjects directory,
and to use a given number of processors (backslashes denote line returns)::

    python setup_mindboggle --install \
                            --home /home/vagrant \
                            --freesurfer /Applications/freesurfer/subjects \
                            --ants /data/antsCorticalThickness/subjects \
                            --proc 6

For help with other setup options::

    python setup_mindboggle --help

See `README <http://mindboggle.info/users/README.html#set-up-mindboggle>`_
for instructions on how to run Mindboggle with this setup!
