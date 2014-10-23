#!/bin/bash

# This script provides directions for installing Mindboggle
# (http://mindboggle.info) and its dependencies.
# You need a good Internet connection.
# NOTE: We highly recommend installing Mindboggle as a virtual machine
#       (http://mindboggle.info/users/INSTALL.html).
#
# Copy this file into a setup directory and run::
#
#    mkdir /Users/arno/mindboggle_setup
#    cp build_mindboggle.sh /Users/arno/mindboggle_setup/
#    cd /Users/arno/mindboggle_setup/
#    bash build_mindboggle.sh

CWD=$(pwd)

#-----------------------------------------------------------------------------
# Install Anaconda Python distribution:
#
# Scientific Python environment
# -----------------------------
# Easily-installed Python distributions are available such as Anaconda,
# Canopy, and PythonXY that include components necessary to run Mindboggle:
# Python, NumPy, SciPy, Networkx.  These distributions also include
# packages that benefit the user, such as IPython (parallel processing) and
# Sphinx (generating documentation).
#
# Your system must have basic compilation software (cmake_ and a C++ compiler
# such as gcc), pip for installing Python libraries, and the Visualization
# Toolkit (after installing the Anaconda Python distribution).
#-----------------------------------------------------------------------------
wget http://repo.continuum.io/miniconda/Miniconda-2.2.2-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
export PATH=$CWD/anaconda/bin:$PATH

#-----------------------------------------------------------------------------
# Install Nipype and remaining dependencies:
#
# Nipype is a software pipeline framework.
# Nibabel is a medical image read/write library.
#-----------------------------------------------------------------------------
conda install --yes pip cmake
conda install --yes numpy scipy nose traits networkx
conda install --yes dateutil ipython-notebook matplotlib
pip install nibabel --use-mirrors
pip install https://github.com/RDFLib/rdflib/archive/master.zip
pip install https://github.com/satra/prov/archive/rdf.zip
pip install https://github.com/nipy/nipype/archive/master.zip

#-----------------------------------------------------------------------------
# Install VTK (for ANTs and Mindboggle):
#-----------------------------------------------------------------------------
conda install --yes vtk
VTK_DIR=$CWD/anaconda/lib/vtk-5.10  # Needed for ANTs
export PATH=$VTK_DIR:$PATH

#-----------------------------------------------------------------------------
# Install compiling utilities:
#-----------------------------------------------------------------------------
sudo apt-get -y update
sudo apt-get -y install g++
sudo apt-get -y install make
sudo apt-get -y install git
sudo apt-get -y install xorg openbox

#-----------------------------------------------------------------------------
# Install Mindboggle (pip for code, git for C++ code):
#
# To check out the latest development version using git:
# git clone git://github.com/binarybottle/mindboggle.git
#
# Step 1: Install Mindboggle Python software
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Change into the ``mindboggle`` directory and install mindboggle
# (depending on permissions you may need to use ``sudo``)::
#    python setup.py install
#
# Step 2: Compile Mindboggle C++ code
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In the mindboggle_tools/bin directory, run cmake (or ccmake) and make::
#    cd mindboggle/mindboggle_tools/bin/
#    cmake ..
#    make
#
# Step 3: Set environment variables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For your system to recognize the C++ code we just compiled, include
# the following lines in your system's equivalent of a .bash_profile file.
#-----------------------------------------------------------------------------
pip install https://github.com/binarybottle/mindboggle/archive/master.zip
git clone https://github.com/binarybottle/mindboggle.git
MINDBOGGLE_TOOLS=$CWD/mindboggle/mindboggle_tools/bin
export PATH=$MINDBOGGLE_TOOLS:$PATH
cd $MINDBOGGLE_TOOLS
cmake ../ -DVTK_DIR:STRING=$CWD/anaconda/lib/vtk-5.10
make
cd $CWD

#-----------------------------------------------------------------------------
# Non-Python software packages
# ----------------------------
# Use FreeSurfer_ and ANTS_ to process T1-weighted MRI data for use by
# Mindboggle.  See README (http://mindboggle.info/users/README.html)
# for instructions on running recon-all in FreeSurfer and
# antsCorticalThickness.sh in ANTs.
#-----------------------------------------------------------------------------
# Install FreeSurfer:
#
# FreeSurfer (http://surfer.nmr.mgh.harvard.edu)
# provides labeled cortical surfaces and non/cortical volumes.
# (Note: Mindboggle expects the FreeSurfer SUBJECTS_DIR environment
# variable to be set.)
#-----------------------------------------------------------------------------
# http://surfer.nmr.mgh.harvard.edu/fswiki/Download
wget -c ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/5.3.0/freesurfer-Linux-centos4_x86_64-stable-pub-v5.3.0.tar.gz
sudo tar xzvf freesurfer-Linux-centos4_x86_64-stable-pub-v5.3.0.tar.gz
rm freesurfer-Linux-centos4_x86_64-stable-pub-v5.3.0.tar.gz
FREESURFER_HOME=$CWD/freesurfer
SUBJECTS_DIR=$FREESURFER_HOME/subjects
export PATH=$FREESURFER_HOME:$SUBJECTS_DIR:$PATH
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Remove some of the larger directories:
cd freesurfer
sudo rm -r average fsfast matlab mni subjects tktools trctrain
cd lib
sudo rm -r cuda qt tcl tcltktixblt
cd $CWD

#-----------------------------------------------------------------------------
# Install ANTs:
#
# ANTs (http://brianavants.wordpress.com/2012/04/13/updated-ants-compile-instructions-april-12-2012/)
# optionally provides tissue segmentation, affine registration to standard
# space, and nonlinear registration for whole-brain labeling, to improve and
# extend Mindboggle results.
# (Note: After installing ANTS, be sure to copy the files in ANTs/Scripts/
# to the antsbin/bin/ directory, and set the ANTSPATH environment variable.)
#-----------------------------------------------------------------------------
# http://brianavants.wordpress.com/2012/04/13/
#        updated-ants-compile-instructions-april-12-2012/
git clone git://github.com/stnava/ANTs.git
mkdir antsbin
cd antsbin
cmake ../ANTs -DVTK_DIR:STRING=$CWD/anaconda/lib/vtk-5.10
make #-j 4
cp ../ANTs/Scripts/* bin/
cd $CWD
export ANTSPATH=$CWD/antsbin/bin/
