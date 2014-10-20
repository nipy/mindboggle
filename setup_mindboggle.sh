#!/bin/bash

# This script provides directions for installing Mindboggle
# (http://mindboggle.info) and its dependencies.
# You need a good Internet connection.
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
#-----------------------------------------------------------------------------
wget http://repo.continuum.io/miniconda/Miniconda-2.2.2-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
export PATH=$CWD/anaconda/bin:$PATH

#-----------------------------------------------------------------------------
# Install Nipype dependencies:
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
# Install FreeSurfer:
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
