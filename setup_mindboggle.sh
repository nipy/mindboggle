#!/bin/bash
#=============================================================================
# This script provides directions for installing Mindboggle and dependencies
# (http://mindboggle.info).  Running it requires a good Internet connection.
# We highly recommend installing Mindboggle as a virtual machine
# (http://mindboggle.info/users/INSTALL.html).
#
# Usage:
#     ./install_mindboggle <download_dir> <install_dir>
#
# Note:
#     <download_dir> must already exist and <install_dir> must not exist.
#
# Authors:
#     - Daniel Clark, 2014
#     - Arno Klein, 2014  (arno@mindboggle.info)  http://binarybottle.com
#
# Copyright 2014,  Mindboggle team, Apache v2.0 License
#=============================================================================

#-----------------------------------------------------------------------------
# Assign download and installation path arguments:
#-----------------------------------------------------------------------------
DL_PREFIX=$1
INSTALL_PREFIX=$2

#-----------------------------------------------------------------------------
# System-wide dependencies:
#-----------------------------------------------------------------------------
apt-get update
apt-get install -y g++ git make xorg

#-----------------------------------------------------------------------------
# Anaconda's miniconda Python distribution for local installs:
#-----------------------------------------------------------------------------
CONDA_DL=${DL_PREFIX}/Miniconda-3.7.0-Linux-x86_64.sh
wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh -P $DL_PREFIX
chmod +x $CONDA_DL
$CONDA_DL -b -p $INSTALL_PREFIX
# Setup PATH
export PATH=${INSTALL_PREFIX}/bin:$PATH

#-----------------------------------------------------------------------------
# Additional resources for installing packages:
#-----------------------------------------------------------------------------
conda install --yes cmake pip

#-----------------------------------------------------------------------------
# Python packages:
#-----------------------------------------------------------------------------
conda install --yes numpy scipy matplotlib nose networkx traits vtk ipython
# Nipype: software pipeline framework; Nibabel: medical image read/write lib
pip install nibabel nipype
VTK_DIR=${INSTALL_PREFIX}/lib/vtk-5.10
#pip install https://github.com/RDFLib/rdflib/archive/master.zip
#pip install https://github.com/satra/prov/archive/rdf.zip

# nilearn.image.resample_img (which calls scipy.ndimage.affine_transform):
NL_DL=${DL_PREFIX}/nilearn
git clone https://github.com/nilearn/nilearn $NL_DL
cd $NL_DL
python setup.py install --prefix=${INSTALL_PREFIX}
cd $DL_PREFIX

#-----------------------------------------------------------------------------
# Mindboggle:
#-----------------------------------------------------------------------------
MB_DL=${DL_PREFIX}/mindboggle
git clone https://github.com/binarybottle/mindboggle.git $MB_DL
cd $MB_DL
python setup.py install --prefix=${INSTALL_PREFIX}
cd ${MB_DL}/mindboggle_tools/bin
cmake ../ -DVTK_DIR:STRING=${VTK_DIR}
make
cd $DL_PREFIX
cp -r ${MB_DL}/mindboggle_tools ${INSTALL_PREFIX}

#-----------------------------------------------------------------------------
# ANTs (http://brianavants.wordpress.com/2012/04/13/
#              updated-ants-compile-instructions-april-12-2012/)
# antsCorticalThickness.h pipeline optionally provides tissue segmentation,
# affine registration to standard space, and nonlinear registration for
# whole-brain labeling, to improve and extend Mindboggle results.
#-----------------------------------------------------------------------------
ANTS_DL=${DL_PREFIX}/ants
git clone https://github.com/stnava/ANTs.git $ANTS_DL
cd $ANTS_DL
git checkout tags/v2.1.0rc2
mkdir ${INSTALL_PREFIX}/ants
cd ${INSTALL_PREFIX}/ants
cmake $ANTS_DL -DVTK_DIR:STRING=${VTK_DIR}
make
cp -r ${ANTS_DL}/Scripts/* ${INSTALL_PREFIX}/ants/bin

#-----------------------------------------------------------------------------
# Create a global environment sourcing script and set environment variables:
#-----------------------------------------------------------------------------
MB_ENV=/etc/profile.d/mb_env.sh
touch $MB_ENV

# -- Local install --
echo "# Local install prefix" >> $MB_ENV
echo "export PATH=${INSTALL_PREFIX}/bin:\$PATH" >> $MB_ENV

# -- ANTs --
echo "# ANTs" >> $MB_ENV
echo "export ANTSPATH=${INSTALL_PREFIX}/ants/bin" >> $MB_ENV
echo "export PATH=\$ANTSPATH:\$PATH" >> $MB_ENV

# -- Mindboggle --
echo "# Mindboggle" >> $MB_ENV
echo "export MINDBOGGLE_TOOLS=${INSTALL_PREFIX}/mindboggle_tools/bin" >> $MB_ENV
echo "export PATH=\$MINDBOGGLE_TOOLS:\$PATH" >> $MB_ENV

#-----------------------------------------------------------------------------
# Finally, remove downloads directory and other non-essential directories:
#-----------------------------------------------------------------------------
rm ${DL_PREFIX}/* -rf

rm_extras=1
if [ $rm_extras == 1 ]

    mv ${INSTALL_PREFIX}/ants/bin ${INSTALL_PREFIX}/freesurfer/ants_bin
    rm -r ${INSTALL_PREFIX}/ants/*
    mv ${INSTALL_PREFIX}/freesurfer/ants_bin ${INSTALL_PREFIX}/ants/bin
fi