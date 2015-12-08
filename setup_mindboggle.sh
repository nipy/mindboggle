#!/bin/bash
#=============================================================================
# This script provides directions for installing Mindboggle and dependencies
# (http://mindboggle.info).  Running it requires a good Internet connection.
# We highly recommend installing Mindboggle as a virtual machine
# (http://mindboggle.info/users/INSTALL.html).
#
# Usage:
#     ./install_mindboggle <download_dir> <install_dir> <env_file> <os> <vtk> <ants>
#
#     For example:
#     ./install_mindboggle /home/arno/downloads \
#                          /home/arno/software/install \
#                          /home/arno/.bash_profile \
#                          osx \
#                          vtk-6.3.0-py27_0
#                          0
#
# Note:
#     <download_dir>, <install_dir>, <env_file> must exist and be full paths.
#     <env_file> is a global environment sourcing script
#                to set environment variables, such as .bash_profile.
#     <os> can be any one of {linux, osx}.
#     <vtk> is the name of the vtk library, with VTK and Python versions.
#     <ants> can be set to 1 or 0, to install ANTS or not.
#
# Authors:
#     - Daniel Clark, 2014
#     - Arno Klein, 2014-2015  (arno@mindboggle.info)  http://binarybottle.com
#
# Copyright 2015,  Mindboggle team, Apache v2.0 License
#=============================================================================

#-----------------------------------------------------------------------------
# Assign download and installation path arguments:
#-----------------------------------------------------------------------------
#DL_PREFIX=/homedir/downloads
#INSTALL_PREFIX=/homedir/software/install
#MB_ENV=/homedir/.bash_profile
DL_PREFIX=$1
INSTALL_PREFIX=$2
MB_ENV=$3
OS=$4
VTK=$5
ANTS=$6
OS_XTRA=x86_64

#-----------------------------------------------------------------------------
# System-wide dependencies:
#-----------------------------------------------------------------------------
if [ $OS = "linux" ]; then
    apt-get update
    apt-get install -y g++ git make xorg
    OS_STR=Linux
else
    OS_STR=MacOSX
fi

#-----------------------------------------------------------------------------
# Anaconda's miniconda Python distribution for local installs:
#-----------------------------------------------------------------------------
CONDA_FILE=Miniconda-latest-${OS_STR}-${OS_XTRA}.sh
CONDA_DL=${DL_PREFIX}/${CONDA_FILE}
if [ $OS = "linux" ]; then
    wget -O $CONDA_DL http://repo.continuum.io/miniconda/${CONDA_FILE}
else
    curl -o $CONDA_DL http://repo.continuum.io/miniconda/${CONDA_FILE}

#chmod +x $CONDA_DL
bash $CONDA_DL -b -p $INSTALL_PREFIX/miniconda

# Setup PATH
export PATH=${INSTALL_PREFIX}/bin:$PATH

#-----------------------------------------------------------------------------
# Additional resources for installing packages:
#-----------------------------------------------------------------------------
conda install --yes cmake pip

#-----------------------------------------------------------------------------
# Python packages:
#-----------------------------------------------------------------------------
conda install --yes numpy scipy matplotlib pandas nose networkx traits vtk ipython
pip install nibabel nipype
VTK_DIR=${INSTALL_PREFIX}/miniconda/pkgs/${VTK}

#-----------------------------------------------------------------------------
# Mindboggle:
#-----------------------------------------------------------------------------
MB_DL=${DL_PREFIX}/mindboggle
git clone https://github.com/binarybottle/mindboggle.git $MB_DL
mv ${MB_DL} ${INSTALL_PREFIX}
cd ${INSTALL_PREFIX}/mindboggle
python setup.py install --prefix=${INSTALL_PREFIX}
cd ${INSTALL_PREFIX}/mindboggle/surface_cpp_tools/bin
cmake ../ -DVTK_DIR:STRING=${VTK_DIR}
make
cd ${INSTALL_PREFIX}

#-----------------------------------------------------------------------------
# ANTs (http://brianavants.wordpress.com/2012/04/13/
#              updated-ants-compile-instructions-april-12-2012/)
# antsCorticalThickness.h pipeline optionally provides tissue segmentation,
# affine registration to standard space, and nonlinear registration for
# whole-brain labeling, to improve and extend Mindboggle results.
#-----------------------------------------------------------------------------
if [ $ANTS = 1 ]; then
    ANTS_DL=${DL_PREFIX}/ants
    git clone https://github.com/stnava/ANTs.git $ANTS_DL
    cd $ANTS_DL
    git checkout tags/v2.1.0rc2
    mkdir ${INSTALL_PREFIX}/ants
    cd ${INSTALL_PREFIX}/ants
    cmake $ANTS_DL #-DVTK_DIR:STRING=${VTK_DIR}
    make
    cp -r ${ANTS_DL}/Scripts/* ${INSTALL_PREFIX}/ants/bin
    # Remove non-essential directories:
    mv ${INSTALL_PREFIX}/ants/bin ${INSTALL_PREFIX}/ants_bin
    rm -rf ${INSTALL_PREFIX}/ants/*
    mv ${INSTALL_PREFIX}/ants_bin ${INSTALL_PREFIX}/ants/bin
fi

#-----------------------------------------------------------------------------
# Set environment variables:
#-----------------------------------------------------------------------------
# -- Local install --
echo "# Local install prefix" >> $MB_ENV
echo "export PATH=${INSTALL_PREFIX}/bin:\$PATH" >> $MB_ENV

# -- Mindboggle --
echo "# Mindboggle" >> $MB_ENV
echo "export surface_cpp_tools=${INSTALL_PREFIX}/mindboggle/surface_cpp_tools/bin" >> $MB_ENV
echo "export PATH=\$surface_cpp_tools:\$PATH" >> $MB_ENV

# -- ANTs --
if [ $ANTS = 1 ]; then
    echo "# ANTs" >> $MB_ENV
    echo "export ANTSPATH=${INSTALL_PREFIX}/ants/bin" >> $MB_ENV
    echo "export PATH=\$ANTSPATH:\$PATH" >> $MB_ENV
fi

#-----------------------------------------------------------------------------
# Finally, remove non-essential directories:
#-----------------------------------------------------------------------------
rm_extras=0
if [ $rm_extras == 1 ]
    rm ${DL_PREFIX}/* -rf
fi
