#!/bin/bash
#=============================================================================
# This script provides directions for installing Mindboggle and dependencies
# (http://mindboggle.info).  Running it requires a good Internet connection.
# We highly recommend installing Mindboggle as a virtual machine
# (http://mindboggle.info/users/INSTALL.html).
# Tested on Ubuntu 14.04 and MacOSX 10.11 x86_64 machines.
#
# Usage:
#     bash setup_mindboggle.sh <download_dir> <install_dir> <env> <os> <ants>
#
#     For example:
#     bash setup_mindboggle.sh /home/vagrant/downloads \
#                              /home/vagrant/install \
#                              /home/vagrant/.bash_profile \
#                              linux  0
# Note:
#     <download_dir>, <install_dir>, and <env> will be created locally
#                                              if they don't exist.
#     <env> is a global environment sourcing script
#           to set environment variables, such as .bash_profile.
#     <os> is the operating system, either "Linux" or "MacOSX".
#     <ants> is set to 1 or 0, to run ANTS or not.
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
ANTS=$5

#-----------------------------------------------------------------------------
# Create folders and file if they don't exist:
#-----------------------------------------------------------------------------
if [ -z "$DL_PREFIX" ]; then
    DL_PREFIX="${HOME}/downloads"
fi
if [ ! -d $DL_PREFIX ]; then
    mkdir -p $DL_PREFIX;
fi

if [ -z "$INSTALL_PREFIX" ]; then
    INSTALL_PREFIX="${HOME}/install"
fi
if [ ! -d $INSTALL_PREFIX ]; then
    mkdir -p $INSTALL_PREFIX;
fi

if [ -z "$MB_ENV" ]; then
    MB_ENV="${HOME}/.bash_profile"
fi
if [ ! -e "$MB_ENV" ] ; then
    touch "$MB_ENV"
fi
if [ ! -w "$MB_ENV" ] ; then
    echo cannot write to $MB_ENV
    exit 1
fi
if [ -z "$OS" ]; then
    OS="Linux"
fi
if [ -z "$ANTS" ]; then
    ANTS=0
fi

#-----------------------------------------------------------------------------
# System-wide dependencies:
#-----------------------------------------------------------------------------
if [ $OS = "Linux" ]; then
    apt-get update
    apt-get install -y g++ git make xorg
fi

#-----------------------------------------------------------------------------
# Anaconda's miniconda Python distribution for local installs:
#-----------------------------------------------------------------------------
CONDA_URL="http://repo.continuum.io/miniconda"
CONDA_FILE="Miniconda-latest-${OS}-x86_64.sh"
CONDA_DL="${DL_PREFIX}/${CONDA_FILE}"
CONDA_PATH="${INSTALL_PREFIX}/miniconda"
if [ $OS = "Linux" ]; then
    wget -O $CONDA_DL ${CONDA_URL}/$CONDA_FILE
else
    curl -o $CONDA_DL ${CONDA_URL}/$CONDA_FILE
fi

bash $CONDA_DL -b -p $CONDA_PATH

#-----------------------------------------------------------------------------
# Set up PATH:
#-----------------------------------------------------------------------------
export PATH=${INSTALL_PREFIX}/bin:$PATH
export PATH=${CONDA_PATH}/bin:$PATH

#-----------------------------------------------------------------------------
# Additional resources for installing packages:
#-----------------------------------------------------------------------------
conda install --yes cmake pip

# To avoid the following errors:
# "No rule to make target `/usr/lib/x86_64-linux-gnu/libGLU.so'"
# ...
# http://techtidings.blogspot.com/2012/01/problem-with-libglso-on-64-bit-ubuntu.html
if [ $OS = "Linux" ]; then
    mkdir /usr/lib64
    ln -s /usr/lib/x86_64-linux-gnu/libGLU.so.1 /usr/lib64/libGLU.so
    ln -s /usr/lib/x86_64-linux-gnu/libSM.so.6 /usr/lib64/libSM.so
    ln -s /usr/lib/x86_64-linux-gnu/libICE.so.6 /usr/lib64/libICE.so
    ln -s /usr/lib/x86_64-linux-gnu/libX11.so.6 /usr/lib64/libX11.so
    ln -s /usr/lib/x86_64-linux-gnu/libXext.so.6 /usr/lib64/libXext.so
    ln -s /usr/lib/x86_64-linux-gnu/libXt.so.6 /usr/lib64/libXt.so
    ln -s /usr/lib/x86_64-linux-gnu/mesa/libGL.so.1 /usr/lib64/libGL.so
fi

#-----------------------------------------------------------------------------
# Python packages:
#-----------------------------------------------------------------------------
conda install --yes numpy scipy matplotlib pandas nose networkx traits vtk ipython
pip install nibabel nipype
VTK_DIR=$CONDA_PATH

#-----------------------------------------------------------------------------
# Mindboggle:
#-----------------------------------------------------------------------------
MB_DL=${DL_PREFIX}/mindboggle
git clone https://github.com/nipy/mindboggle.git $MB_DL
mv $MB_DL $INSTALL_PREFIX
cd ${INSTALL_PREFIX}/mindboggle
python setup.py install --prefix=$INSTALL_PREFIX
cd ${INSTALL_PREFIX}/mindboggle/surface_cpp_tools/bin
cmake ../ -DVTK_DIR:STRING=$VTK_DIR
make
cd $INSTALL_PREFIX

#-----------------------------------------------------------------------------
# ANTs (http://brianavants.wordpress.com/2012/04/13/
#              updated-ants-compile-instructions-april-12-2012/)
# antsCorticalThickness.h pipeline optionally provides tissue segmentation,
# affine registration to standard space, and nonlinear registration for
# whole-brain labeling, to improve and extend Mindboggle results.
#-----------------------------------------------------------------------------
if [ $ANTS -eq 1 ]; then
    ANTS_DL=${DL_PREFIX}/ants
    git clone https://github.com/stnava/ANTs.git $ANTS_DL
    cd $ANTS_DL
    git checkout tags/v2.1.0rc2
    mkdir ${INSTALL_PREFIX}/ants
    cd ${INSTALL_PREFIX}/ants
    cmake $ANTS_DL #-DVTK_DIR:STRING=$VTK_DIR
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
if [ $ANTS -eq 1 ]; then
    echo "# ANTs" >> $MB_ENV
    echo "export ANTSPATH=${INSTALL_PREFIX}/ants/bin" >> $MB_ENV
    echo "export PATH=\$ANTSPATH:\$PATH" >> $MB_ENV
fi

source $MB_ENV

#-----------------------------------------------------------------------------
# Finally, remove non-essential directories:
#-----------------------------------------------------------------------------
rm_extras=1
if [ $rm_extras -eq 1 ]; then
    rm -r ${DL_PREFIX}/*
fi

