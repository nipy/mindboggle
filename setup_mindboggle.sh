#!/bin/bash
#=============================================================================
# This script provides directions for installing Mindboggle and dependencies
# (http://mindboggle.info).  Running it requires a good Internet connection.
# We highly recommend installing Mindboggle as a virtual machine
# (http://mindboggle.info/users/INSTALL.html).
# Tested on an Ubuntu 14.04 machine.
#
# Usage:
#     bash setup_mindboggle.sh
#
#     Or with arguments:
#     bash setup_mindboggle.sh <download_dir> <install_dir> <env>
#                              <os> <ants> <sudo>
#
#     For example:
#     bash setup_mindboggle.sh /home/vagrant/downloads \
#                              /home/vagrant/install \
#                              /home/vagrant/.bash_profile \
#                              Linux  0  1
# Note:
#     All arguments are optional:
#     <download_dir>, <install_dir>, and <env> will be created if required.
#     <env> is a global environment sourcing script
#           to set environment variables, such as .bash_profile.
#     <os> is the operating system, either "Linux" or "MacOSX".
#     <ants> is set to 1 or 0, to install/enable ANTS or not.
#     <sudo> is set to 1 or 0, to install system-wide dependencies on Linux
#            using sudo or not.
#
# Authors:
#     - Daniel Clark, 2014
#     - Arno Klein, 2014-2016  (arno@mindboggle.info)  http://binarybottle.com
#
# Copyright 2016,  Mindboggle team, Apache v2.0 License
#=============================================================================

#-----------------------------------------------------------------------------
# Assign download and installation path arguments:
#-----------------------------------------------------------------------------
DOWNLOAD=$1
INSTALL=$2
ENV=$3
OS=$4
ANTS=$5
SUDO=$6

#-----------------------------------------------------------------------------
# Create folders and file if they don't exist:
#-----------------------------------------------------------------------------
if [ -z "$DOWNLOAD" ]; then
    DOWNLOAD="${HOME}/downloads"
fi
if [ ! -d $DOWNLOAD ]; then
    mkdir -p $DOWNLOAD;
fi

if [ -z "$INSTALL" ]; then
    INSTALL="${HOME}/install"
fi
if [ ! -d $INSTALL ]; then
    mkdir -p $INSTALL;
fi

if [ -z "$ENV" ]; then
    ENV="${HOME}/.bash_profile"
fi
if [ ! -e "$ENV" ] ; then
    touch "$ENV"
fi
if [ ! -w "$ENV" ] ; then
    echo cannot write to $ENV
    exit 1
fi
if [ -z "$OS" ]; then
    OS="Linux"
fi
if [ -z "$ANTS" ]; then
    ANTS=0
fi
if [ -z "$SUDO" ]; then
    SUDO=1
fi

#-----------------------------------------------------------------------------
# System-wide dependencies:
#-----------------------------------------------------------------------------
if [ $OS = "Linux" ]; then
    if [ $SUDO -eq 1 ]; then
        sudo apt-get update
        sudo apt-get install -y g++ git make xorg
    else
        apt-get update
        apt-get install -y g++ git make xorg
    fi
fi

#-----------------------------------------------------------------------------
# Anaconda's miniconda Python distribution for local installs:
#-----------------------------------------------------------------------------
CONDA_URL="http://repo.continuum.io/miniconda"
CONDA_FILE="Miniconda-latest-${OS}-x86_64.sh"
CONDA_DL="${DOWNLOAD}/${CONDA_FILE}"
CONDA_PATH="${INSTALL}/miniconda2"
if [ $OS = "Linux" ]; then
    wget -O $CONDA_DL ${CONDA_URL}/$CONDA_FILE
else
    curl -o $CONDA_DL ${CONDA_URL}/$CONDA_FILE
fi

bash $CONDA_DL -b -p $CONDA_PATH

#-----------------------------------------------------------------------------
# Set up PATH:
#-----------------------------------------------------------------------------
export PATH=${INSTALL}/bin:$PATH
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
    if [ $SUDO -eq 1 ]; then
        sudo mkdir /usr/lib64
        sudo ln -s /usr/lib/x86_64-linux-gnu/libGLU.so.1 /usr/lib64/libGLU.so
        sudo ln -s /usr/lib/x86_64-linux-gnu/libSM.so.6 /usr/lib64/libSM.so
        sudo ln -s /usr/lib/x86_64-linux-gnu/libICE.so.6 /usr/lib64/libICE.so
        sudo ln -s /usr/lib/x86_64-linux-gnu/libX11.so.6 /usr/lib64/libX11.so
        sudo ln -s /usr/lib/x86_64-linux-gnu/libXext.so.6 /usr/lib64/libXext.so
        sudo ln -s /usr/lib/x86_64-linux-gnu/libXt.so.6 /usr/lib64/libXt.so
        sudo ln -s /usr/lib/x86_64-linux-gnu/mesa/libGL.so.1 /usr/lib64/libGL.so
    else
        mkdir /usr/lib64
        ln -s /usr/lib/x86_64-linux-gnu/libGLU.so.1 /usr/lib64/libGLU.so
        ln -s /usr/lib/x86_64-linux-gnu/libSM.so.6 /usr/lib64/libSM.so
        ln -s /usr/lib/x86_64-linux-gnu/libICE.so.6 /usr/lib64/libICE.so
        ln -s /usr/lib/x86_64-linux-gnu/libX11.so.6 /usr/lib64/libX11.so
        ln -s /usr/lib/x86_64-linux-gnu/libXext.so.6 /usr/lib64/libXext.so
        ln -s /usr/lib/x86_64-linux-gnu/libXt.so.6 /usr/lib64/libXt.so
        ln -s /usr/lib/x86_64-linux-gnu/mesa/libGL.so.1 /usr/lib64/libGL.so
    fi
fi

#-----------------------------------------------------------------------------
# Python packages:
#-----------------------------------------------------------------------------
conda install --yes numpy scipy matplotlib pandas nose networkx traits vtk ipython
pip install nibabel nipype

#-----------------------------------------------------------------------------
# Mindboggle:
#-----------------------------------------------------------------------------
MB_DL=${DOWNLOAD}/mindboggle
git clone https://github.com/nipy/mindboggle.git $MB_DL
mv $MB_DL $INSTALL
cd ${INSTALL}/mindboggle
python setup.py install --prefix=$INSTALL
cd ${INSTALL}/mindboggle/surface_cpp_tools/bin
cmake ../  # -DVTK_DIR:STRING=$VTK_DIR
make
cd $INSTALL

#-----------------------------------------------------------------------------
# ANTs (http://brianavants.wordpress.com/2012/04/13/
#              updated-ants-compile-instructions-april-12-2012/)
# antsCorticalThickness.h pipeline optionally provides tissue segmentation,
# affine registration to standard space, and nonlinear registration for
# whole-brain labeling, to improve and extend Mindboggle results.
#-----------------------------------------------------------------------------
if [ $ANTS -eq 1 ]; then
    ANTS_DL=${DOWNLOAD}/ants
    git clone https://github.com/stnava/ANTs.git $ANTS_DL
    cd $ANTS_DL
    git checkout tags/v2.1.0rc2
    mkdir ${INSTALL}/ants
    cd ${INSTALL}/ants
    cmake $ANTS_DL  # -DVTK_DIR:STRING=$VTK_DIR
    make
    cp -r ${ANTS_DL}/Scripts/* ${INSTALL}/ants/bin
    # Remove non-essential directories:
    mv ${INSTALL}/ants/bin ${INSTALL}/ants_bin
    rm -rf ${INSTALL}/ants/*
    mv ${INSTALL}/ants_bin ${INSTALL}/ants/bin
fi

#-----------------------------------------------------------------------------
# Set environment variables:
#-----------------------------------------------------------------------------

# -- Local install --
echo "# Local install prefix" >> $ENV
echo "export PATH=${INSTALL}/bin:\$PATH" >> $ENV

# -- Mindboggle --
echo "# Mindboggle" >> $ENV
echo "export surface_cpp_tools=${INSTALL}/mindboggle/surface_cpp_tools/bin" >> $ENV
echo "export PATH=\$surface_cpp_tools:\$PATH" >> $ENV
echo "export PYTHONPATH=\$PYTHONPATH:\${INSTALL}/mindboggle" >> $ENV

# -- ANTs --
if [ $ANTS -eq 1 ]; then
    echo "# ANTs" >> $ENV
    echo "export ANTSPATH=${INSTALL}/ants/bin" >> $ENV
    echo "export PATH=\$ANTSPATH:\$PATH" >> $ENV
fi

source $ENV

#-----------------------------------------------------------------------------
# Finally, remove non-essential directories:
#-----------------------------------------------------------------------------
rm_extras=0
if [ $rm_extras -eq 1 ]; then
    rm -r ${DOWNLOAD}/*
fi

