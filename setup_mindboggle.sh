#!/bin/bash
#=============================================================================
# This script provides directions for installing Mindboggle and dependencies
# (http://mindboggle.info) on a Linux machine.
# Running it requires a good Internet connection.
# We highly recommend installing Mindboggle as a virtual machine
# (http://mindboggle.info/users/INSTALL.html).
# This script has been tested on an Ubuntu 14.04 machine.
#
# This script assumes that if it is not installing ANTs, then the global
# environment file already includes the path to the ANTs bin directory.
# ANTs is used to perform affine registration to standard (MNI152) space,
# refine gray/white matter segmentation, and perform nonlinear volume
# registration for whole-brain labeling.
#
# Usage:
#     source ./setup_mindboggle.sh
#
#     Or with arguments:
#     source ./setup_mindboggle.sh
#               <absolute path to download directory (create if empty)>
#               <absolute path to install directory (create if empty)>
#               <absolute path to environment file (or create .bash_profile)>
#               <install ants? set to 1 (default is 0)>
#
#     Example:
#     source ./setup_mindboggle.sh /home/vagrant/downloads \
#            /home/vagrant/install /home/vagrant/.bash_profile 0
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
ANTS=$4
SUDO=1
OS=Linux

#-----------------------------------------------------------------------------
# Create folders and file if they don't exist:
#-----------------------------------------------------------------------------
if [ -z "$DOWNLOAD" ]; then
    DOWNLOAD="$HOME/downloads"
fi
if [ ! -d $DOWNLOAD ]; then
    mkdir -p $DOWNLOAD;
fi

if [ -z "$INSTALL" ]; then
    INSTALL="$HOME/install"
fi
if [ ! -d $INSTALL ]; then
    mkdir -p $INSTALL;
fi

if [ -z "$ENV" ]; then
    ENV="$HOME/.bash_profile"
fi
if [ ! -e "$ENV" ] ; then
    touch "$ENV"
fi
if [ ! -w "$ENV" ] ; then
    echo cannot write to $ENV
    exit 1
fi
if [ -z "$ANTS" ]; then
    ANTS=0
fi
#if [ -z "$SUDO" ]; then
#    SUDO=1
#fi
#if [ -z "$OS" ]; then
#    OS="Linux"
#fi

#-----------------------------------------------------------------------------
# Install system-wide dependencies in linux:
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
# Install Anaconda's latest miniconda Python distribution:
#-----------------------------------------------------------------------------
CONDA_URL="http://repo.continuum.io/miniconda"
CONDA_FILE="Miniconda-latest-$OS-x86_64.sh"
CONDA_DL="$DOWNLOAD/$CONDA_FILE"
CONDA_PATH="$INSTALL/miniconda2"
if [ $OS = "Linux" ]; then
    wget -O $CONDA_DL $CONDA_URL/$CONDA_FILE
else
    curl -o $CONDA_DL $CONDA_URL/$CONDA_FILE
fi
bash $CONDA_DL -b -p $CONDA_PATH

# Set environment variables:
echo "# Conda" >> $ENV
#echo "export PATH=\$INSTALL/bin:\$PATH" >> $ENV
echo "export PATH=$CONDA_PATH/bin:\$PATH" >> $ENV
source $ENV

#-----------------------------------------------------------------------------
# Fix paths to linux libraries using symbolic links:
#-----------------------------------------------------------------------------
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
# Use conda to install additional resources for installing packages:
#-----------------------------------------------------------------------------
conda install --yes cmake pip

#-----------------------------------------------------------------------------
# Use conda and pip to install the latest Python packages:
#-----------------------------------------------------------------------------
conda install --yes numpy scipy matplotlib pandas networkx vtk ipython
pip install nibabel nipype

#-----------------------------------------------------------------------------
# Install the latest Mindboggle:
#-----------------------------------------------------------------------------
surface_cpp_tools=$INSTALL/mindboggle/surface_cpp_tools/bin
git clone https://github.com/nipy/mindboggle.git $INSTALL/mindboggle
cd $INSTALL/mindboggle
python setup.py install  #--prefix=$INSTALL
cd $surface_cpp_tools
cmake ../  # -DVTK_DIR:STRING=$VTK_DIR
make

# Set environment variables:
echo "# Mindboggle" >> $ENV
echo "export PATH=$surface_cpp_tools:\$PATH" >> $ENV
source $ENV

#-----------------------------------------------------------------------------
# Install ANTs v2.1.0rc3
# The antsCorticalThickness.h pipeline optionally provides gray/white matter
# segmentation, affine registration to standard space, and nonlinear volume
# registration for whole-brain labeling, to improve Mindboggle results.
#-----------------------------------------------------------------------------
if [ $ANTS -eq 1 ]; then
    ANTS_DL=$DOWNLOAD/ants
    git clone https://github.com/stnava/ANTs.git $ANTS_DL
    cd $ANTS_DL
    git checkout tags/v2.1.0rc3
    mkdir $INSTALL/ants
    cd $INSTALL/ants
    cmake $ANTS_DL  # -DVTK_DIR:STRING=$VTK_DIR
    make
    cp -r $ANTS_DL/Scripts/* $INSTALL/ants/bin
    # Remove non-essential directories:
    mv $INSTALL/ants/bin $INSTALL/ants_bin
    rm -rf $INSTALL/ants/*
    mv $INSTALL/ants_bin $INSTALL/ants/bin

    # Set environment variables:
    echo "# ANTs" >> $ENV
    echo "export ANTSPATH=$INSTALL/ants/bin" >> $ENV
    echo "export PATH=$ANTSPATH:\$PATH" >> $ENV
    source $ENV
fi

#-----------------------------------------------------------------------------
# Remove non-essential directories:
#-----------------------------------------------------------------------------
rm_extras=0
if [ $rm_extras -eq 1 ]; then
    rm -r $DOWNLOAD/*
fi

