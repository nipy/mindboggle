#!/bin/bash
#=============================================================================
# This script provides directions for installing Mindboggle and dependencies
# (http://mindboggle.info) on a Linux machine (tested on Ubuntu 12.04, 14.04)
# with Python 3. Running it requires a good Internet connection.
#
# This script assumes that if it isn't installing ANTs, the global environment
# file .bash_profile already includes the path to the ANTs bin directory.
# ANTs is used to perform affine registration to standard (MNI152) space,
# refine gray/white matter segmentation, and perform nonlinear volume
# registration for whole-brain labeling. Installing and running ANTs
# requires considerable (> 1 GB) RAM.
#
# Usage:
#     source ./install_mindboggle.sh
#
#     Or with arguments:
#     source ./install_mindboggle.sh
#               <absolute path to download directory (create if empty)>
#               <absolute path to install directory (create if empty)>
#               <absolute path to environment file (or create .bash_profile)>
#               <absolute path to home directory (default is /home/vagrant)>
#               <install ants: yes or no? (default is yes)>
#
#     Example:
#     source ./install_mindboggle.sh /home/vagrant/downloads \
#            /home/vagrant/install /home/vagrant/.bash_profile yes
#
# Authors:
#     - Daniel Clark, 2014
#     - Arno Klein, 2014-2016  (arno@mindboggle.info)  http://binarybottle.com
#
# Copyright 2016,  Mindboggle team, Apache v2.0 License
#=============================================================================

#-----------------------------------------------------------------------------
# Path arguments:
#-----------------------------------------------------------------------------
DOWNLOAD=$1
INSTALL=$2
ENV=$3
HOME=$4
ANTS=$5
VTK_DIR=$6

#-----------------------------------------------------------------------------
# OS and sudo:
#-----------------------------------------------------------------------------
OS=Linux
SUDO=1

#-----------------------------------------------------------------------------
# Create folders and file if they don't exist:
#-----------------------------------------------------------------------------
if [ -z "$HOME" ]; then
    HOME="/home/vagrant"
fi
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
    ANTS="no"
fi
if [ -z "$OS" ]; then
    OS="Linux"
fi
if [ -z "$SUDO" ]; then
    SUDO=1
fi

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
CONDA_URL="https://repo.continuum.io/miniconda"
CONDA_FILE="Miniconda3-latest-$OS-x86_64.sh"
CONDA_DL="$DOWNLOAD/$CONDA_FILE"
CONDA_PATH="$INSTALL/miniconda3"
if [ $OS = "Linux" ]; then
    wget -O $CONDA_DL $CONDA_URL/$CONDA_FILE
else
    curl -o $CONDA_DL $CONDA_URL/$CONDA_FILE
fi
bash $CONDA_DL -b -p $CONDA_PATH

# Set environment variables:
echo "# Conda" >> $ENV
echo "export PATH=$CONDA_PATH/bin:\$PATH" >> $ENV
source $ENV

#-----------------------------------------------------------------------------
# Fix paths to Linux libraries using symbolic links:
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
pip install --upgrade pip

#-----------------------------------------------------------------------------
# Install VTK:
#-----------------------------------------------------------------------------
conda install -c https://conda.anaconda.org/clinicalgraphics vtk

#-----------------------------------------------------------------------------
# Use conda and pip to install the latest Python packages:
#-----------------------------------------------------------------------------
conda install --yes numpy scipy matplotlib pandas networkx ipython
pip install nibabel

#-----------------------------------------------------------------------------
# Install nipype:
# prov requires lxml which requires libxml
#-----------------------------------------------------------------------------
conda install lxml
pip install prov traits nose future simplejson
git clone https://github.com/nipy/nipype.git $INSTALL/nipype
cd $INSTALL/nipype
python setup.py install

# pip won't install all nipype's dependencies, so use conda
# or you will have to do it manually. prov requires lxml which requires libxml
#pip install nipype
#pip install --upgrade https://github.com/nipy/nipype/archive/master.zip

#-----------------------------------------------------------------------------
# Install the latest Mindboggle:
#-----------------------------------------------------------------------------
vtk_cpp_tools=$INSTALL/mindboggle/vtk_cpp_tools/bin
git clone https://github.com/nipy/mindboggle.git $INSTALL/mindboggle
cd $INSTALL/mindboggle
python setup.py install
mkdir $vtk_cpp_tools
cd $vtk_cpp_tools
if [ -z "$VTK_DIR" ]; then
    #VTK_DIR="$INSTALL/VTK"
    VTK_DIR="$INSTALL/lib/cmake/vtk-7.1/Modules"
    cmake ../ -DVTK_DIR:STRING=$VTK_DIR
else
    cmake ../
fi
make

# Set environment variables:
echo "# Mindboggle" >> $ENV
echo "export vtk_cpp_tools=$vtk_cpp_tools" >> $ENV
echo "export PATH=$vtk_cpp_tools:\$PATH" >> $ENV
source $ENV

#-----------------------------------------------------------------------------
# Install ANTs v2.1.0rc3
# The antsCorticalThickness.h pipeline optionally provides gray/white matter
# segmentation, affine registration to standard space, and nonlinear volume
# registration for whole-brain labeling, to improve Mindboggle results.
#-----------------------------------------------------------------------------
if [ $ANTS = "yes" ]; then
    ANTS_DL=$DOWNLOAD/ants
    ANTSPATH=$INSTALL/ants/bin
    git clone https://github.com/stnava/ANTs.git $ANTS_DL
    cd $ANTS_DL
    git checkout tags/v2.1.0rc3
    mkdir $INSTALL/ants
    cd $INSTALL/ants
    cmake $ANTS_DL  # -DVTK_DIR:STRING=$VTK_DIR
    make
    cp -r $ANTS_DL/Scripts/* $ANTSPATH

    # Set environment variables:
    echo "# ANTs" >> $ENV
    echo "export ANTSPATH=$ANTSPATH" >> $ENV
    echo "export PATH=$ANTSPATH:\$PATH" >> $ENV
    source $ENV
fi

#-----------------------------------------------------------------------------
# Remove non-essential directories
# (set to 0 to keep a complete box for easy git updates of ANTs):
#-----------------------------------------------------------------------------
rm_extras=1
if [ $rm_extras -eq 1 ]; then
    if [ $ANTS = "yes" ]; then
        if [ $SUDO -eq 1 ]; then
            sudo mv $ANTSPATH $INSTALL/ants_bin
            sudo rm -rf $INSTALL/ants/*
            sudo mv $INSTALL/ants_bin $ANTSPATH
        else
            mv $ANTSPATH $INSTALL/ants_bin
            rm -rf $INSTALL/ants/*
            mv $INSTALL/ants_bin $ANTSPATH
        fi
    fi
    #rm -r $DOWNLOAD/*
fi

