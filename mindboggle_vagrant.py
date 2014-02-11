#!/usr/bin/python
"""
Mindboggle setup script.

This program installs Mindboggle and its dependencies
in a Vagrant box (virtual machine).
You must have Vagrant (http://www.vagrantup.com) and
VirtualBox (http://www.virtualbox.org) installed.
The Mindboggle Python package automates shape analysis of
anatomical labels and features extracted from human brain
MR image data. See http://mindboggle.info
for up-to-date documentation.

Further setup is necessary within the Vagrant box:
FS, A...

Authors:
Arno Klein, 2014  .  arno@mindboggle.info  .  www.binarybottle.com

Copyright 2014,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
import os
import argparse
from string import Template

#=============================================================================
# Command-line arguments
#=============================================================================
parser = argparse.ArgumentParser(description="""
                    The Mindboggle Python package automates shape analysis of
                    anatomical labels and features extracted from human brain
                    MR image data. See http://mindboggle.info
                    for up-to-date documentation.
                    This program installs Mindboggle and its dependencies
                    in a Vagrant box (virtual machine).
                    You must have Vagrant (http://www.vagrantup.com) and
                    VirtualBox (http://www.virtualbox.org) installed.
                    Ex: python mindboggle_vagrant.py --num 1 --mem 4096
                        --cpu 50 --ants /data/antsCorticalThickness""",
                     formatter_class = lambda prog:
                     argparse.HelpFormatter(prog, max_help_position=40))

parser.add_argument("-o", "--out",
                    help='Output directory (default: $HOME/mindboggled)',
                    default=os.path.join(os.environ['HOME'],
                                         'mindboggled'), metavar='STR')
parser.add_argument("--freesurfer",
                    help=("FreeSurfer subjects directory (default if not "
                          "set: $SUBJECTS_DIR environment variable)"),
                    metavar='STR')
parser.add_argument("--ants",
                    help=("Optional directory containing "
                          "antsCorticalThickness.sh "
                          "output in individual subject directories"),
                    metavar='STR')
parser.add_argument("--atlases",
                    help=("Optional directory containing additional atlases"),
                    metavar='STR')
parser.add_argument("--cpu",
                    help=('Maximum CPU (percent): "--cpu 50" (default)'),
                    type=int, default=50, metavar='INT')
parser.add_argument("--mem",
                    help=('Maximum memory (MB): "--mem 4096" (default)'),
                    type=int, default=4096, metavar='INT')
parser.add_argument("--num",
                    help=('Number of processors: "--num 1" (default)'),
                    type=int, default=1, metavar='INT')

args = parser.parse_args()

print(args)

if not os.path.exists(args.out):
    os.mkdir(args.out)
out = 'config.vm.synced_folder "{0}", "mindboggled"'.format(args.out)
if args.freesurfer:
    freesurfer = 'config.vm.synced_folder "{0}", ' \
                 '"freesurfer_subjects"'.format(args.freesurfer)
else:
    freesurfer = 'config.vm.synced_folder "{0}", ' \
                 '"freesurfer_subjects"'.format(os.environ['SUBJECTS_DIR'])
if args.ants:
    ants = 'config.vm.synced_folder "{0}", ' \
           '"ants_subjects"'.format(args.ants)
else:
    ants = ''
if args.atlases:
    atlases = 'config.vm.synced_folder "{0}", ' \
              '"atlases"'.format(args.atlases)
else:
    atlases = ''
if args.cpu:
    cpu = args.cpu
else:
    cpu = '50'
if args.mem:
    mem = args.mem
else:
    mem = '4096'
if args.num:
    num = args.num
else:
    num = '1'

#=============================================================================
# Configurations to include in Vagrant file
#=============================================================================
vagrant_file = """# Vagrant file (http://www.vagrantup.com)
# (after https://github.com/nipy/nipype/blob/master/Vagrantfile)
# Vagrantfile API/syntax version:
VAGRANTFILE_API_VERSION = "2"

$script

# Configure Vagrant:
Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|

    # Every Vagrant virtual environment requires a box to build off of:
    config.vm.box = "mindboggle"
    config.vm.hostname = "boggler"

    # URL from where the 'config.vm.box' box will be fetched
    # if it doesn't already exist on the user's system:
    #config.vm.box_url = "http://domain.com/path/to/above.box"
    config.vm.box_url = "http://files.vagrantup.com/precise64.box"

    # Create a forwarded port mapping which allows access to a specific port
    # within the machine from a port on the host machine. In the example below,
    # accessing "localhost:8080" will access port 80 on the guest machine:
    config.vm.network :forwarded_port, guest: 80, host: 8080

    # Create a private network, which allows host-only access to the machine
    # using a specific IP:
    #config.vm.network :private_network, ip: "192.168.33.10"

    # Create a public network, which generally matched to bridged network.
    # Bridged networks make the machine appear as another physical device on
    # your network:
    config.vm.network :public_network

    # If true, then any SSH connections made will enable agent forwarding.
    # Default value: false
    #config.ssh.forward_agent = true

    # Share an additional folder to the guest VM. The first argument is
    # the path on the host to the actual folder. The second argument is
    # the path on the guest to mount the folder. And the optional third
    # argument is a set of non-required options:
    #config.vm.synced_folder "../data", "/vagrant_data"
    $out
    $freesurfer
    $ants
    $atlases

    # Provider-specific configuration so you can fine-tune various
2    # backing providers for Vagrant. These expose provider-specific options.
    # Example for VirtualBox:
    config.vm.provider :virtualbox do |vb|
      vb.customize ["modifyvm", :id, "--cpuexecutioncap", "$cpu"]
      vb.customize ["modifyvm", :id, "--memory", "$mem"]
      vb.customize ["modifyvm", :id, "--cpus", "$num"]
      vb.customize ["modifyvm", :id, "--ioapic", "on"]
    end

    config.vm.provision "shell", :privileged => false, inline: $$script

end
"""

#=============================================================================
# Script to include in Vagrant file
#=============================================================================
script = """$script = <<SCRIPT

# Install Anaconda Python distribution:
wget http://repo.continuum.io/miniconda/Miniconda-2.2.2-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
#echo "export PATH=$HOME/anaconda/bin:\$PATH" >> .bashrc
PATH=$HOME/anaconda/bin:$PATH

# Install nipype and nipype dependencies:
$HOME/anaconda/bin/conda install --yes pip
$HOME/anaconda/bin/conda install --yes numpy scipy nose traits networkx
$HOME/anaconda/bin/conda install --yes dateutil ipython-notebook matplotlib
$HOME/anaconda/bin/pip install nibabel --use-mirrors
$HOME/anaconda/bin/pip install https://github.com/RDFLib/rdflib/archive/master.zip
$HOME/anaconda/bin/pip install https://github.com/satra/prov/archive/enh/rdf.zip
$HOME/anaconda/bin/pip install https://github.com/nipy/nipype/archive/master.zip

# Install compiling utilities:
sudo apt-get install g++
sudo apt-get install cmake  # $HOME/anaconda/bin/conda install --yes cmake
sudo apt-get install make
sudo apt-get install git
sudo apt-get install xorg openbox

# Install VTK and set environment variables:
$HOME/anaconda/bin/conda install --yes vtk
export VTK_DIR=$HOME/anaconda/lib/vtk-5.10/

# Install Mindboggle:
#PREFIX_PATH=$HOME/Mindboggle
#$HOME/anaconda/bin/pip install --install-option="--prefix=$PREFIX_PATH" https://github.com/binarybottle/mindboggle/archive/master.zip
$HOME/anaconda/bin/pip install https://github.com/binarybottle/mindboggle/archive/master.zip

# Set Mindboggle environment variables:
MINDBOGGLE_TOOLS=/vagrant/mindboggle_tools/bin/
PATH=$MINDBOGGLE_TOOLS:$PATH
#export DYLD_LIBRARY_PATH=$HOME/anaconda/lib/vtk-5.10:${DYLD_LIBRARY_PATH}
export MINDBOGGLE_TOOLS

# Install FreeSurfer:
# http://surfer.nmr.mgh.harvard.edu/fswiki/Download
#wget -c ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/5.3.0/freesurfer-Linux-centos4_x86_64-stable-pub-v5.3.0.tar.gz
#mv freesurfer-Linux-centos4_x86_64-stable-pub-v5.3.0.tar.gz /usr/local
#cd /usr/local
#tar xzvf freesurfer-Linux-centos4_x86_64-stable-pub-v5.3.0.tar.gz
#cd $HOME

# Set FreeSurfer environment variables:
FREESURFER_HOME=/usr/local/freesurfer
SUBJECTS_DIR=${FREESURFER_HOME}/subjects
PATH=$SUBJECTS_DIR:$FREESURFER_HOME:$PATH
export FREESURFER_HOME SUBJECTS_DIR
#source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Install ANTs:
# http://brianavants.wordpress.com/2012/04/13/updated-ants-compile-instructions-april-12-2012/
#git clone git://github.com/stnava/ANTs.git
#mkdir antsbin
#cd antsbin
#cmake ../ANTs
#make # -j 4
#cp ../ANTs/Scripts/* bin/
#cd $HOME

# Set ANTs environment variables:
ANTSPATH=$HOME/antsbin/bin/
PATH=$ANTSPATH:$PATH
export ANTSPATH

echo "export $PATH" >> .bashrc

#source .bashrc

# Compile Mindboggle C++ code:
cd /vagrant/mindboggle_tools/bin/
$HOME/anaconda/bin/cmake /vagrant/mindboggle_tools
make
cd /

SCRIPT
"""

#=============================================================================
# Write FreeSurfer .license file
#=============================================================================
freesurfer_license = """arno@mindboggle.info
18192
 *Cr4e1z13elAY"""
license_file = os.path.join(os.environ['FREESURFER_HOME'], '.license')
f = open(license_file, 'w')
f.write(freesurfer_license)
f.close()

#=============================================================================
# Write Vagrantfile with substitutions, and run "vagrant up"
#=============================================================================
template = Template(vagrant_file)
new_vagrant_file = template.substitute(script=script,
                                       out=out,
                                       freesurfer=freesurfer,
                                       ants=ants,
                                       atlases=atlases,
                                       cpu=cpu,
                                       mem=mem,
                                       num=num)
f = open('Vagrantfile', 'w')
f.write(new_vagrant_file)
f.close()
