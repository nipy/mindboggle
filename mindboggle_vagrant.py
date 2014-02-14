#!/usr/bin/python
"""
Mindboggle virtual machine setup script.

This program installs Mindboggle and dependencies, such as
FreeSurfer and ANTs, in a VirtualBox virtual machine (vm)
using Vagrant software, which can package the vm.
You must have Vagrant (http://www.vagrantup.com) and
VirtualBox (http://www.virtualbox.org) installed,
as well as a good internet connection.

The Mindboggle software automates shape analysis of anatomical labels
and features extracted from human brain MR image data.
See http://mindboggle.info for up-to-date documentation.

To install ANTs and FreeSurfer in a local mounted directory::

    python mindboggle_vagrant.py --install
    vagrant up install

To configure mindboggle's local input and output directories
(type "python mindboggle_vagrant.py -h" for descriptions of the arguments)::

    python mindboggle_vagrant.py --dependencies XXX --ants XXX --out XXX

And to run mindboggle::

    vagrant up run


ADVANCED: To build the mindboggle.box Vagrant box from scratch::

    git clone https://github.com/binarybottle/mindboggle.git
    mv mindboggle_vagrant.py mindboggle/
    cd mindboggle
    python mindboggle_vagrant.py --install --build_from_scratch
    vagrant up build
    vagrant package --base build --output minboggle.box
    <Upload mindboggle.box to http://mindboggle.info/>


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
The Mindboggle Python package automates shape analysis of anatomical labels
and features extracted from human brain MR image data. http://mindboggle.info
has up-to-date documentation. This program installs Mindboggle and
dependencies, such as FreeSurfer and ANTs, in a VirtualBox virtual machine
using Vagrant software, which can package the virtual machine. You must have
Vagrant (http://www.vagrantup.com) and VirtualBox (http://www.virtualbox.org)
installed, as well as a good internet connection. Ex: python
mindboggle_vagrant.py
    --dependencies /data/mindboggle_dependencies
    --ants /data/antsCorticalThickness
    --num 6 --mem 4096 --cpu 75 """,
                     formatter_class = lambda prog:
                     argparse.HelpFormatter(prog, max_help_position=40))

parser.add_argument("-o", "--out",
                    help='Output directory (default: $HOME/mindboggled)',
                    default=os.path.join(os.environ['HOME'],
                                         'mindboggled'), metavar='STR')
parser.add_argument("--freesurfer",
                    help=("FreeSurfer subjects directory (default: "
                          "$SUBJECTS_DIR)"),
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
parser.add_argument("--dependencies",
                    help=("Local directory with the ANTs and FreeSurfer "
                          "software (same as the --install argument, "
                          "but --install is only set when installing)"),
                    metavar='STR')
parser.add_argument("--install",
                    help=("Install ANTs and FreeSurfer software in a local "
                          "mounted directory (--install is used only once "
                          "before running mindboggle for the first time)"),
                    action='store_true')
parser.add_argument("--build_from_scratch",
                    help=("Build mindboggle.box Vagrant box from scratch, "
                          "which contains everything for mindboggle except "
                          "ANTs and FreeSurfer (ADVANCED)"),
                    action='store_true')

args = parser.parse_args()

print(args)

home = "/home/vagrant/"
if not os.path.exists(args.out):
    os.mkdir(args.out)
out_string = 'config.vm.synced_folder "{0}", ' \
             '"{1}mindboggled"'.format(args.out, home)
if args.freesurfer:
    freesurfer_string = 'config.vm.synced_folder "{0}", ' \
        '"{1}freesurfer_subjects"'.format(args.freesurfer, home)
else:
    freesurfer_string = 'config.vm.synced_folder "{0}", ' \
        '"{1}freesurfer_subjects"'.format(os.environ['SUBJECTS_DIR'], home)
if args.ants:
    ants_string = 'config.vm.synced_folder "{0}", ' \
                  '"{1}ants_subjects"'.format(args.ants, home)
else:
    ants_string = ''
if args.atlases:
    atlases_string = 'config.vm.synced_folder "{0}", ' \
                     '"{1}atlases"'.format(args.atlases, home)
else:
    atlases_string = ''
dependencies = args.dependencies
if dependencies:
    dependencies_string = 'config.vm.synced_folder "{0}", ' \
                          '"{1}dependencies"'.format(dependencies, home)
else:
    dependencies_string = ''

#=============================================================================
# Create Vagrantfile
#=============================================================================
vagrant_file = """# Vagrant file (http://www.vagrantup.com)
# Vagrantfile API/syntax version:
VAGRANTFILE_API_VERSION = "2"

$script

#=============================================================================
# Configure Vagrantfile
#=============================================================================
# Configure Vagrant:
Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|

    #-------------------------------------------------------------------------
    # Build the mindboggle.box Vagrant base box -- ONLY FOR THE INITIAL BUILD!
    #-------------------------------------------------------------------------
    config.vm.define :build do |build_config|

        # Build from an existing Vagrant virtual box:
        build_config.vm.box = "precise64"
        build_config.vm.box_url = "http://files.vagrantup.com/precise64.box"

        # Create a forwarded port mapping which allows access to a specific port
        # within the machine from a port on the host machine. In the example below,
        # accessing "localhost:8080" will access port 80 on the guest machine:
        build_config.vm.network :forwarded_port, guest: 80, host: 8080
        build_config.vm.network :public_network
        #build_config.vm.network :private_network, ip: "192.168.33.10"

        build_config.vm.provider :virtualbox do |vb|
          vb.customize ["modifyvm", :id, "--cpuexecutioncap", "$cpu"]
          vb.customize ["modifyvm", :id, "--memory", "$mem"]
          vb.customize ["modifyvm", :id, "--cpus", "$num"]
          vb.customize ["modifyvm", :id, "--ioapic", "on"]
        end

        build_config.vm.provision "shell", :privileged => false, inline: $$script
    end

    #-------------------------------------------------------------------------
    # Configure the mindboggle.box Vagrant base box to install dependencies
    #-------------------------------------------------------------------------
    config.vm.define :install do |install_config|

        # Build from the mindboggle.box Vagrant virtual box
        # (VirtualBox caches it so you only download once):
        install_config.vm.box = "mindboggle_install"
        install_config.vm.box_url = "http://mindboggle.info/mindboggle_20140214.box"
        install_config.vm.hostname = "boggler"

        # Create a forwarded port mapping which allows access to a specific port
        # within the machine from a port on the host machine. In the example below,
        # accessing "localhost:8080" will access port 80 on the guest machine:
        install_config.vm.network :forwarded_port, guest: 80, host: 8000 #8080
        install_config.vm.network :public_network
        #install_config.vm.network :private_network, ip: "192.168.33.10"

        install_config.vm.provider :virtualbox do |vb|
          vb.customize ["modifyvm", :id, "--cpuexecutioncap", "$cpu"]
          vb.customize ["modifyvm", :id, "--memory", "$mem"]
          vb.customize ["modifyvm", :id, "--cpus", "$num"]
          vb.customize ["modifyvm", :id, "--ioapic", "on"]
        end

        install_config.vm.provision "shell", :privileged => false, inline: $$script
    end

    #-------------------------------------------------------------------------
    # Configure the mindboggle.box Vagrant base box to mount directories
    #-------------------------------------------------------------------------
    config.vm.define :run do |run_config|

        # Build from the mindboggle.box Vagrant virtual box
        # (VirtualBox caches it so you only download once):
        run_config.vm.box = "mindboggle"
        run_config.vm.box_url = "http://mindboggle.info/mindboggle_20120214.box"
        run_config.vm.hostname = "boggler"

        # Share an additional folder to the guest VM. The first argument is
        # the path on the host to the actual folder. The second argument is
        # the path on the guest to mount the folder.
        #run_config.vm.synced_folder "../data", "/vagrant_data"
        $dependencies_string
        $freesurfer_string
        $ants_string
        $atlases_string
        $out_string

        # Create a forwarded port mapping which allows access to a specific port
        # within the machine from a port on the host machine. In the example below,
        # accessing "localhost:8080" will access port 80 on the guest machine:
        run_config.vm.network :forwarded_port, guest: 80, host: 8080
        run_config.vm.network :public_network
        #run_config.vm.network :private_network, ip: "192.168.33.10"

        run_config.vm.provider :virtualbox do |vb|
          vb.customize ["modifyvm", :id, "--cpuexecutioncap", "$cpu"]
          vb.customize ["modifyvm", :id, "--memory", "$mem"]
          vb.customize ["modifyvm", :id, "--cpus", "$num"]
          vb.customize ["modifyvm", :id, "--ioapic", "on"]
        end

        run_config.vm.provision "shell", :privileged => false, inline: $$script
    end
end
"""

#=============================================================================
# Build script to include in Vagrant file
#=============================================================================
script = ''
if args.install and args.build_from_scratch:
    script = """$script = <<SCRIPT

#-----------------------------------------------------------------------------
# Install Anaconda Python distribution:
#-----------------------------------------------------------------------------
wget http://repo.continuum.io/miniconda/Miniconda-2.2.2-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
export PATH=/home/vagrant/anaconda/bin:$PATH

#-----------------------------------------------------------------------------
# Install nipype and nipype dependencies:
#-----------------------------------------------------------------------------
conda install --yes pip cmake
conda install --yes numpy scipy nose traits networkx
conda install --yes dateutil ipython-notebook matplotlib
pip install nibabel --use-mirrors
pip install https://github.com/RDFLib/rdflib/archive/master.zip
pip install https://github.com/satra/prov/archive/enh/rdf.zip
pip install https://github.com/nipy/nipype/archive/master.zip

#-----------------------------------------------------------------------------
# Install compiling utilities:
#-----------------------------------------------------------------------------
sudo apt-get update
sudo apt-get install -y g++
sudo apt-get install -y make
sudo apt-get install -y git
sudo apt-get install -y xorg openbox

#-----------------------------------------------------------------------------
# Install VTK:
#-----------------------------------------------------------------------------
conda install --yes vtk
#export VTK_DIR=/home/vagrant/anaconda/lib/vtk-5.10/  # Needed for ANTs
echo "export VTK_DIR=/home/vagrant/anaconda/lib/vtk-5.10/" >> .bashrc

#-----------------------------------------------------------------------------
# Install Mindboggle:
#-----------------------------------------------------------------------------
#git clone https://github.com/binarybottle/mindboggle.git
pip install https://github.com/binarybottle/mindboggle/archive/master.zip
echo "export MINDBOGGLE_TOOLS=/vagrant/mindboggle_tools/bin/" >> .bashrc
echo "export PATH=\$MINDBOGGLE_TOOLS:\$ANTSPATH:\$FREESURFER:\$PATH" >> .bashrc
cd /vagrant/mindboggle_tools/bin/
cmake /vagrant/mindboggle_tools
make
cd /home/vagrant

source .bashrc

SCRIPT
"""

#=============================================================================
# Install ANTs and FreeSurfer in mounted directory
#=============================================================================
if args.install and not args.build_from_scratch:
    script = """$script = <<SCRIPT

#-----------------------------------------------------------------------------
# Install FreeSurfer in mounted directory:
#-----------------------------------------------------------------------------
# http://surfer.nmr.mgh.harvard.edu/fswiki/Download
wget -c ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/5.3.0/freesurfer-Linux-centos4_x86_64-stable-pub-v5.3.0.tar.gz
sudo mv freesurfer-Linux-centos4_x86_64-stable-pub-v5.3.0.tar.gz $dependencies
cd $dependencies
sudo tar xzvf freesurfer-Linux-centos4_x86_64-stable-pub-v5.3.0.tar.gz
rm freesurfer-Linux-centos4_x86_64-stable-pub-v5.3.0.tar.gz
cd /home/vagrant
echo "export FREESURFER_HOME=$dependencies/freesurfer" >> .bashrc
echo "export SUBJECTS_DIR=\$FREESURFER_HOME/subjects" >> .bashrc
echo "source \$FREESURFER_HOME/SetUpFreeSurfer.sh" >> .bashrc

#-----------------------------------------------------------------------------
# Install ANTs in mounted directory:
#-----------------------------------------------------------------------------
# http://brianavants.wordpress.com/2012/04/13/
#        updated-ants-compile-instructions-april-12-2012/
cd $dependencies
git clone git://github.com/stnava/ANTs.git
mkdir antsbin
cd antsbin
cmake ../ANTs
make #-j 4
cp ../ANTs/Scripts/* bin/
cd /home/vagrant
echo "export ANTSPATH=$dependencies/antsbin/bin/" >> .bashrc

source .bashrc

SCRIPT
"""

if args.install and not args.build_from_scratch:

    freesurfer_license = """arno@mindboggle.info
    18192
     *Cr4e1z13elAY"""
    license_file = os.path.join(os.environ['FREESURFER_HOME'], '.license')
    f = open(license_file, 'w')
    f.write(freesurfer_license)
    f.close()

#=============================================================================
# Write Vagrantfile with substitutions
#=============================================================================
template = Template(vagrant_file)
new_vagrant_file = template.substitute(script=script,
                                       out_string=out_string,
                                       freesurfer_string=freesurfer_string,
                                       ants_string=ants_string,
                                       atlases_string=atlases_string,
                                       cpu=args.cpu,
                                       mem=args.mem,
                                       num=args.num,
                                       dependencies=dependencies,
                                       dependencies_string=dependencies_string)
f = open('Vagrantfile', 'w')
f.write(new_vagrant_file)
f.close()
