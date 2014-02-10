# Vagrant file (http://www.vagrantup.com)
# (after https://github.com/nipy/nipype/blob/master/Vagrantfile)
# Vagrantfile API/syntax version:
VAGRANTFILE_API_VERSION = "2"

$script = <<SCRIPT

# Install Anaconda Python distribution:
wget http://repo.continuum.io/miniconda/Miniconda-2.2.2-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b
echo "export PATH=$HOME/anaconda/bin:\$PATH" >> .bashrc

# Install nipype and nipype dependencies:
$HOME/anaconda/bin/conda install --yes pip numpy scipy nose traits networkx
$HOME/anaconda/bin/conda install --yes dateutil ipython-notebook matplotlib
$HOME/anaconda/bin/pip install nibabel --use-mirrors
$HOME/anaconda/bin/pip install https://github.com/RDFLib/rdflib/archive/master.zip
$HOME/anaconda/bin/pip install https://github.com/satra/prov/archive/enh/rdf.zip
$HOME/anaconda/bin/pip install https://github.com/nipy/nipype/archive/master.zip

# Install compiling utilities:
$HOME/anaconda/bin/conda install --yes cmake
sudo apt-get install make

# Install Mindboggle and dependencies (besides FreeSurfer and ANTs):
$HOME/anaconda/bin/conda install --yes vtk
#PREFIX_PATH=$HOME/Mindboggle
#$HOME/anaconda/bin/pip install --install-option="--prefix=$PREFIX_PATH" https://github.com/binarybottle/mindboggle/archive/master.zip
$HOME/anaconda/bin/pip install https://github.com/binarybottle/mindboggle/archive/master.zip

# Set Mindboggle environment variables:
MINDBOGGLE_TOOLS=/vagrant/mindboggle_tools/bin/
echo "export $MINDBOGGLE_TOOLS:\$PATH" >> .bashrc
#export DYLD_LIBRARY_PATH=$HOME/anaconda/lib/vtk-5.10:${DYLD_LIBRARY_PATH}

# Set FreeSurfer environment variables:
FREESURFER_HOME=$HOME/freesurfer
SUBJECTS_DIR=$HOME/freesurfer/subjects
echo "export $FREESURFER_HOME:\$PATH" >> .bashrc
echo "export $SUBJECTS_DIR:\$PATH" >> .bashrc
#source $FREESURFER_HOME/SetUpFreeSurfer.sh

# Set ANTs environment variables:
ANTSPATH=$HOME/antsbin/bin/
echo "export $ANTSPATH:\$PATH" >> .bashrc

# Compile Mindboggle C++ code:
cd /vagrant/mindboggle_tools/bin/
$HOME/anaconda/bin/cmake /vagrant/mindboggle_tools
make
cd /

SCRIPT


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
    #config.vm.network :public_network

    # If true, then any SSH connections made will enable agent forwarding.
    # Default value: false
    #config.ssh.forward_agent = true

    # Share an additional folder to the guest VM. The first argument is
    # the path on the host to the actual folder. The second argument is
    # the path on the guest to mount the folder. And the optional third
    # argument is a set of non-required options:
    #config.vm.synced_folder "../data", "/vagrant_data"
    config.vm.synced_folder "/appsdir/freesurfer/subjects", "/freesurfer_subjects"
    config.vm.synced_folder "/data/antsCorticalThickness", "/ants_subjects"
    

    # Provider-specific configuration so you can fine-tune various
    # backing providers for Vagrant. These expose provider-specific options.
    # Example for VirtualBox:
    config.vm.provider :virtualbox do |vb|
      vb.customize ["modifyvm", :id, "--cpuexecutioncap", "80"]
      vb.customize ["modifyvm", :id, "--memory", "8192"]
      vb.customize ["modifyvm", :id, "--cpus", "4"]
      vb.customize ["modifyvm", :id, "--ioapic", "on"]
    end

    config.vm.provision "shell", :privileged => false, inline: $script

end
