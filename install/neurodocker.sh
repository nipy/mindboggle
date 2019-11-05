#!/bin/bash

###############################################################################
# Generate a Dockerfile and Singularity recipe for building a Mindboggle container
# (https://mindboggle.info).
# The Dockerfile and/or Singularity recipe installs most of Mindboggle's dependencies,
# including preprocessing software packages FreeSurfer and ANTs,
# and visualization software roygbiv:
#   - Set up neurodebian and miniconda and install conda packages
#   - Install FreeSurfer and ANTS.
#   - Install Mindboggle templates from the Open Science Framework
#   - Install Mindboggle
#   - Install roygbiv for output visualization
#
# Steps to build, upload, and deploy the Mindboggle docker and/or singularity image:
#
# 1. Create or update the Dockerfile and Singuarity recipe:
# bash neurodocker.sh
#
# 2. Build the docker image:
# docker build -t mindboggle -f Dockerfile .
#
#    and/or singularity image:
# singularity build mindboggle.simg Singularity 
#
# 3. Push to Docker hub:
# (https://docs.docker.com/docker-cloud/builds/push-images/)
# export DOCKER_ID_USER="nipy"
# docker login
# docker tag mindboggle nipy/mindboggle  # See: https://docs.docker.com/engine/reference/commandline/tag/
# docker push nipy/mindboggle
#
# 4. Pull from Docker hub (or use the original):
# docker pull nipy/mindboggle
#
# In the following, the Docker container can be the original (mindboggle)
# or the pulled version (nipy/mindboggle), and is given access to /Users/arno
# on the host machine.
#
# 5. Enter the bash shell of the Docker container, and add port mappings:
# docker run --rm -ti -v /Users/arno:/home/jovyan/work -p 8888:8888 -p 5000:5000 nipy/mindboggle bash
#
# 6. Run the Docker container as an executable (variables set for clarity):
# HOST=/Users/binarybottle # path on host to access input and output
# DOCK=/home/jovyan/work # path to HOST from Docker container
# IMAGE=$DOCK/example_mri_data/T1.nii.gz # input image (from container)
# ID=arno # ID for brain image
# OUT=$DOCK/mindboggle123_output # '--out $OUT' is OPTIONAL
# docker run --rm -ti -v $HOST:/home/jovyan/work nipy/mindboggle $IMAGE --id $ID --out $OUT
#
###############################################################################

image="kaczmarj/neurodocker:master@sha256:936401fe8f677e0d294f688f352cbb643c9693f8de371475de1d593650e42a66"

# Generate a dockerfile for building a mindboggle container
docker run --rm ${image} generate docker \
  --base neurodebian:stretch \
  --pkg-manager apt \
  --install graphviz tree git-annex-standalone vim \
    emacs-nox nano less ncdu tig sed build-essential \
    libsm-dev libx11-dev libxt-dev libxext-dev libglu1-mesa \
  --run 'ln -s /usr/lib/x86_64-linux-gnu /usr/lib64' \
  --miniconda \
    conda_install="python=3.6 pip jupyter cmake nipype>=1.1.4 mesalib vtk=8.2.0=py36ha8e561a_201 pandas
      matplotlib colormath nilearn tbb-devel nose etelemetry" \
    pip_install="datalad[full] duecredit" \
    create_env="mb" \
    activate=true \
  --workdir /opt \
  --run 'mkdir -p /opt/data && cd /opt/data && \
    curl -sSL https://osf.io/download/rh9km/?revision=2 -o templates.zip && \
    unzip templates.zip && \
    rm -f /opt/data/templates.zip && \
    curl -sSL https://files.osf.io/v1/resources/hvc52/providers/osfstorage/57c1a8f06c613b01f98d68a9/?zip= -o OASIS-TRT-20_brains.zip && \
    mkdir OASIS-TRT-20_brains && \
    cd OASIS-TRT-20_brains && \
    unzip ../OASIS-TRT-20_brains.zip && \
    cd .. && \
    rm OASIS-TRT-20_brains.zip && \
    curl -sSL https://files.osf.io/v1/resources/zevma/providers/osfstorage/5783dfcab83f6901f963735c/?zip= -o cmalabels.zip && \
    mkdir OASIS-TRT-20_DKT31_CMA_labels_v2 && \
    cd OASIS-TRT-20_DKT31_CMA_labels_v2 && \
    unzip ../cmalabels.zip && \
    cd .. && \
    rm cmalabels.zip && \
    curl -sSL https://osf.io/download/d2cmy/ -o OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_OASIS-30_v2.nii.gz && \
    rm -rf __MACOSX' \
  --run-bash 'source /opt/miniconda-latest/etc/profile.d/conda.sh && \
        conda activate mb && \
        git clone https://github.com/nipy/mindboggle.git && \
        cd /opt/mindboggle && \
        python setup.py install && \
        mkdir /opt/vtk_cpp_tools && \
        cd /opt/vtk_cpp_tools && \
        cmake /opt/mindboggle/vtk_cpp_tools && \
        make' \
  --env vtk_cpp_tools=/opt/vtk_cpp_tools \
  --run-bash 'source /opt/miniconda-latest/etc/profile.d/conda.sh && \
        conda activate mb && \
        conda install -y flask && \
        git clone https://github.com/akeshavan/roygbiv && \
        cd /opt/roygbiv && \
        git checkout fbbf31c29952d0ea22ed05d98e0a5a7e7d0827f9 && \
        python setup.py install && \
        cd /opt && \
        rm -rf /opt/roygbiv' \
  --run 'mkdir -p /.jupyter && echo c.NotebookApp.ip = \"0.0.0.0\" > /.jupyter/jupyter_notebook_config.py' \
  --ants version=b43df4bfc8 method=source cmake_opts='-DBUILD_SHARED_LIBS=ON' make_opts='-j 4' \
  --freesurfer version=6.0.0-min \
  --run 'curl -sSL https://osf.io/download/n3ud2/?revision=1 -o /opt/freesurfer-6.0.0-min/license.txt' \
> Dockerfile


# Generate a singularity recipe for building a mindboggle container
docker run --rm ${image} generate singularity \
  --base neurodebian:stretch \
  --pkg-manager apt \
  --install graphviz tree git-annex-standalone vim \
    emacs-nox nano less ncdu tig sed build-essential \
    libsm-dev libx11-dev libxt-dev libxext-dev libglu1-mesa \
  --run 'ln -s /usr/lib/x86_64-linux-gnu /usr/lib64' \
  --miniconda \
    conda_install="python=3.6 pip jupyter cmake mesalib vtk=8.2.0=py36ha8e561a_201 pandas
      matplotlib colormath nipype>=1.1.4 nilearn tbb-devel nose etelemetry" \
    pip_install="datalad[full] duecredit" \
    create_env="mb" \
    activate=true \
  --workdir /opt \
  --run 'mkdir -p /opt/data && cd /opt/data && \
    curl -sSL https://osf.io/download/rh9km/?revision=2 -o templates.zip && \
    unzip templates.zip && \
    rm -f /opt/data/templates.zip && \
    curl -sSL https://files.osf.io/v1/resources/hvc52/providers/osfstorage/57c1a8f06c613b01f98d68a9/?zip= -o OASIS-TRT-20_brains.zip && \
    mkdir OASIS-TRT-20_brains && \
    cd OASIS-TRT-20_brains && \
    unzip ../OASIS-TRT-20_brains.zip && \
    cd .. && \
    rm OASIS-TRT-20_brains.zip && \
    curl -sSL https://files.osf.io/v1/resources/zevma/providers/osfstorage/5783dfcab83f6901f963735c/?zip= -o cmalabels.zip && \
    mkdir OASIS-TRT-20_DKT31_CMA_labels_v2 && \
    cd OASIS-TRT-20_DKT31_CMA_labels_v2 && \
    unzip ../cmalabels.zip && \
    cd .. && \
    rm cmalabels.zip && \
    curl -sSL https://osf.io/download/d2cmy/ -o OASIS-TRT-20_jointfusion_DKT31_CMA_labels_in_OASIS-30_v2.nii.gz && \
    rm -rf __MACOSX' \
  --run-bash 'source /opt/miniconda-latest/etc/profile.d/conda.sh && \
        conda activate mb && \
        git clone https://github.com/nipy/mindboggle.git && \
        cd /opt/mindboggle && \
        python setup.py install && \
        mkdir /opt/vtk_cpp_tools && \
        cd /opt/vtk_cpp_tools && \
        cmake /opt/mindboggle/vtk_cpp_tools && \
        make' \
  --env vtk_cpp_tools=/opt/vtk_cpp_tools \
  --run-bash 'source /opt/miniconda-latest/etc/profile.d/conda.sh && \
        conda activate mb && \
        conda install -y flask && \
        git clone https://github.com/akeshavan/roygbiv && \
        cd /opt/roygbiv && \
        git checkout fbbf31c29952d0ea22ed05d98e0a5a7e7d0827f9 && \
        python setup.py install && \
        cd /opt && \
        rm -rf /opt/roygbiv' \
  --run 'mkdir -p /.jupyter && echo c.NotebookApp.ip = \"0.0.0.0\" > /.jupyter/jupyter_notebook_config.py' \
  --ants version=b43df4bfc8 method=source cmake_opts='-DBUILD_SHARED_LIBS=ON' make_opts='-j 4' \
  --freesurfer version=6.0.0-min \
  --run 'curl -sSL https://osf.io/download/n3ud2/?revision=1 -o /opt/freesurfer-6.0.0-min/license.txt' \
> Singularity
