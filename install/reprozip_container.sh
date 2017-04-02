# Reduce the size of the mindboggle docker container.
#-----------------------------------------------------------------------------
# First build the docker image:
# docker build -t mindboggle -f Dockerfile.mindboggle.complete .
#
# Then run the docker container:
# docker run --rm -ti -v /homedir:/home/jovyan/work -p 8888:8888 -p 5000:5000 mindboggle bash
#
# In this script, use reprozip to reduce the mindboggle docker container
# to just the necessary components.
#
# (c) 2017 by Arno Klein <arno@mindboggle.info> (CC-BY license)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Retain only those portions of the docker container that are required
#-----------------------------------------------------------------------------
DATA=/data
T1=$DATA/mindboggle_input_example/raw/T1.nii.gz
T2=$DATA/mindboggle_input_example/raw/T2.nii.gz
SUBJECT=arno
TEMPLATE=$DATA/OASIS-30_Atropos_template
FREESURFER_SUBJECT=$DATA/mindboggle_input_example/freesurfer/subjects/arno
ANTS_SUBJECT=$DATA/mindboggle_input_example/ants/subjects/arno
MINDBOGGLECACHE=$DATA/mindboggle_cache
MINDBOGGLING=$DATA/mindboggling
MINDBOGGLED=$DATA/mindboggled

#-----------------------------------------------------------------------------
# Retain only those portions of the docker container that are required
#-----------------------------------------------------------------------------
# reprozip trace -d mindboggle \
#     recon-all -all -i $T1 -i $T2 -s $SUBJECT && \
#     antsCorticalThickness.sh -d 3 -a $T1 -o ants \
#         -e $TEMPLATE/T_template0.nii.gz \
#         -t $TEMPLATE/T_template0_BrainCerebellum.nii.gz \
#         -m $TEMPLATE/T_template0_BrainCerebellumProbabilityMask.nii.gz \
#         -f $TEMPLATE/T_template0_BrainCerebellumExtractionMask.nii.gz \
mindboggle $FREESURFER_SUBJECT --cpus 8 \
         --ants $ANTS_SUBJECT/antsBrainSegmentation.nii.gz \
         --out $MINDBOGGLED \
         --working $MINDBOGGLING \
         --cache $MINDBOGGLECACHE \
         --no_volumes \
         --no_surfaces \
         --no_labels \
         --no_shapes \
         --no_sulci \
         --no_points \
         --no_moments \
         --no_spectra \
         --no_thickness \
         --fundi \
         --moments 10 \
         --spectra 10 \
         --graph hier \
         --roygbiv
#        --help --version \
#        --my_atlas --my_atlases --my_graywhite --my_transform \
#        --plugin --plugin_args


# reprozip pack -d mindboggle-trace mindboggle_20170402.rpz