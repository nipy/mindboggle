# Reduce the size of the mindboggle docker container.
# NOTE: The reprozip command is still under construction!
#-----------------------------------------------------------------------------
# First build the docker image:
# docker build -t mindboggle -f Dockerfile.mindboggle.complete .
#
# Then run the docker container:
# docker run --rm -ti -v /Users/arno:/home/jovyan/work -p 8888:8888 -p 5000:5000 mindboggle bash
#
# In this script, use reprozip to reduce the mindboggle docker container
# to just the necessary components.
#
# (c) 2017 by Arno Klein <arno@mindboggle.info> (CC-BY license)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Retain only those portions of the docker container that are required
#-----------------------------------------------------------------------------
DATA=/home/jovyan/work/Data
T1=$DATA/example_mri_data/T1.nii.gz
T2=$DATA/example_mri_data/T2.nii.gz
SUBJECT=arno
TEMPLATE=$DATA/OASIS-30_Atropos_template
FREESURFER_SUBJECTS=$DATA/freesurfer_subjects
FREESURFER_SUBJECT=$FREESURFER_SUBJECTS/$SUBJECT
ANTS_SUBJECT=$DATA/ants_subjects/$SUBJECT
#FREESURFER_SUBJECT=$DATA/mindboggle_input_example/freesurfer/subjects/arno
#ANTS_SUBJECT=$DATA/mindboggle_input_example/ants/subjects/arno
MINDBOGGLED=$DATA/mindboggled

#-----------------------------------------------------------------------------
# Retain only those portions of the docker container that are required
#-----------------------------------------------------------------------------
mkdir $ANTS_SUBJECT
reprozip trace -d mindboggle \
    recon-all -all -i $T1 -T2 $T2 -s $SUBJECT -sd $FREESURFER_SUBJECTS -cw256 && \
    antsCorticalThickness.sh -d 3 -a $T1 -o $ANTS_SUBJECT/ants \
        -e $TEMPLATE/T_template0.nii.gz \
        -t $TEMPLATE/T_template0_BrainCerebellum.nii.gz \
        -m $TEMPLATE/T_template0_BrainCerebellumProbabilityMask.nii.gz \
        -f $TEMPLATE/T_template0_BrainCerebellumExtractionMask.nii.gz \
        -p $TEMPLATE/Priors2/priors%d.nii.gz && \
    mindboggle $FREESURFER_SUBJECT --cpus 8 \
         --ants $ANTS_SUBJECT/antsBrainSegmentation.nii.gz \
         --out $MINDBOGGLED \
         --fundi \
         --moments 10 \
         --spectra 10 \
         --graph hier \
         --roygbiv
         #--no_volumes \
         #--no_surfaces \
         #--no_labels \
         #--no_shapes \
         #--no_sulci \
         #--no_points \
         #--no_moments \
         #--no_spectra \
         #--no_thickness \
         #--help --version \
         #--my_atlas --my_atlases --my_graywhite --my_transform \
         #--plugin --plugin_args \

#reprozip pack -d mindboggle-trace mindboggle.rpz
