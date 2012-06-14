#!/bin/bash
#
# Code for extracting a brain from a brain image, using SyN and Atropos.
#
# Usage:   "bash extract_brain.sh INPUT_FILE OUT_PATH"
# Example: "bash extract_brain.sh ./subjects/S01.nii.gz ./output/S01"
#
# NOTE: Make sure to set the ANTSPATH environment variable
#       either locally or on the cluster:
#         export ANTSPATH="/Users/arno/ANTS-build/"
#         export PATH="${PATH}:${ANTSPATH}"
#
#       Set fix_headers=1 below only if subject and template
#       headers do not correspond.
#
# Output files:
# mask_template.nii.gz   # mask1. mask: deformed template priors
# mask_atropos.nii.gz    # mask2. mask: postprocessed Atropos output
# mask.nii.gz            # mask3. mask1 x mask2
# mask2.nii.gz           # mask4. mask: second Atropos output (5 below)
# head_n3.nii.gz         # 1. original after n3 bias field correction
# brain_masked.nii.gz    # 2. original x mask1
# brain.nii.gz           # 3. original x mask3 (final)
# brain_n3.nii.gz        # 4.        1 x mask3 (final+n3)
# brain2.nii.gz          # 5. original x mask4 (final2)
# brain2_n3.nii.gz       # 6.        1 x mask4 (final2 n3)
# brain_segmented.nii.gz # 7.        4 after Atropos segmentation
#
# Steps of the algorithm:
# 0. Run N3 bias field correction.
# 1. Deform a template image to a subject image that we want to extract,
#    and apply the transform to template segmentation class images.
# 2. Create a mask about the boundary of the deformed template.
# 3. Run Atropos with the masked, transformed segmented images as priors.
# 4. Extract the largest components of class 1, 2, and 3 to extract the brain.
# 5. Mask the resulting brain with a transformed template brain mask.
#
# Note: Atropos should be used in a data-specific manner:
# Build a local template, model relevant brain tissues, and apply to new data.
# Ideally, one would localize errors using training data and refine the template.
#
# (c) Brian Avants, arno klein (2010) arno@mindboggle.info

####################
# SET STAGES TO RUN:
####################
run_n3=1
run_register=1
run_transform=1
run_mask=1
run_atropos=1

run_postprocess1=0
skip_postprocess1=1
run_postprocess2=1

fix_headers=0  # set fix_headers=0 unless there is a gross misregistration in ITK-Snap

###################################
# EDIT INFORMATION BELOW WITH CARE:
###################################
#export ANTSPATH="/Users/arno/ANTS-build/"
#export ANTSPATH="/data/BI/Toolbox/software/ants_mac_svn460/bin/"
#export PATH="${PATH}:${ANTSPATH}"
ITERS='30x90x20'
FORMAT='.nii.gz' # format of the template and output files (NIfTI: .nii.gz)
DIM=3

# Sensitive parameters:
#   questions:
#   1. ROI width setting?
#   2. PRIORWT = 0, 0.1, 0.5, 0.9, 0.95?
#   3. Incorporate distance variables? What should the variance be (value of $dist)?
CLUSTERSIZE=250    # reject isolated objects of this size or less in classes 1 2 3
NTHBIGGESTCLUST=10 # keep the Nth biggest clusters for classes 1 2 3
MRF=0  #0.25  # we need to optimize this term ...
PRIORWT=0.9
DIST=5,1
DILATE=10 # increasing this includes more tissue exterior to brain
ERODE=10  # increasing includes more tissue interior to brain

# Arguments:
INPUT_FILE=$1
OUT_PATH=$2
OUT=${OUT_PATH}/

# Template files in: T1, segmented, and partial-volum/smoothed individual classes:
TEMPLATEname='template18'
NCLASSES=8
#TEMPLATEname='B2template'
#NCLASSES=11
CLASSname=${TEMPLATEname}_class  # + "0[1,2,3...,#classes]"
LABELSname=${TEMPLATEname}_${NCLASSES}classes
TEMPLATE_PATH=$(pwd)/template_atropos/
TEMPLATE_IN=${TEMPLATE_PATH}${TEMPLATEname}${FORMAT}
CLASS_IN=${TEMPLATE_PATH}${CLASSname}  # + "0[1,2,3...,#classes]"
LABELS_IN=${TEMPLATE_PATH}${LABELSname}${FORMAT}

# Template files out:
TEMPLATE=${OUT}template${FORMAT}
CLASS=${OUT}class  # + "0[1,2,3...,#classes]"
LABELS=${OUT}labels${FORMAT}

# Output files:
MASK_TEMPLATE=${OUT}mask_template${FORMAT}  # mask1. mask: deformed template priors
MASK_ATROPOS=${OUT}mask_atropos${FORMAT}    # mask2. mask: postprocessed Atropos output
OUT_MASK=${OUT}mask${FORMAT}                # mask3. mask1 x mask2
OUT_MASK2=${OUT}mask2${FORMAT}              # mask4. mask: second Atropos output (5 below)
OUT_HEAD_N3=${OUT}head_n3${FORMAT}          # 1. original after n3 bias field correction
OUT_MASKED=${OUT}brain_masked${FORMAT}      # 2. original x mask1
OUT_BRAIN=${OUT}brain${FORMAT}              # 3. original x mask3 (final)
OUT_BRAIN_N3=${OUT}brain_n3${FORMAT}        # 4.        1 x mask3 (final n3)
# Post-postprocess output:
OUT_SEG=${OUT}brain_segmented${FORMAT}      # 5.        4 after Atropos segmentation
OUT_BRAIN2=${OUT}brain2${FORMAT}            # 6. original x mask4 (final2)
OUT_BRAIN2_N3=${OUT}brain2_n3${FORMAT}      # 7.        1 x mask4 (final2 n3)

# Intermediate files (classes, masks):
TEMP=${OUT}temp${FORMAT}
PRIORname=${OUT}prior  # + "0[1,2,3...,#classes]"
PRIOR=${PRIORname}%02d${FORMAT}
POST=${OUT}post%02d${FORMAT}
MASKD=${OUT}maskD${FORMAT}  # a larger mask containing brain and some exterior tissues
MASKE=${OUT}maskE${FORMAT}  # a smaller mask eroded within brain
MASKB=${OUT}maskB${FORMAT}  # the mask of brain-boundary
SEG=${OUT}seg${FORMAT}      # Atropos segmentation
SEG2=${OUT}seg2${FORMAT}    # Atropos segmentation #2 (post-postprocessing)
GM=${OUT}GM${FORMAT}        # Segmented grey matter   (post-postprocessing)
WM=${OUT}WM${FORMAT}        # Segmented white matter  (post-postprocessing)
CSF=${OUT}CSF${FORMAT}      # Segmented CSF           (post-postprocessing)

########
# Begin:
########

# Make output directory:
if [ -d $OUT ]; then
  echo "$OUT already exists."
else
  mkdir $OUT
fi

# Smooth the priors (>3) to accommodate individual residual variability (done once):
run_smooth_priors=0
if [ $run_smooth_priors -eq 1 ]; then
  if [ $NCLASSES -eq 11 ]; then
    MINCLASS=1
  elif [ $NCLASSES -eq 8 ]; then
    MINCLASS=4
  fi
  x=$MINCLASS
  let "d=$NCLASSES+1"
  while [ $x -lt $d ]; do
    ${ANTSPATH}ThresholdImage $DIM ${LABELS_IN} ${CLASS_IN}0${x}${FORMAT} ${x} ${x}
    ${ANTSPATH}SmoothImage    $DIM ${CLASS_IN}0${x}${FORMAT} 1 ${CLASS_IN}0${x}${FORMAT}
    ${ANTSPATH}ImageMath      $DIM ${CLASS_IN}0${x}${FORMAT} Normalize ${CLASS_IN}0${x}${FORMAT}
    let "x+=1"
  done
fi

##################################
# 0. Run N3 bias field correction:
##################################
if [ $run_n3 -eq 1 ]; then
  ${ANTSPATH}N3BiasFieldCorrection $DIM $INPUT_FILE  $OUT_HEAD_N3 8
  ${ANTSPATH}N3BiasFieldCorrection $DIM $OUT_HEAD_N3 $OUT_HEAD_N3 4
  ${ANTSPATH}N3BiasFieldCorrection $DIM $OUT_HEAD_N3 $OUT_HEAD_N3 2
fi

########################################################################
# 1. Deform a template image to a subject image that we want to extract,
#    and apply the transform to template segmentation class images:
########################################################################

# 1a. Compute the registration transform:
if [ $run_register -eq 1 ]; then
  # Create new images/headers of the template images to be in accord with the subject:
  if [ $fix_headers -eq 1 ]; then
    exe="${ANTSPATH}ImageMath $DIM $TEMPLATE CompareHeadersAndImages $INPUT_FILE $TEMPLATE_IN"
    echo $exe; $exe
  else
    exe="cp $TEMPLATE_IN $TEMPLATE"
    echo $exe; $exe
  fi
  # Register template to the subject:
  ANTSINTRO=0
  if [ $ANTSINTRO -eq 1 ]; then
    XFM='GR'; METRIC='CC'
    PARAMS="-m $ITERS -o $OUT -s $METRIC -t $XFM"
    exe="${ANTSPATH}antsIntroduction.sh -d $DIM -r ${OUT_HEAD_N3} -i ${TEMPLATE} $PARAMS"
    echo $exe; $exe
  else
    AFFINE_ITERS="--number-of-affine-iterations 10000x10000x10000x10000x10000"
    RIGID="--rigid-affine false"
    PARAMS="-r Gauss[3,0] -t SyN[0.25] -i $ITERS --use-Histogram-Matching $AFFINE_ITERS $RIGID"
    exe="${ANTSPATH}ANTS $DIM -m CC[$OUT_HEAD_N3,$TEMPLATE,1,4] -o ${OUT}${FORMAT} $PARAMS"
    echo $exe; $exe
  fi
fi

# 1b. Apply the registration transform to the segmented class images of the template:
if [ $run_transform -eq 1 ]; then
  x=1
  let "d=$NCLASSES+1"
  while [ $x -lt $d ]; do
    # Create new images/headers of the template class images to be in accord with the subject:
    if [ $fix_headers -eq 1 ]; then
      exe="${ANTSPATH}ImageMath $DIM ${CLASS}0${x}${FORMAT} CompareHeadersAndImages $OUT_HEAD_N3 ${CLASS_IN}0${x}${FORMAT}"
      echo $exe; $exe
    else
      exe="cp ${CLASS_IN}0${x}${FORMAT} ${CLASS}0${x}${FORMAT}"
      echo $exe; $exe
    fi
    # Transform template class images to the subject:
    c_start="${ANTSPATH}WarpImageMultiTransform $DIM"
    c_end="-R $OUT_HEAD_N3 ${OUT}Warp${FORMAT} ${OUT}Affine.txt"
    exe="$c_start ${CLASS}0${x}${FORMAT} ${PRIORname}0${x}${FORMAT} $c_end"
    echo $exe; $exe
    let "x+=1"
  done
fi

###############################################################
# 2. Create a mask about the boundary of the deformed template:
###############################################################
if [ $run_mask -eq 1 ]; then
  ${ANTSPATH}MultiplyImages $DIM $OUT_HEAD_N3 0 $TEMP
  x=1
  d=4  # use three classes to create the mask
  while [ $x -lt $d ]; do
    ${ANTSPATH}ImageMath $DIM $TEMP + ${PRIORname}0${x}${FORMAT} $TEMP
    let "x+=1"
  done
  ${ANTSPATH}ThresholdImage $DIM $TEMP $MASK_TEMPLATE 0.5 $NCLASSES

  # Brains that come very close to the edge of the image volume
  # result in a mask that erodes improperly, so the outermost voxels
  # are zeroed to ensure that the mask erodes properly.
  ${ANTSPATH}ImageMath $DIM $MASK_TEMPLATE PadImage $MASK_TEMPLATE -1
  ${ANTSPATH}ImageMath $DIM $MASK_TEMPLATE PadImage $MASK_TEMPLATE 1

  ${ANTSPATH}ImageMath $DIM $MASKD MD $MASK_TEMPLATE $DILATE
  ${ANTSPATH}ImageMath $DIM $MASKE ME $MASK_TEMPLATE $ERODE
  ${ANTSPATH}ImageMath $DIM $MASKB - $MASKD $MASKE
fi

#################################################################
# 3. Run Atropos with the transformed segmented images as priors:
#################################################################
if [ $run_atropos -eq 1 ]; then
  PARAMS1="-c [2,0.001] -x $MASKB -m [$MRF,1] -u false -h 0 -e 0"
  PARAMS2="-l 1[${DIST}] -l 2[${DIST}] -l 3[${DIST}]"
  PARAMS3="-l 4[${DIST}] -l 5[${DIST}] -l 6[${DIST}]"
  PARAMS4="-l 7[${DIST}] -l 8[${DIST}]"
  if [ $NCLASSES -eq 11 ]; then
    PARAMS5=" -l 9[${DIST}] -l 10[${DIST}] -l 11[${DIST}]"
  else
    PARAMS5=""
  fi
  PARAMS="$PARAMS1 $PARAMS2 $PARAMS3 $PARAMS4 $PARAMS5"
  exe="${ANTSPATH}Atropos $DIM -a $OUT_HEAD_N3 -i PriorProbabilityImages[$NCLASSES,$PRIOR,$PRIORWT] -o [$SEG,$POST] $PARAMS"
  echo $exe; $exe
fi

##############################################################################
# 4. Extract the largest components of class 1, 2, and 3 to extract the brain:
##############################################################################
if [ $run_postprocess1 -eq 1 ]; then
  ${ANTSPATH}MultiplyImages $DIM $MASKE 1 $MASK_ATROPOS
  for x in 1 2 3 ; do
    #if [ $x -eq 1 ] ; then CLUSTERSIZE=100 ; fi
    #if [ $x -eq 2 ] ; then CLUSTERSIZE=200 ; fi
    #if [ $x -eq 3 ] ; then CLUSTERSIZE=100 ; fi
    ${ANTSPATH}ThresholdImage $DIM $SEG $TEMP $x $x
    ${ANTSPATH}LabelClustersUniquely $DIM $TEMP $TEMP $CLUSTERSIZE
    ${ANTSPATH}ThresholdImage $DIM $TEMP $TEMP 1 $NTHBIGGESTCLUST
    ${ANTSPATH}ImageMath $DIM $MASK_ATROPOS + $TEMP $MASK_ATROPOS
  done
  ${ANTSPATH}ImageMath $DIM $MASK_ATROPOS GetLargestComponent $MASK_ATROPOS
  ${ANTSPATH}ImageMath $DIM $MASK_ATROPOS MD $MASK_ATROPOS 1
  ${ANTSPATH}ImageMath $DIM $MASK_ATROPOS FillHoles $MASK_ATROPOS
  ${ANTSPATH}ImageMath $DIM $MASK_ATROPOS ME $MASK_ATROPOS 1
  ${ANTSPATH}MultiplyImages $DIM $MASK_ATROPOS $MASK_TEMPLATE $OUT_MASK

  ##############################
  # 5. Mask the resulting brain:
  ##############################
  ${ANTSPATH}MultiplyImages $DIM $MASK_TEMPLATE $INPUT_FILE $OUT_MASKED
  ${ANTSPATH}MultiplyImages $DIM $OUT_MASK $INPUT_FILE $OUT_BRAIN
  ${ANTSPATH}MultiplyImages $DIM $OUT_MASK $OUT_HEAD_N3 $OUT_BRAIN_N3

fi

##############################
# 6. Run Atropos segmentation:
##############################
if [ $run_postprocess2 -eq 1 ]; then

  # Skip postprocess1 routine
  if [ $skip_postprocess1 -eq 1 ]; then
    ${ANTSPATH}MultiplyImages $DIM $MASK_TEMPLATE $INPUT_FILE $OUT_MASKED
    ${ANTSPATH}Atropos $DIM -a $OUT_MASKED -x $MASK_TEMPLATE -m [0.25,1] -i Kmeans[3] -o $SEG2 -c [2,0.001]
  else
    ${ANTSPATH}Atropos $DIM -a $OUT_BRAIN_N3 -x $OUT_MASK -m [0.25,1] -i Kmeans[3] -o $SEG2 -c [2,0.001]
  fi

  # 3-tissue segmentation
  ${ANTSPATH}ThresholdImage $DIM $SEG2 $WM  3 3
  ${ANTSPATH}ThresholdImage $DIM $SEG2 $GM  2 2
  ${ANTSPATH}ThresholdImage $DIM $SEG2 $CSF 1 1
  # Get the largest connected component of WM and of GM and make sure each has no holes
  ${ANTSPATH}ImageMath $DIM $WM GetLargestComponent $WM
  ${ANTSPATH}ImageMath $DIM $GM GetLargestComponent $GM
  ${ANTSPATH}ImageMath $DIM $WM FillHoles $WM
  ${ANTSPATH}ImageMath $DIM $GM FillHoles $GM
  # Recreate the segmentation with these "repaired" images
  ${ANTSPATH}MultiplyImages $DIM $WM 3 $WM
  ${ANTSPATH}MultiplyImages $DIM $GM 2 $GM
  ${ANTSPATH}ImageMath $DIM $OUT_SEG + $CSF $WM
  ${ANTSPATH}ImageMath $DIM $OUT_SEG + $GM $OUT_SEG

  ###############################################
  # 7. Mask the brain with the segmentation mask:
  ###############################################
  ${ANTSPATH}ThresholdImage $DIM $OUT_SEG $OUT_MASK2 1 3
  ${ANTSPATH}ImageMath $DIM $OUT_MASK2 FillHoles $OUT_MASK2
  ${ANTSPATH}MultiplyImages $DIM $OUT_MASK2 $INPUT_FILE $OUT_BRAIN2
  ${ANTSPATH}MultiplyImages $DIM $OUT_MASK2 $OUT_HEAD_N3 $OUT_BRAIN2_N3

fi

