#!/bin/sh
# hi arno -- a minor but interesting addition to the lohmann script.  
# realized later that it's easy to locate sulci on the surface by the procedure outlined below.
# (c) brian avants, 2009/02

NUMPARAMS=$#

if [ $NUMPARAMS -lt 2  ]
then 
echo " USAGE ::  "
echo "  sh  $0  binarywhitemattersegmentation.nii outputnamingprefix "
exit
fi


WM=$1 # should be binary image of white matter 

OUT=$2 # the output naming prefix 

if  [ ${#WM} -lt 1 -o  ! -f $WM ]
then
echo " Problem with specified Image => User Specified Value = $WM " 
exit
fi

# we would prefer to use topology models in all below 

# the closing operation - arachnoid surface
 ImageMath 3 ${OUT}closed.nii MD $WM  15
 ImageMath 3 ${OUT}closed.nii ME ${OUT}closed.nii  15

# unfolded or lissencephalic surface 
# we would prefer to do this in a topology preserving way
 ImageMath 3 ${OUT}opened.nii ME $WM 5
 ImageMath 3 ${OUT}opened.nii MD ${OUT}opened.nii 5

# distance transform from arachnoid surface to inside of brain
# perhaps should be a geodesic distance transform 
 ImageMath 3 ${OUT}closed.nii Neg ${OUT}closed.nii
 ImageMath 3 ${OUT}distance.nii D ${OUT}closed.nii

# multiply the distance transform by the "unfolded" surface
# this gives the "feature" surface/image that could be input 
# to an image registration similarity metric 
ImageMath 3 ${OUT}lohmann.nii m ${OUT}distance.nii ${OUT}opened.nii 

# surface stuff for visualization 
ImageMath 3 ${OUT}surface.nii ME ${OUT}opened.nii 2
ImageMath 3 ${OUT}surface.nii - ${OUT}opened.nii ${OUT}surface.nii 
ImageMath 3 ${OUT}maxvals.nii m ${OUT}lohmann.nii ${OUT}surface.nii 

# find sulci using curvature of the WM/lohmann surface 
SurfaceCurvature ${OUT}lohmann.nii ${OUT}lohmann_curvature.nii 1.5
    ImageMath 3 ${OUT}lohmann_curvature.nii m ${OUT}lohmann_curvature.nii -1
   ThresholdImage 3 ${OUT}lohmann_curvature.nii ${OUT}lohmann_threshold.nii 0.0001 1.e9
# get sulci near the surface
        ImageMath 3 ${OUT}sulci.nii m ${OUT}surf.nii ${OUT}lohmann_threshold.nii 
# get sulci longer than $N voxels     
     N=100
LabelClustersUniquely 3 ${OUT}sulci_clustered.nii ${OUT}sulci.nii $N
echo " Found # Sulci "
        MeasureMinMaxMean 3 ${OUT}sulci_clustered.nii 

