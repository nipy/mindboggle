#!/bin/bash
#
# Inputs
#
source=$1
target=$2
output=$3
landmark_string="_feature"
number_of_features=2

#
# Registration parameters
#
dim=3
transform_gradient_step_size=0.5
iterations="30x100x10"
affine_iterations="10000x10000x10000x10000x10000"
regularize_similarity_gradient_sigma=3
regularize_deformation_field_sigma=0
measure_radius=2
measure_weight=0.75

#
# Landmark-driven registration parameters
#
boundary_only=0
kNeighborhood=5
$partial_matching_iters=100000
landmark_weight1=0.75
landmark_percent1=100
landmark_sigma1=5
landmark_weight2=0.75
landmark_percent2=100
landmark_sigma2=5

#
# Arguments
#
intensity_difference="CC[$target, $source, $measure_weight, $measure_radius]"

transform="SyN[$transform_gradient_step_size] -i $iterations --use-Histogram-Matching"

initialize="--number-of-affine-iterations $affine_iterations"

regularize="Gauss[$regularize_similarity_gradient_sigma, $regularize_deformation_field_sigma]"

landmark_guidance1="PSE[$target,$source,$target_landmarks1,$source_landmarks1,$landmark_weight1,$landmark_percent1,$landmark_sigma1,$boundary_only,$kNeighborhood,$partial_matching_iters]"

landmark_guidance2="PSE[$target,$source,$target_landmarks2,$source_landmarks2,$landmark_weight2,$landmark_percent2,$landmark_sigma2,$boundary_only,$kNeighborhood,$partial_matching_iters]"

$count=0;
while [ $count -lt $number_of_features ]; do
    Value of count is: $COUNT
    let count=count+1
done 

#
# Commands
#
${ANTSPATH}/ANTS $dim -m $intensity_difference -m $landmark_guidance -o $output -r $regularize -t $transform $initialize

${ANTSPATH}/WarpImageMultiTransform $dim $source $output -R $target ${output}Warp.nii.gz ${output}Affine.txt

${ANTSPATH}/WarpImageMultiTransform $dim $source_landmarks ${output}_landmarks -R $target_landmarks ${output}Warp.nii.gz ${output}Affine.txt --use-NN