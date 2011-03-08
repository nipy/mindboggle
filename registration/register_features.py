"""
Feature-based registration

(c) Arno Klein (2011)  .  arno@mindboggle.info  .  http://www.mindboggle.info
"""
                                                                               
#                                                                        Inputs
source = 1
target = 2
output_stem = 3
landmark_string = "_feature"
number_of_features = 1

#                                                       Registration parameters
dim = 3
transform_gradient_step_size = 0.5
iterations = "30x100x10"
affine_iterations = "10000x10000x10000x10000x10000"
regularize_similarity_gradient_sigma = 3
regularize_deformation_field_sigma = 0
measure_radius = 2
measure_weight = 0.75

#                                              Landmark registration parameters
boundary_only = 0
kNeighborhood = 5
$partial_matching_iters = 100000
landmark_weight1 = 0.75
landmark_percent1 = 100
landmark_sigma1 = 5
landmark_weight2 = 0.75
landmark_percent2 = 100
landmark_sigma2 = 5

#                                                                     Arguments

command1 = ANTSPATH + "ANTS " + dim
command2 = ANTSPATH + "WarpImageMultiTransform " + dim

intensity_based = "-m CC[$target, $source, $measure_weight, $measure_radius]"

landmark_based = "-m PSE[$target,$source,$target_landmarks1,$source_landmarks1,$landmark_weight1,$landmark_percent1,$landmark_sigma1,$boundary_only,$kNeighborhood,$partial_matching_iters]"

transform = "-t SyN[$transform_gradient_step_size] -i $iterations --use-Histogram-Matching"

initialize = "--number-of-affine-iterations $affine_iterations"

regularize = "-r Gauss[$regularize_similarity_gradient_sigma, $regularize_deformation_field_sigma]"

output = "-o " + output_stem

"""
$count = 0;
while [ $count -lt $number_of_features ]; do
    Value of count is: $COUNT
    let count = count+1
done 
"""

#                                                                      Commands

cmd = ", ".join([command1, intensity_based, landmark_based, output, regularize, transform, initialize])
print(cmd)

cmd = ", ".join([command2, source, output, "-R "+target, output+"Warp.nii.gz", output+"Affine.txt"])
print(cmd)