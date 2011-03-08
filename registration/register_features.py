"""
Feature-based registration

(c) Arno Klein (2011)  .  arno@mindboggle.info  .  http://www.mindboggle.info
"""
                                                                               
import shlex
from subprocess import Popen

#                                                                        Inputs
ANTSPATH = "/Users/arno/Software/ANTS-build/"
source = "sourceFILE"
target = "targetFILE"
output_stem = "outFILE"
landmark_string = "_feature"
source_landmarks = source + landmark_string
target_landmarks = target + landmark_string
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
partial_matching_iters = 100000
landmark_weight = 0.75
landmark_percent = 100
landmark_sigma = 5

#                                                                     Arguments

command1 = ANTSPATH + "ANTS " + str(dim)
command2 = ANTSPATH + "WarpImageMultiTransform " + str(dim)

intensity_based = "-m CC[" + ", ".join([target, source, str(measure_weight), str(measure_radius)]) + "]"

landmark_based = "-m PSE[" + ", ".join([target, source, target_landmarks, source_landmarks, str(landmark_weight), str(landmark_percent), str(landmark_sigma), str(boundary_only), str(kNeighborhood), str(partial_matching_iters)]) + "]"

transform = "-t SyN[" + str(transform_gradient_step_size) +"] -i " + str(iterations) + " --use-Histogram-Matching"

initialize = "--number-of-affine-iterations " + str(affine_iterations)

regularize = "-r Gauss[" + str(regularize_similarity_gradient_sigma) + ", " + str(regularize_deformation_field_sigma) + "]"

output = "-o " + output_stem

"""
$count = 0;
while [ $count -lt $number_of_features ]; do
    Value of count is: $COUNT
    let count = count+1
done 
"""

#                                                                      Commands
args = [command1, intensity_based, landmark_based, output, regularize, transform, initialize]
print(args)
#p = subprocess.Popen(args)

args = [command2, source, output, '-R '+target, output+'Warp.nii.gz', output+'Affine.txt']
print(args)
#p = subprocess.Popen(args)
