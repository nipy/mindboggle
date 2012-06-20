"""
Feature-based registration

Command:

Example:

Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2011  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""
                                                                               
import os
from subprocess import Popen

#
# Inputs
#
ANTSPATH = os.environ.get("ANTSPATH")
source = "test/S01"
target = "test/S02"
output = "S01toS02_pits075fundi05int1_100x100x100"
ext = ".nii.gz"
source_file = source + ext
target_file = target + ext
output_file = output + ext

#
# Intensity registration parameters
#
dim = 3
gradient_step_size = 0.25
iterations = "100x100x100" #"100x100x100"
similarity_gradient_sigma = 3
deformation_field_sigma = 0
intensity_radius = 4
intensity_weight = 1 #0.1, 1
options = " --use-Histogram-Matching"
initialize = " --number-of-affine-iterations 10000x10000x10000x10000x10000"

#
# Landmark registration parameters
#
labels = ["pits", "fundi"]
weights = [0.75, 0.5] #[0.1,0.1] #[1.0,0.01] #[0.75, 0.5] # adds up to 1.25, which is > intensity_weight
percents = [0.99, 0.99]  # real number
sigmas = [25, 25]  # need experiments with parzen models of the data
                   # (big numbers are nearly uniform distributions)
boundaries = [0, 0]
neighbors = [10, 10]
matching_iters = [100000, 100000]  # partial matching iterations

#
# Arguments
#
warp = ANTSPATH + "ANTS " + str(dim)
apply_warp = ANTSPATH + "WarpImageMultiTransform " + str(dim)

transform = "-t SyN[" + str(gradient_step_size) +"] -i " + \
            str(iterations) + options

regularize = "-r Gauss[" + str(similarity_gradient_sigma) + ", " + \
             str(deformation_field_sigma) + "]"

out = "-o " + output_file

intensity = [target_file, source_file, intensity_weight, intensity_radius]
intensity = "-m CC[" + ", ".join([str(s) for s in intensity]) + "]"

#
# Arguments for multiple features as landmarks
#
landmarks = ""
for i, label in enumerate(labels):
    args = [target_file, source_file, target+label+ext, source+label+ext,
            weights[i], percents[i], sigmas[i],
            boundaries[i], neighbors[i], matching_iters[i]]
    landmarks = " ".join([landmarks, "-m PSE[" + ", ".join([str(s) for s in args]) + "]"])

#
# Run commands
#
if os.path.exists(source_file) and os.path.exists(target_file):
	args = [warp, intensity, landmarks, out, regularize, transform, initialize]
	print(" ".join(args))
	#p = Popen(args)
	
	args = [apply_warp, source_file, output_file, '-R ' + target_file, output+'Warp'+ext, output+'Affine.txt']
	print(" ".join(args))
	#p = Popen(args)
