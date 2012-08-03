#!/usr/bin/python
"""
Register landmarks in brain images.

This program uses ANTs SyN registration method.

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

"""
Example code:

from os.path import exists
from subprocess import call
from numpy import float, isnan
dim = 3

#=========================
# Transform via a template
#=========================
if not os.path.exists(xfm_brain_dir+output_stem+ext):
    args = " ".join([ANTSPATH+'WarpImageMultiTransform', str(dim), \
            source, xfm_brain_dir+output_stem+ext, '-R', target, \
            '-i', inv_xfm_stem+'Affine.txt', \
            inv_xfm_stem+'InverseWarp.nii.gz', \
            xfm_stem+'Warp.nii.gz', \
            xfm_stem+'Affine.txt'])
    if verbose: print(args); print(''); p = call(args, shell="True")

#============================================================
# Fill target atlas mask with transformed source atlas labels
#============================================================
args = " ".join(['ImageMath', str(dim), xfm_atlas_dir+output_stem+'_filled'+ext, \
                 'PropagateLabelsThroughMask', brainmask_dir+file2+ext, \
                                               xfm_atlas_dir+output_stem+ext])
if verbose: print(args); print(''); p = call(args, shell="True")

#====================================================================
# Measure overlap of target atlas and transformed source atlas labels
#====================================================================
results_file = results_dir+output_stem+'.txt'
f_eval = open(results_file, 'w');
average_dice = 0
average_jacc = 0

for label in labels:
    args = " ".join(['c3d', xfm_atlas_dir+output_stem+'_filled'+ext, \
                            atlas_dir+file2+ext, '-overlap', str(label), \
                            '>'+results_dir+'temp_overlap.txt'])
    p = call(args, shell="True")
    f = open(results_dir+'temp_overlap.txt','r')
    temp = f.read()
    dice = 0
    jacc = 0
    if temp != '':
        dice = float(temp.split()[-2].split(',')[0])
        jacc = float(temp.split()[-1].split(',')[0])
    print_out = " ".join(['Label:', str(label), 'Dice:', str(dice), \
                                             'Jaccard:', str(jacc)])
    print(print_out)
    f_eval.close()
    f_eval = open(results_file, 'a')
    f_eval.write(print_out + '\n')
    if isnan(dice):
        dice = 0
    if isnan(jacc):
        jacc = 0
    average_dice += dice
    average_jacc += jacc
average_dice = average_dice/len(labels)
average_jacc = average_jacc/len(labels)
print_out1 = 'Average Dice: ' + str(average_dice)
print_out2 = 'Average Jacc: ' + str(average_jacc)
print(print_out1);
print(print_out2)
f_eval.close()
f_eval = open(results_file, 'a')
f_eval.write('\n' + print_out1 + '\n' + print_out2 + '\n\n')
f_eval.close()
f_avg.close()
f_avg = open(avg_results_file, 'a');
f_avg.write(output_stem + ' ' + str(average_dice) + ' ' + str(average_jacc) + '\n')
"""
