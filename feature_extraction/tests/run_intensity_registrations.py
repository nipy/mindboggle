#!/usr/bin/python
"""
Run registration commands

(c) 2011, @rno klein
"""

register_to_template = 0
register_pairs_via_template = 1

import os
from os.path import exists
from subprocess import call
from numpy import random, float, isnan

from measure_overlap import measure_overlap

#
# Inputs
#
verbose = 1
register = 0
register_extras = 0
fill_target_mask = 1
measure_overlap = 1
dim = 3

#------------
# Directories
#------------
source_files       = ['m1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m11','m12']
target_files       = source_files
ANTSPATH           = os.environ.get("ANTSPATH")
FSLPATH            = '/usr/local/fsl/bin/'
out_path           = '/hd2/Archive/registration_evaluation_2011_output/'
xfm_dir            = os.path.join( out_path, 'Transforms/')
xfm_brain_dir      = os.path.join( out_path, 'Transformed_Brains/')
xfm_atlas_dir      = os.path.join( out_path, 'Transformed_Atlases/')
atlas_dir          = '/hd2/Brains/CUMC12/Atlases/'
brain_dir          = '/hd2/Brains/CUMC12/Brains/'
ext                = '.nii.gz'
template           = '/hd2/Brains/CUMC12/CUMC12template.nii.gz'
xfm_via_template   = 1

use_landmarks      = 1
#landmarks          = 'pits_kiho_im_binary'
#landmarks          = 'fundi_brain_visa_binary'
#landmarks          = 'fundi_gang_li_binary'

#
# Registration parameters
#
gradient_step_size = 0.5
iterations = "30x100x10"
options = " --use-Histogram-Matching"
initialize =  " --number-of-affine-iterations 10000x10000x10000x10000x10000"
intensity_weight = 1.0
landmark_weight1 = 1.0
landmark_weight2 = 1.0
#
# Regularization parameters
#
regularizer = "Gauss"
regularizer_setting = 3
deformation_field_sigma = 0
#
# Intensity parameters
#
intensity_measure = "CC"
intensity_setting = 3
#
# Landmark parameters
#
landmark_measure1 = "PSE"
landmark_measure2 = "MSQ"
percent = 1.0 # real number: 1.0 = 100%
boundary = 0  # 0: not only boundaries
sigma = 10
neighbors = 1
matching_iter = 100000 # partial matching iterations

#
# Run!
#
if use_landmarks:


# Intensity-based registration:
else:
	for file in source_files:
		cmd = ANTSPATH + 'ANTS ' + dim + ' -m CC[' + template + ',' + \
		      brain_dir+file+ext+',1,3] -o ' + xfm_dir+file+ext + \
		      ' -r Gauss[3,0] -t SyN['+gradient_step_size+']  -i ' + iterations + ' --use-Histogram-Matching ' + \
		      ' --number-of-affine-iterations 10000x10000x10000x10000x10000 '
		print cmd; #os.system(cmd)


if os.path.exists(source) and os.path.exists(target) and os.path.exists(source_labels) and os.path.exists(target_labels):
  results_file = results_dir+output_stem+"_"+regularizer+"_"+intensity_measure+"_"+landmark_measure1+"_"+landmark_measure2+".txt"
  if save_results and os.path.exists(results_file) and overwrite == 0:
      raise NameError('File already exists: '+results_file+' -- set "overwrite to 1."')
  else:
    if save_results:
      f_eval = open(results_file, 'w')
    print('Number of tests: ' + str(len(gradient_step_sizes)*len(regularizer_settings)*len(intensity_weights)*len(landmark_weights1)*len(landmark_weights2)*len(sigmas)*len(neighbors)))

                    output_stem = file + '_to_' + file2
                    output = outpath + output_stem

                    #
                    # Arguments for ANTS
                    #

                    # Output:
                    out = "-o " + output+ext  #_file

                    # Transform:
                    warp = ANTSPATH + "ANTS " + str(dim)
                    apply_warp = ANTSPATH + "WarpImageMultiTransform " + str(dim)

                    transform = "-t SyN[" + str(gradient_step_size) +"] -i " + \
                                str(iterations) + options + initialize

                    # Regularize:
                    if regularizer == "Gauss":
                        regularize = "-r Gauss[" + str(regularizer_setting) + ", " + \
                                     str(deformation_field_sigma) + "]"
                    else:
                        regularize = "-r DMFFD["+ regularizer_setting + ", 0, " + str(bspline_order) + "]"

                    # Intensity similarity:
                    intensity = [target, source, intensity_weight, intensity_setting]
                    intensity = " -m "+intensity_measure+"[" + ", ".join([str(s) for s in intensity]) + "]"

                    # Landmark similarity:
                    if multiple_landmarks:
                        landmarks1 = ""
                        landmarks2 = ""
                        for ilabel, label in enumerate(labels):
                            source_landmarks = source_landmarks_path + str(label) + ext
                            target_landmarks = target_landmarks_path + str(label) + ext
                            if landmark_measure1 == "PSE":
                                lm_args1 = [target, source, target_landmarks, source_landmarks,
                                        landmark_weight1, percent, sigma,
                                        boundary, neighbor, matching_iter]
                                landmarks1 += ", ".join([" -m PSE[" + ", ".join([str(s) for s in lm_args1]) + "]"])
                            if landmark_measure2 == "MSQ":
                                lm_args2 = [target_landmarks, source_landmarks, landmark_weight2, 0]
                                landmarks2 += " ".join([" -m MSQ[" + ", ".join([str(s) for s in lm_args2]) + "]"])
                    else:
                        lm_args1 = [target, source, target_landmarks, source_landmarks,
                                landmark_weight1, percent, sigma,
                                boundary, neighbor, matching_iter]
                        landmarks1 = ", ".join([" -m PSE[" + ", ".join([str(s) for s in lm_args1]) + "]"])
                        lm_args2 = [target_landmarks, source_landmarks, landmark_weight2, 0]
                        landmarks2 = " ".join([" -m MSQ[" + ", ".join([str(s) for s in lm_args2]) + "]"])

                    #
                    # Run commands
                    #
                    args1 = " ".join([warp, out, transform, regularize, intensity, landmarks1, landmarks2])
                    if register:
                        if verbose: print(args1); print('')
                        p = call(args1, shell="True")

                        args2 = " ".join([apply_warp, source_labels, output+'_labels'+ext, '-R ' + target_labels, output+'Warp'+ext, output+'Affine.txt', '--use-NN'])
                        if verbose: print(args2); print('')
                        p = call(args2, shell="True")

                    if register_extras:
                        args3 = " ".join([apply_warp, source, output+ext, '-R ' + target, output+'Warp'+ext, output+'Affine.txt'])
                        if verbose: print(args3); print('')
                        p = call(args3, shell="True")

                        args4 = " ".join([apply_warp, source_landmarks, output+'_landmarks'+ext, '-R ' + target, output+'Warp'+ext, output+'Affine.txt', '--use-NN'])
                        if verbose: print(args4); print('')
                        p = call(args4, shell="True")

                    if fill_target_mask:
                        target_mask_file = target_labels
                        args5 = " ".join(['ImageMath',str(dim),output+'_labels_filled'+ext,'PropagateLabelsThroughMask',target_mask_file,output+'_labels'+ext])
                        if verbose: print(args5); print('')
                        p = call(args5, shell="True")

                    #
                    # Measure overlap
                    #
                    if measure_overlap:
                        print(output)
                        if save_results:
                            f_eval.close()
                            f_eval = open(results_file, 'a')
                            f_eval.write(output + '\n\n' + args1 + '\n\n')
                            
                        average_dice = 0
                        average_jacc = 0
                        for label in labels:
                            args6 = " ".join(['c3d', output+'_labels'+ext, target_labels, '-overlap', str(label), '>'+output_dir+'temp_overlap.txt'])
                            p = call(args6, shell="True")
                            f = open(output_dir+'temp_overlap.txt','r')
                            temp = f.read()
                            dice = float(temp.split()[-2].split(',')[0])
                            jacc = float(temp.split()[-1].split(',')[0])
                            print_out = ' '.join(['Label:', str(label), 'Dice:', str(dice), 'Jaccard:', str(jacc)])
                            print(print_out)
                            if save_results:
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
                        print(print_out1)
                        print(print_out2)
                        if save_results:
                            f_eval.close()
                            f_eval = open(results_file, 'a')
                            f_eval.write('\n' + output_stem + '\n' + print_out1 + '\n' + print_out2 + '\n\n')
    if save_results:
        f_eval.close()
else:
    if not os.path.exists(source):
        raise NameError('Check input file ' + source)
    elif not os.path.exists(target):
        raise NameError('Check input file ' + target)
    elif not os.path.exists(source_labels):
        raise NameError('Check input file ' + source_labels)
    elif not os.path.exists(target_labels):
        raise NameError('Check input file ' + target_labels)



#------------------------------
# Register subjects to template
#------------------------------
if register_to_template:
	for file in source_files:
	
        # Landmark + intensity-based registration:
        if use_landmarks:
            for file2 in target_files:
            
				warp = ANTSPATH + "ANTS " + str(dim)
				apply_warp = ANTSPATH + "WarpImageMultiTransform " + str(dim)
				
				transform = "-t SyN[" + str(gradient_step_size) +"] -i " + \
				            str(iterations) + options
				
				regularize = "-r Gauss[" + str(similarity_gradient_sigma) + ", " + \
				             str(deformation_field_sigma) + "]"
				
				out = "-o " + output_file
				
				intensity = [target_file, source_file, intensity_weight, intensity_radius]
				intensity = "-m CC[" + ", ".join([str(s) for s in intensity]) + "]"
                
				args = [warp, intensity, landmarks, out, regularize, transform, initialize]
				print(" ".join(args))
				#p = Popen(args)

			    cmd = ANTSPATH + 'ANTS ' + dim + ' -m CC[' + ants_template + ',' + \
			          brain_dir+file+ext+',1,3] -o ' + xfm_dir+file+ext + \
			          ' -r Gauss[3,0] -t SyN['+gradient_step_size+']  -i ' + iterations + ' --use-Histogram-Matching ' + \
			          ' --number-of-affine-iterations 10000x10000x10000x10000x10000 '
			    print cmd; #os.system(cmd)

        # Intensity-based registration:
        else:
        
		    cmd = ANTSPATH + 'ANTS ' + dim + ' -m CC[' + ants_template + ',' + \
		          brain_dir+file+ext+',1,3] -o ' + xfm_dir+file+ext + \
		          ' -r Gauss[3,0] -t SyN['+gradient_step_size+']  -i ' + iterations + ' --use-Histogram-Matching ' + \
		          ' --number-of-affine-iterations 10000x10000x10000x10000x10000 '
		    print cmd; #os.system(cmd)

#----------------------
# Pairwise registration
#----------------------
if register_pairs_via_template:
	for file in source_files:
	    for file2 in target_files:
	        if file2 != file:
		        
		        # Transform brains
		        cmd = ANTSPATH + 'WarpImageMultiTransform ' + dim + ' ' + brain_dir+file+ext + ' ' + \
		              xfm_brain_dir+file+'_to_'+file2+ext + \
		              ' -R ' + brain_dir+file2+ext + \
		              ' -i ' + xfm_dir+file2+'Affine.txt ' + \
		              xfm_dir+file2+'InverseWarp.nii.gz ' + \
		              xfm_dir+file+'Warp.nii.gz ' + \
		              xfm_dir+file+'Affine.txt --use-NN '
		        print cmd; os.system(cmd)
		        #cmd = FSLPATH + 'fslswapdim ' + xfm_brain_dir+file+'_to_'+file2+ext + \
                #      ' -x y z ' + xfm_brain_dir+file+'_to_'+file2+ext
		        #print cmd; os.system(cmd)
		        
		        # Transform atlases
		        cmd = ANTSPATH + 'WarpImageMultiTransform ' + dim + ' ' + atlas_dir+file+ext + ' ' + \
		              xfm_atlas_dir+file+'_to_'+file2+ext + \
		              ' -R ' + atlas_dir+file2+ext + \
		              ' -i ' + xfm_dir+file2+'Affine.txt ' + \
		              xfm_dir+file2+'InverseWarp.nii.gz ' + \
		              xfm_dir+file+'Warp.nii.gz ' + \
		              xfm_dir+file+'Affine.txt --use-NN '
		        print cmd; os.system(cmd)
		        #cmd = ANTSPATH + 'PermuteFlipImageOrientationAxes 3 ' + xfm_atlas_dir+file+'_to_'+file2+ext + \
                #      ' ' + xfm_atlas_dir+file+'_to_'+file2+ext + ' 0 1 2 1 0 0 '
		        #print cmd; os.system(cmd)
		        #cmd = FSLPATH + 'fslswapdim ' + xfm_atlas_dir+file+'_to_'+file2+ext + \
                #      ' -x y z ' + xfm_atlas_dir+file+'_to_'+file2+ext
		        #print cmd; os.system(cmd)
		        oarisnet
