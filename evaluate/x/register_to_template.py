#!/usr/bin/python
"""
Register brains, landmarks, and labels to a template.

(c) 2011, @rno klein
"""

import os
from os.path import exists
from subprocess import call
from numpy import float, isnan

# Run intensity-based registration
# 1. Register brains to template
# 2. Transform brains to each other via template
# 3. Transform landmarks to template
register_to_template = 1
transform_pairs_via_template = 1
transform_landmarks_to_template = 0

# Run landmark-driven registration to template:
register_landmarks_to_template = 0
transform_landmarks_via_template = 0

# Atlas-based evaluation for the above settings:
# 1. prepare target atlas mask
# 2. transform source atlas
# 3. fill #1 with #2
# 4. measure overlap of #3 with target atlas labels
prepare_target_mask = 0
evaluate_with_atlases = 1

verbose = 1
dim = 3

#
# Files
#
source_files       = ['m1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m11','m12']
target_files       = ['m1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m11','m12']
#source_files       = ['m1','m2','m3','m4']#,'m5','m6']
#target_files       = ['m1','m2','m3','m4']#,'m5','m6']
ANTSPATH           = os.environ.get("ANTSPATH")
FSLPATH            = '/usr/local/fsl/bin/'
out_path           = '/hd2/Archive/registration_evaluation_2011_output/'
xfm_dir            = os.path.join( out_path, 'Transforms/')
xfm_brain_dir      = os.path.join( out_path, 'Transformed_Brains/')
xfm_landmarks_dir   = os.path.join( out_path, 'Transformed_Landmarks/')
xfm_atlas_dir      = os.path.join( out_path, 'Transformed_Atlases/')
atlas_dir          = '/hd2/Brains/CUMC12/Atlases/'
brain_dir          = '/hd2/Brains/CUMC12/Brains/'
brainmask_dir      = '/hd2/Brains/CUMC12/BrainMasks/'
ext                = '.nii.gz'
template           = '/hd2/Brains/CUMC12/CUMC12template.nii.gz'

landmarks_dir      = '/hd2/Brains/CUMC12/Landmarks/pits_kiho_im_binary/'
landmark_type      = 'pits_kiho_im'
landmarks_dir      = '/hd2/Brains/CUMC12/Landmarks/pits_yrjo_hame_binary/'
landmark_type      = 'pits_yrjo_hame'
landmarks_dir      = '/hd2/Brains/CUMC12/Landmarks/pits_forrest_bao_binary/'
landmark_type      = 'pits_forrest_bao'
landmarks_dir      = '/hd2/Brains/CUMC12/Landmarks/ribbons_brain_visa_binary/'
landmark_type      = 'ribbons_brain_visa'
landmarks_dir      = '/hd2/Brains/CUMC12/Landmarks/fundi_gang_li_binary/'
landmark_type      = 'fundi_gang_li'
landmarks_dir      = '/hd2/Brains/CUMC12/Landmarks/fundi_brain_visa_binary/'
landmark_type      = 'fundi_brain_visa'
landmarks_dir      = '/hd2/Brains/CUMC12/Landmarks/fundi_forrest_bao_binary/'
landmark_type      = 'fundi_forrest_bao'

results_dir        = os.path.join( out_path, 'Results/')
label_file         = 'CUMC12_labels_regions.txt'

#
# Registration parameters
#
gradient_step_size = 0.5
iterations = "30x100x10"
options = " --use-Histogram-Matching"
initialize =  " --number-of-affine-iterations 10000x10000x10000x10000x10000"
warp = ANTSPATH + "ANTS " + str(dim) + " -t SyN[" + str(gradient_step_size) +"] -i " + \
                            str(iterations) + options + initialize
apply_warp = ANTSPATH + "WarpImageMultiTransform " + str(dim)
#
# Regularization parameters
#
regularizer = "Gauss"
regularizer_setting = 3
deformation_field_sigma = 0
regularize = "-r Gauss[" + str(regularizer_setting) + ", " + \
                           str(deformation_field_sigma) + "]"
#
# Intensity parameters
#
intensity_measure = "CC"
intensity_weight = 1.0
intensity_setting = 3
#
# Landmark parameters
#
landmark_measure1 = "PSE"
landmark_measure2 = "MSQ"
landmark_weight1 = 0.1
landmark_weight2 = 0.1
percent = 1.0 # real number: 1.0 = 100%
boundary = 0  # 0: not only boundaries
sigma = 10
neighbor = 100
matching_iter = 100000 # partial matching iterations

if evaluate_with_atlases:
    f = open(label_file,'r')
    label_table = f.readlines()
    f.close()
    labels = []
    for row in label_table:
        labels.append(int(row.split()[0]))

#------------------------------------------
# Register brains and landmarks to template
#------------------------------------------
if register_to_template + transform_landmarks_to_template + \
   prepare_target_mask > 0:
    for file in source_files:
        source = brain_dir+file+ext
        output = xfm_dir+file+'_to_template'
        out = '-o ' + output+ext
        if os.path.exists(source) and os.path.exists(template) and os.path.exists(xfm_dir):

            # Intensity-based registration to template:
            if register_to_template:
                intensity = [template, source, intensity_weight, intensity_setting]
                intensity = "-m "+intensity_measure+"[" + ", ".join([str(s) for s in intensity]) + "]"
                args = " ".join([warp, regularize, intensity, out])
                if verbose: print(args); print(''); p = call(args, shell="True")

            # Prepare binary (target atlas) masks for filling with labels:
            if prepare_target_mask:
                args = " ".join(['c3d', atlas_dir+file+ext, '-binarize -o', brainmask_dir+file+ext])
                if verbose: print(args); print(''); p = call(args, shell="True")

            # Transform landmarks to template space:
            if transform_landmarks_to_template:
                source_landmarks = landmarks_dir+file+ext
                output_landmarks = xfm_landmarks_dir+file+'_to_template_'+landmark_type+ext
                try:
                    os.path.exists(source_landmarks) and os.path.exists(xfm_landmarks_dir)
                except:
                    raise NameError('Check ' + source_landmarks + ' and ' + xfm_landmarks_dir)
                args = " ".join([apply_warp, source_landmarks, output_landmarks, \
                                 '-R', template, output+'Warp'+ext, output+'Affine.txt', '--use-NN'])
                if verbose: print(args); print(''); p = call(args, shell="True")
        else:
            if not os.path.exists(source):
                raise NameError('Check input file ' + source)
            elif not os.path.exists(template):
                raise NameError('Check input file ' + template)
            elif not os.path.exists(xfm_dir):
                raise NameError('Check input file ' + xfm_dir)

#--------------------------------------------------------------
# Register landmarks to transformed landmarks in template space
#--------------------------------------------------------------
if register_landmarks_to_template:
    for file in source_files:
        source = brain_dir+file+ext
        source_landmarks = landmarks_dir+file+ext
        for file2 in target_files:
            if file2 != file:
                template_landmarks = xfm_landmarks_dir+file2+'_to_template_'+landmark_type+ext
                output_xfm = xfm_dir+file+'_to_'+file2+'_in_template_space_'+landmark_type+ext
                if os.path.exists(source) and os.path.exists(template) and \
                   os.path.exists(source_landmarks) and os.path.exists(template_landmarks):

                    # Intensity similarity:
                    intensity = [template, source, intensity_weight, intensity_setting]
                    intensity = " -m "+intensity_measure+"[" + ", ".join([str(s) for s in intensity]) + "]"

                    # Landmark similarity:
                    lm_args1 = [template, source, template_landmarks, source_landmarks,
                                landmark_weight1, percent, sigma, boundary, neighbor, matching_iter]
                    landmarks1 = ", ".join([" -m PSE[" + ", ".join([str(s) for s in lm_args1]) + "]"])
                    lm_args2 = [template_landmarks, source_landmarks, landmark_weight2, 0]
                    landmarks2 = " ".join([" -m MSQ[" + ", ".join([str(s) for s in lm_args2]) + "]"])

                    #
                    # Run command
                    #
                    args = " ".join([warp, '-o', output_xfm, regularize, intensity, landmarks1, landmarks2])
                    if verbose: print(args); print(''); p = call(args, shell="True")
                else:
                    if not os.path.exists(source):
                        raise NameError('Check input file ' + source)
                    elif not os.path.exists(template):
                        raise NameError('Check input file ' + template)
                    elif not os.path.exists(source_landmarks):
                        raise NameError('Check input file ' + source_landmarks)
                    elif not os.path.exists(template_landmarks):
                        raise NameError('Check input file ' + template_landmarks)

#----------------------------------------------
# Apply intensity-based registration transforms
# to register brains to each other via template
#----------------------------------------------
if transform_pairs_via_template:
    if evaluate_with_atlases:
        avg_results_file = results_dir+'dice_jacc_overlaps.txt'
        f_avg = open(avg_results_file, 'w');
    for file in source_files:
        source = brain_dir+file+ext
        for file2 in target_files:
            if file2 != file:
                target = brain_dir+file2+ext
                if os.path.exists(brain_dir+file+ext) and \
                    os.path.exists(brain_dir+file2+ext) and \
                    os.path.exists(xfm_dir+file+'_to_templateWarp.nii.gz'):
                    output_stem = file + '_to_' + file2

                    # Transform brains
                    args = " ".join([ANTSPATH+'WarpImageMultiTransform', str(dim), \
                            source, xfm_brain_dir+output_stem+ext, '-R',target, \
                            '-i', xfm_dir+file2+'_to_templateAffine.txt', \
                            xfm_dir+file2+'_to_templateInverseWarp.nii.gz', \
                            xfm_dir+file+'_to_templateWarp.nii.gz', \
                            xfm_dir+file+'_to_templateAffine.txt'])
                    #if verbose: print(args); print(''); p = call(args, shell="True")

                    if evaluate_with_atlases:

                        # Transform atlases
                        source_labels = atlas_dir+file+ext
                        target_labels = atlas_dir+file2+ext
                        args = " ".join([ANTSPATH+'WarpImageMultiTransform', str(dim), \
                                source_labels, xfm_atlas_dir+output_stem+ext, '-R', target_labels, \
                                '-i', xfm_dir+file2+'_to_templateAffine.txt', \
                                xfm_dir+file2+'_to_templateInverseWarp.nii.gz', \
                                xfm_dir+file+'_to_templateWarp.nii.gz', \
                                xfm_dir+file+'_to_templateAffine.txt','--use-NN'])
                        #if verbose: print(args); print(''); p = call(args, shell="True")

                        # Fill target atlas mask with transformed source atlas labels
                        args = " ".join(['ImageMath', str(dim), xfm_atlas_dir+output_stem+'_filled'+ext, \
                                         'PropagateLabelsThroughMask', brainmask_dir+file2+ext, \
                                                                       xfm_atlas_dir+output_stem+ext])
                        #if verbose: print(args); print(''); p = call(args, shell="True")

                        # Measure overlap of target atlas and transformed source atlas labels
                        results_file = results_dir+output_stem+'.txt'
                        f_eval = open(results_file, 'w');
                        average_dice = 0
                        average_jacc = 0
                        print(results_file)

                        for label in labels:
                            args = " ".join(['c3d', xfm_atlas_dir+output_stem+'_filled'+ext, \
                                                    atlas_dir+file2+ext, '-overlap', str(label), \
                                                    '>'+results_dir+'temp_overlap.txt'])
                            p = call(args, shell="True")
                            f = open(results_dir+'temp_overlap.txt','r')
                            temp = f.read()
                            if temp != '':
                                dice = float(temp.split()[-2].split(',')[0])
                                jacc = float(temp.split()[-1].split(',')[0])
                            else:
                                dice = 0.0
                                jacc = 0.0
                            if isnan(dice):
                                dice = 0.0
                            if isnan(jacc):
                                jacc = 0.0
                            print_out = ' '.join(['Label:', str(label), 'Dice:', str(dice), \
                                                                     'Jaccard:', str(jacc)])
                            print(print_out)
                            f_eval.close()
                            f_eval = open(results_file, 'a')
                            f_eval.write(print_out + '\n')
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
                        f_eval.write(print_out1 + '\n' + print_out2 + '\n')
                        f_eval.close()
                        f_avg.close()
                        f_avg = open(avg_results_file, 'a');
                        f_avg.write(output_stem + ' ' + str(average_dice) + ' ' + str(average_jacc) + '\n')
                else:
                    if not os.path.exists(brain_dir+file+ext):
                        raise NameError('Check input file ' + brain_dir+file+ext)
                    elif not os.path.exists(brain_dir+file2+ext):
                        raise NameError('Check input file ' + brain_dir+file2+ext)
                    elif not os.path.exists(xfm_dir+file+'Warp.nii.gz'):
                        raise NameError('Check input file ' + xfm_dir+file+'Warp.nii.gz')
    if evaluate_with_atlases:
        f_avg.close()

#----------------------------------------------
# Apply landmark-driven registration transforms
# to register brains to each other via template
#----------------------------------------------
if transform_landmarks_via_template:
    if evaluate_with_atlases:
        avg_results_file = results_dir+'dice_jacc_overlaps_'+landmark_type+'.txt'
        f_avg = open(avg_results_file, 'w');
    for file in source_files:
        source = brain_dir+file+ext
        source_landmarks = landmarks_dir+file+ext
        for file2 in target_files:
            if file2 != file:
                target = brain_dir+file2+ext
                target_landmarks = landmarks_dir+file2+ext
                if os.path.exists(source) and \
                   os.path.exists(target) and \
                   os.path.exists(source_landmarks) and \
                   os.path.exists(target_landmarks):

                    pair = file+'_to_'+file2
                    inv_pair = file2+'_to_'+file
                    output_stem = pair+'_'+landmark_type
                    xfm_stem = xfm_dir+pair+'_in_template_space_'+landmark_type
                    inv_xfm_stem = xfm_dir+inv_pair+'_in_template_space_'+landmark_type
                    # Transform brains
                    if not os.path.exists(xfm_brain_dir+output_stem+ext):
	                    args = " ".join([ANTSPATH+'WarpImageMultiTransform', str(dim), \
	                            source, xfm_brain_dir+output_stem+ext, '-R', target, \
	                            '-i', inv_xfm_stem+'Affine.txt', \
	                            inv_xfm_stem+'InverseWarp.nii.gz', \
	                            xfm_stem+'Warp.nii.gz', \
	                            xfm_stem+'Affine.txt'])
	                    if verbose: print(args); print(''); p = call(args, shell="True")

                    # Transform landmarks
                    if not os.path.exists(xfm_landmarks_dir+output_stem+ext):
                        args = " ".join([ANTSPATH+'WarpImageMultiTransform', str(dim), \
	                            source_landmarks, xfm_landmarks_dir+output_stem+ext, '-R',target_landmarks, \
	                            '-i', inv_xfm_stem+'Affine.txt', \
	                            inv_xfm_stem+'InverseWarp.nii.gz', \
	                            xfm_stem+'Warp.nii.gz', \
	                            xfm_stem+'Affine.txt','--use-NN'])
                        if verbose: print(args); print(''); p = call(args, shell="True")

                    if evaluate_with_atlases:
                        if not os.path.exists(xfm_atlas_dir+output_stem+ext):
                          if not os.path.exists(results_dir+output_stem+'.txt'):

	                        # Transform atlases
	                        source_labels = atlas_dir+file+ext
	                        target_labels = atlas_dir+file2+ext
	                        args = " ".join([ANTSPATH+'WarpImageMultiTransform', str(dim), \
	                                source_labels, xfm_atlas_dir+output_stem+ext, '-R', target_labels, \
	                            '-i', inv_xfm_stem+'Affine.txt', \
	                            inv_xfm_stem+'InverseWarp.nii.gz', \
	                            xfm_stem+'Warp.nii.gz', \
	                            xfm_stem+'Affine.txt','--use-NN'])
	                        if verbose: print(args); print(''); p = call(args, shell="True")
	 
	                        # Fill target atlas mask with transformed source atlas labels
	                        args = " ".join(['ImageMath', str(dim), xfm_atlas_dir+output_stem+'_filled'+ext, \
	                                         'PropagateLabelsThroughMask', brainmask_dir+file2+ext, \
	                                                                       xfm_atlas_dir+output_stem+ext])
	                        if verbose: print(args); print(''); p = call(args, shell="True")
	
	                        # Measure overlap of target atlas and transformed source atlas labels
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
                else:
                    if not os.path.exists(source_landmarks):
                        raise NameError('Check input file ' + source_landmarks)
                    elif not os.path.exists(target_landmarks):
                        raise NameError('Check input file ' + target_landmarks)
    if evaluate_with_atlases:
        f_avg.close()
