#!/usr/bin/python
"""
Run registration commands

(c) 2011, @rno klein
"""

import os
from os.path import exists
from subprocess import call
from numpy import random, float, isnan

from measure_overlap import measure_overlap

#
# Inputs
#
register_to_template = 0
register_landmarks_to_template = 0
register_pairs_via_template = 0
register_landmarks_via_template = 1

fill_target_mask = 0
measure_overlap = 0

verbose = 1
dim = 3

#
# Files
#
source_files       = ['m1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m11','m12']
target_files       = source_files
ANTSPATH           = os.environ.get("ANTSPATH")
FSLPATH            = '/usr/local/fsl/bin/'
out_path           = '/hd2/Archive/registration_evaluation_2011_output/'
xfm_dir            = os.path.join( out_path, 'Transforms/')
xfm_brain_dir      = os.path.join( out_path, 'Transformed_Brains/')
xfm_landmark_dir   = os.path.join( out_path, 'Transformed_Landmarks/')
xfm_atlas_dir      = os.path.join( out_path, 'Transformed_Atlases/')
atlas_dir          = '/hd2/Brains/CUMC12/Atlases/'
brain_dir          = '/hd2/Brains/CUMC12/Brains/'
ext                = '.nii.gz'
template           = '/hd2/Brains/CUMC12/CUMC12template.nii.gz'
landmarks_dir      = '/hd2/Brains/CUMC12/Landmarks/ribbons_brain_visa_binary/'
landmark_type      = 'ribbons_brain_vis'
#landmarks_dir      = '/hd2/Brains/CUMC12/Landmarks/fundi_gang_li_binary/'
#landmark_type      = 'fundi_gang_li'
#landmarks_dir      = '/hd2/Brains/CUMC12/Landmarks/fundi_brain_visa_binary/'
#landmark_type      = 'fundi_brain_visa'

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
landmark_weight1 = 1.0
landmark_weight2 = 1.0
percent = 1.0 # real number: 1.0 = 100%
boundary = 0  # 0: not only boundaries
sigma = 10
neighbors = 1
matching_iter = 100000 # partial matching iterations

#------------------------------------------
# Register brains and landmarks to template
#------------------------------------------
if register_to_template + register_landmarks_to_template > 0:
    for file in source_files:
        source = brain_dir+file+ext
        output = xfm_dir+file+'_to_template'
        out = '-o ' + output+ext
        if os.path.exists(source) and os.path.exists(template) and os.path.exists(xfm_dir):

            # Intensity:
            if register_to_template:
                intensity = [template, source, intensity_weight, intensity_setting]
                intensity = "-m "+intensity_measure+"[" + ", ".join([str(s) for s in intensity]) + "]"
                args = " ".join([warp, regularize, intensity, out])
                if verbose: print(args); print(''); p = call(args, shell="True")

            # Landmarks:
            if register_landmarks_to_template:
                source_landmarks = landmarks_dir+file+ext
                output_landmarks = xfm_landmark_dir+file+'_to_template_'+landmark_type+ext
                try:
                    os.path.exists(source_landmarks) and os.path.exists(xfm_landmark_dir)
                except:
                    raise NameError('Check ' + source_landmarks + ' and ' + xfm_landmark_dir)
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

#----------------------------------------------------------------------------
# Apply registration transforms to register brains to each other via template
#----------------------------------------------------------------------------
if register_pairs_via_template:
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
                    args = ANTSPATH + 'WarpImageMultiTransform ' + str(dim) + ' ' + source + ' ' + \
                      xfm_brain_dir+output_stem+ext + \
                      ' -R ' + target + \
                      ' -i ' + xfm_dir+file2+'_to_templateAffine.txt ' + \
                      xfm_dir+file2+'_to_templateInverseWarp.nii.gz ' + \
                      xfm_dir+file+'_to_templateWarp.nii.gz ' + \
                      xfm_dir+file+'_to_templateAffine.txt --use-NN '
                    if verbose: print(args); print(''); p = call(args, shell="True")

                    # Transform atlases
                    args = ANTSPATH + 'WarpImageMultiTransform ' + str(dim) + ' ' + \
                      atlas_dir+file+ext + ' ' + \
                      xfm_atlas_dir+output_stem+ext + \
                      ' -R ' + atlas_dir+file2+ext + \
                      ' -i ' + xfm_dir+file2+'_to_templateAffine.txt ' + \
                      xfm_dir+file2+'_to_templateInverseWarp.nii.gz ' + \
                      xfm_dir+file+'_to_templateWarp.nii.gz ' + \
                      xfm_dir+file+'_to_templateAffine.txt --use-NN '
                    if verbose: print(args); print(''); p = call(args, shell="True")

                else:
                    if not os.path.exists(brain_dir+file+ext):
                        raise NameError('Check input file ' + brain_dir+file+ext)
                    elif not os.path.exists(brain_dir+file2+ext):
                        raise NameError('Check input file ' + brain_dir+file2+ext)
                    elif not os.path.exists(xfm_dir+file+'Warp.nii.gz'):
                        raise NameError('Check input file ' + xfm_dir+file+'Warp.nii.gz')

#----------------------------------------------------------------
# Register brains to each other via template, driven by landmarks
#----------------------------------------------------------------
if register_landmarks_via_template:
    for file in source_files:
        source = brain_dir+file+ext
        for file2 in target_files:
            if file2 != file:
                target = brain_dir+file2+ext
                if os.path.exists(source_landmarks) and \
                    os.path.exists(target_landmarks):

                    source_landmarks = landmarks_dir+file+ext
                    target_landmarks = landmarks_dir+file2+ext
                    output_stem = file + '_to_' + file2 + '_' + landmark_type
                    output = output_stem + ext

                    # Intensity similarity:
                    intensity = [target, source, intensity_weight, intensity_setting]
                    intensity = " -m "+intensity_measure+"[" + ", ".join([str(s) for s in intensity]) + "]"

                    # Landmark similarity:
                    lm_args1 = [target, source, target_landmarks, source_landmarks,
                                landmark_weight1, percent, sigma, boundary, neighbor, matching_iter]
                    landmarks1 = ", ".join([" -m PSE[" + ", ".join([str(s) for s in lm_args1]) + "]"])
                    lm_args2 = [target_landmarks, source_landmarks, landmark_weight2, 0]
                    landmarks2 = " ".join([" -m MSQ[" + ", ".join([str(s) for s in lm_args2]) + "]"])

                    #
                    # Run commands
                    #
                    args = " ".join([warp, -o, xfm_dir+output, \
                                     transform, regularize, intensity, landmarks1, landmarks2])
                    if verbose: print(args); print(''); p = call(args, shell="True")

                    args = " ".join([apply_warp, source, xfm_brain_dir+output, '-R ' + target, \
                                     xfm_dir+output_stem+'Warp'+ext, xfm_dir+output_stem+'Affine.txt'])
                    if verbose: print(args); print(''); p = call(args, shell="True")

                    args = " ".join([apply_warp, source_labels, xfm_atlas_dir+output, '-R ' + target_labels, \
                                     xfm_dir+output_stem+'Warp'+ext, xfm_dir+output_stem+'Affine.txt', '--use-NN'])
                    if verbose: print(args); print(''); p = call(args, shell="True")

                    args = " ".join([apply_warp, source_landmarks, xfm_landmark_dir+output, '-R ' + target, \
                                     xfm_dir+output_stem+'Warp'+ext, xfm_dir+output_stem+'Affine.txt', '--use-NN'])
                    if verbose: print(args); print(''); p = call(args, shell="True")

                else:
                  if not os.path.exists(source_landmarks):
                      raise NameError('Check input file ' + source_landmarks)
                  elif not os.path.exists(target_landmarks):
                      raise NameError('Check input file ' + target_landmarks)


"""


  results_file = results_dir+output_stem+"_"+regularizer+"_"+intensity_measure+"_"+landmark_measure1+"_"+landmark_measure2+".txt"
  if save_results and os.path.exists(results_file) and overwrite == 0:
      raise NameError('File already exists: '+results_file+' -- set "overwrite to 1."')
  else:
    if save_results:
      f_eval = open(results_file, 'w')
    print('Number of tests: ' + str(len(gradient_step_sizes)*len(regularizer_settings)*len(intensity_weights)*len(landmark_weights1)*len(landmark_weights2)*len(sigmas)*len(neighbors)))

                    output_stem = file + '_to_' + file2
                    output = outpath + output_stem


                    # Intensity similarity:
                    intensity = [target, source, intensity_weight, intensity_setting]
                    intensity = " -m "+intensity_measure+"[" + ", ".join([str(s) for s in intensity]) + "]"

                    # Landmark similarity:
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

"""