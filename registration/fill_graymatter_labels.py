#!/usr/bin/python
# brian avants: "the initial box labels are propagated through the gray matter with 
# gm-probability dependent speed. it uses the fast marching algorithm.   
# you can control how tightly the propagation follows the gray matter label 
# by adjusting the speed image -- e.g. a binary speed image will constrain 
# the propagated label only to the gm."

"""

(c) 2010, @rno klein
"""
import sys
import os

#reg_method = 'ANTS_via_MNI152'
#reg_method = 'ANTS_via_templates'
#reg_method = 'ANTS_via_Colin'
#reg_method = 'ANTS_via_LPBA40_FLIRT'
#reg_method = 'ANTS_via_LPBA40_AIR'
#reg_method = 'ANTS_via_LPBA40_SPM'
#reg_method = 'ANTS_via_LPBA40_ANTS'
reg_method = 'ANTS_via_LPBA40_ANTS2'
#reg_method = 'ANTS1.9'
path_in    = '/hd2/data/Projects/registration_evaluation_2010_output/LPBA40/Transformed_Atlases/'+reg_method+'/'
path_out   = '/hd2/data/Projects/registration_evaluation_2010_output/LPBA40/Transformed_Atlases/'+reg_method+'_freesurferGM/'
path_mask  = '/hd2/data/Atlases/LONI/LPBA40_renamed/native_freesurfer_labels_cortex/'    # ANTS and ANTS_via_templates

subjects = ['S01','S02','S03','S04','S05','S06','S07','S08','S09','S10',
            'S11','S12','S13','S14','S15','S16','S17','S18','S19','S20',
            'S21','S22','S23','S24','S25','S26','S27','S28','S29','S30',
            'S31','S32','S33','S34','S35','S36','S37','S38','S39','S40']
subjects1 = subjects
subjects2 = subjects
#subjects1 = ['S01','S02','S03','S04','S05','S06','S07','S08','S09','S10']
#subjects1 = ['S11','S12','S13','S14','S15','S16','S17','S18','S19','S20']
#subjects1 = ['S21','S22','S23','S24','S25','S26','S27','S28','S29','S30']
#subjects1 = ['S31','S32','S33','S34','S35','S36','S37','S38','S39','S40']
#subjects1.reverse()

for subject1 in subjects1:
    for subject2 in subjects2:

        if subject1 != subject2:

            input_file = path_in+subject1+'_to_'+subject2+'.nii'
            output_file = path_out+subject1+'_to_'+subject2+'.nii'
            propagation_mask_file = path_mask+subject2+'.nii'

            if os.path.isfile(output_file):
              pass
            else:
              print('No '+output_file)
              cmd = 'ImageMath 3 ' + output_file + ' PropagateLabelsThroughMask ' + propagation_mask_file + ' ' + input_file
              print cmd; os.system(cmd)
