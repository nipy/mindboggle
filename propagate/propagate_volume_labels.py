#!/usr/bin/python
"""

Command: to do

Example: to do

Brian Avants: "the initial box labels are propagated through the gray matter 
with gm-probability dependent speed. it uses the fast marching algorithm.   
you can control how tightly the propagation follows the gray matter label 
by adjusting the speed image -- e.g. a binary speed image will constrain 
the propagated label only to the gm."

Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2011  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import sys
import os

reg_method = 'ANTS'
path_in    = './Transformed_Atlases/'+reg_method+'/'
path_out   = './'+reg_method+'_freesurferGM/'
path_mask  = './native_freesurfer_labels_cortex/' 

subjects = ['S01','S02','S03']
subjects1 = subjects
subjects2 = subjects

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
