"""
Create a registration template for FreeSurfer following:
http://surfer.nmr.mgh.harvard.edu/fswiki/SurfaceRegAndTemplates

The process of creating a template involves calculating a "mean" 
(across subjects) pattern of curvature-related values, along with the variance 
of these variables. This task is performed by program mris_make_template. 
However, mris_make_template expects the input subjects to be already aligned 
(ie: already have ?h.sphere.reg files). This of course will not be the case 
if you are about to create a template!

The process starts ("Prep") by choosing one reference subject, and producing 
a template (Here, mytemp0.tif). No ?h.sphere is required; this step uses the
subject's own pattern as the reference.

In Round 1, all reference subjects data are registered to the initial template, 
producing ?h.sphere.myreg0.  Using that registration, a new template
(mytemp1.tif) is made using all reference subjects as input.

Round 2 is essentially a repeat of Round 1, but produces an improved template 
(mytemp2.tif) because Round 2's initial registration step is based on a better 
template than in Round 1.  The steps of Round 2 could be iterated again 
if needed, but from what I've learned from Bruce Fischl, this is generally not 
necessary.

Notes

1. The website's diagram assumes that the subjects have all been processed 
once by FreeSurfer's recon-all pipeline or similar to produce all the required 
files needed for the template-creation process. Recon-all will also have 
produced ?h.sphere.reg aligned to the standard reference template 
(?h.average.curvature.filled.buckner40.tif), and corresponding annotation files 
-- for purposes here, these can be ignored. Indeed, once a satisfactory new 
template has been created, the recon-all script can be modified to use the new 
template, and the subjects can be re-run to create new improved annotations.

2. Comment: Since the process uses a particular subject as the starting point, 
the final position of the template pattern is dependent on that first subject. 
The iterations of registration and template creation work to improve the 
variance aspect of the template.

Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2011  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os, sys

make_initial_template = 1
make_template = 1
iterations = 2
subject_path = 'ManualSurfandVolLabels/subjects/'
save_path = 'templates_freesurfer/'
name = 'NKI_Rockland'
hemis = ['lh','rh']
subject_list_from_file = 1  # Get subject list from a file (1) or directory (0)
if subject_list_from_file:
    subject_list_file = 'subject_lists/NKI_Rockland.txt'

    f = open(subject_list_file)
    subject_lines = f.readlines()
    subject_list = []
    subject_string = ""
    for subject in subject_lines:
        subject_list.append(subject.strip("\n"))
        subject_string += subject.strip("\n") + " "
# Get list of subjects from a directory:
else:
    subject_list = os.listdir(subject_path)

# Make initial template:  mris_make_template [options] <hemi> <surface name> 
#                                            <subject> <subject>...<output name>
if make_initial_template:
    for hemi in hemis:
        template = save_path+hemi+'.'+name+'_0.tif'
        if os.path.exists(template):
            pass
        else:
            cmd = 'mris_make_template ' + hemi + ' sphere ' + subject_list[0] +\
                  ' ' + template
            print(cmd); os.system(cmd)

# Iteratively refine template:
if make_template:
    for i in range(0,iterations):
        for hemi in hemis:
            for subject in subject_list:
                src_file = subject_path + subject + '/surf/' + hemi + '.sphere'
                reg_file = src_file + '.reg' + str(i) + '_' + name
                template = save_path + hemi + '.' + name + '_' + str(i) + '.tif'
                if os.path.exists(reg_file):
                    pass
                else:
                    args = ['mris_register -curv', src_file, template, reg_file]
                    print(' '.join(args)); os.system(' '.join(args));

            template_new = save_path+hemi+'.'+name+'_'+str(i+1)+'.tif'
            if os.path.exists(template_new):
                pass
            else:
                args = ['mris_make_template', hemi, 'sphere.reg'+str(i)+'_'+name,
                        subject_string, template_new]
                print(' '.join(args)); os.system(' '.join(args)); # p = Popen(args);    

