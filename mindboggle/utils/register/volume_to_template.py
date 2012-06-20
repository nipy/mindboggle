#!/usr/bin/python
"""
Use ANTs to register a subject's brain image volume to a template.

Command: python <this file name> <subject brain image volume file>
                                 <average template brain image volume file>
                                 <output path>

Example: python volume_to_template.py
                subject_brain.nii.gz template_brain.nii.gz ./output

This program uses ANTs SyN registration method.

Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2011  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os
import sys

template = '../templates/volume_templates/KKI_template.nii.gz'

# Check inputs
if len(sys.argv) == 4:
    subject = sys.argv[1] + '/'
    template = sys.argv[2] + '/'
    output_path = sys.argv[3]
else:
    print("Please check your command-line arguments.")
    sys.exit()

if os.path.exists(subject) and os.path.exists(template) and \
   os.path.exists(output_path):
    pass
else:
    print(subject +', ' + template +', or ' + output_path + " doesn't exist.")
    sys.exit()

output = output_path + 'registered_to_template_.nii.gz',
args = ['ANTS 3 -m', 'CC[' + template + ',' + subject + ',1,2]',
        '-r Gauss[2,0] -t SyN[0.5] -i 30x99x11',
        '-o', output, '--use-Histogram-Matching',
        '--number-of-affine-iterations 10000x10000x10000x10000x10000']

cmd = 'date'; os.system(cmd)
print(' '.join(args)); os.system(' '.join(args)); # p = Popen(args);
cmd = 'date'; os.system(cmd)
