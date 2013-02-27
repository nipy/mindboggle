#!/usr/bin/env python
"""
Run Mindboggle on the Mindboggle-101 distribution.

Author:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2013  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""
import os
from mindboggle.utils.io_file import read_columns

command = 'python ../pipeline.py'
output_path = '/desk/output'
atlases_file = os.path.join(os.environ['MINDBOGGLE_DATA'], 'info', 'atlases101.txt')
atlases = read_columns(atlases_file, n_columns=1)[0]
atlas_strings = ['MMRR','OASIS','NKI-RS','NKI-TRT']

for atlas in atlases:

#   if atlas_strings[0] in atlas:
    if 'HLN' in atlas or 'Twins' in atlas or 'Afterthought' in atlas or 'Colin27' in atlas:
#    if 'Afterthought' in atlas:

        cmd = '{0} {1} {2}'.format(command, output_path, atlas)
        print(cmd); os.system(cmd)
