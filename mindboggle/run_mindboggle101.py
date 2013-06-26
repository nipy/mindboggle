# Script for running Mindboggle on Mindboggle-101 set
import os
from mindboggle.utils.io_table import read_columns

out_path = '/homedir/Data/Mindboggle-101/'
x_path = os.path.join(os.environ['MINDBOGGLE'], 'x')
atlas_list_file = os.path.join(x_path, 'mindboggle101_atlases.txt')
atlas_list = read_columns(atlas_list_file, 1)[0]

for atlas in atlas_list:
<<<<<<< HEAD
    #if 'HLN' in atlas or 'Twins' in atlas or
    #   'Colin' in atlas or 'After' in atlas or
    #   'MMRR-3T7T' in atlas:
    #if 'MMRR-21' in atlas:
    #if 'OASIS-TRT' in atlas:
    if 'NKI-TRT' in atlas:
    #if 'NKI-RS' in atlas:
        cmd = ' '.join(['python pipeline.py', out_path, atlas])
        print(cmd); os.system(cmd)
=======
    #if 'HLN' in atlas or 'Twins' in atlas or\
    #   'Colin' in atlas or 'After' in atlas or\
    #   'MMRR-3T7T' in atlas:
    #if 'MMRR-21' in atlas:
    #if 'OASIS-TRT' in atlas:
    #if 'NKI-TRT' in atlas:
    if 'NKI-RS' in atlas:
        cmd = ' '.join(['python pipeline.py', out_path, atlas])
        print(cmd); os.system(cmd)
>>>>>>> b62c6068b9919b4d6d993aeca1650ea2db52ea63
