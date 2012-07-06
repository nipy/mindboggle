#!/usr/bin/python
"""
Test

Authors:
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

def test_fundi_hmmf(debug = 0):

    import os
    import re
    import numpy as np

    test_data_path = '/home/arno/Dropbox/yrjo_code_io_data'
    inputs = ['mean_curvatures', 'depths', 'vertices', 'faces', 'min_directions',
              'output_sulci', 'output_anchor_points', 'output_L',
              'output_fundi']

    if debug:
        n_vertices = 5000

    outputs = []
    for input_type in inputs:
        print(input_type)
        file = os.path.join(test_data_path, input_type + '.dat')
        if input_type in ['mean_curvatures', 'depths', 'output_sulci', 'output_L']:
            input_array = np.loadtxt(file)
        else:
            f = open(file, 'r')
            input_string = f.read()
            f.close()
            input_split = re.split('[,;\n]', input_string)
            input_list = [float(s) for s in input_split if s!='']
            if input_type in ['faces', 'min_directions']:
                input_list = [int(s) for s in input_list]
                input_array = np.reshape(input_list, (-1, 3))
            elif input_type == 'vertices':
                input_array = np.reshape(input_list, (-1, 3))
            elif input_type in ['output_anchor_points', 'output_fundi']:
                input_array = np.reshape(input_list, (-1, 1))

        if debug:
            if input_type == 'faces':
                i_ignore = np.unique(np.where(input_array >= n_vertices)[0])
                i_use = np.setxor1d(range(len(input_array)), i_ignore)
                input_array = input_array[i_use, :]
            else:
                input_array = input_array[range(n_vertices), :]

        outputs.append(input_array)

    return outputs
