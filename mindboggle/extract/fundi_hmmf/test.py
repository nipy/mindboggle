#!/usr/bin/python
"""
Test

Authors:
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

def test_fundi_hmmf():

    import os
    import re
    import numpy as np

    test_data_path = '/Users/arno/Dropbox/yrjo_code_io_data/input'
    inputs = ['mean_curvatures', 'depths', 'vertices', 'faces', 'min_directions']

    outputs = []
    for input_type in inputs:
        print(input_type)
        file = os.path.join(test_data_path, input_type + '.dat')
        if input_type in ['mean_curvatures', 'depths']:
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

        outputs.append(input_array)

    return outputs