#!/usr/bin/env python
"""
Miscellaneous computations.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def sigmoid(values, gain, shift):
    """
    Map values with sigmoid function to range [0,1].

    Y(t) = 1/(1 + exp(-gain*(values - shift))
    """
    import numpy as np

    tiny = 0.000000001

    # Make sure argument is a numpy array
    if type(values) != np.ndarray:
        values = np.array(values)

    return 1.0 / (1.0 + np.exp(-gain * (values - shift)) + tiny)
