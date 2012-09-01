#!/usr/bin/python
"""
Compute likelihood and cost values for curves on a surface mesh.

Authors:
    - Yrjo Hame  (yrjo.hame@gmail.com)
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#=================
# Sigmoid function
#=================
def sigmoid(values, gain, shift):
    """
    Map values with sigmoid function to range [0,1].

    Y(t) = 1/(1 + exp(-gain*(values - shift))
    """
    import numpy as np

    # Make sure argument is a numpy array
    if type(values) != np.ndarray:
        values = np.array(values)

    return 1.0 / (1.0 + np.exp(-gain * (values - shift)))


#====================
# Likelihood function
#====================
def compute_likelihood(depths, curvatures):
    """
    Compute (fundus) curve likelihood values on a surface mesh.

    The *slope_factor* is used to compute the "gain" of the slope of sigmoidal values.

    Parameters
    ----------
    depths : depth values in [0,1]: [#sulcus vertices x 1] numpy array
    curvatures : mean curvature values in [-1,1] [#sulcus vertices x 1] array

    Returns
    -------
    likelihoods : likelihood values [#sulcus vertices x 1] numpy array

    """

    import numpy as np

    # Make sure arguments are numpy arrays
    if type(depths) != np.ndarray:
        depths = np.array(depths)
    if type(curvatures) != np.ndarray:
        curvatures = np.array(curvatures)

    #==========================================
    # Normalize depth values to interval [0,1]
    # Curvature values retain their values
    # Compute the means and std. deviations
    #=========================================
    depths_norm = depths / max(depths)
    depth_avg = np.mean(depths_norm)
    depth_std = np.std(depths_norm)
    curve_avg = np.mean(curvatures)
    curve_std = np.std(curvatures)

    #==========================================
    # Find slope for depth and curvature values
    #==========================================
    # Factor influencing "gain" or "sharpness" of the sigmoidal function below
    # slope_factor = abs(np.log((1. / x) - 1))  # 2.197224577 for x = 0.9
    # gain_depth = slope_factor / (2 * depth_std)
    gain_depth = 1 / depth_std
    gain_curve = 1 / curve_std
    shift_depth = depth_avg - depth_std
    shift_curve = curve_avg

    #==========================
    # Compute likelihood values
    #==========================
    # Map values with sigmoid function to range [0,1]
    depth_sigmoid = sigmoid(depths_norm, gain_depth, shift_depth)
    curve_sigmoid = sigmoid(curvatures, gain_curve, shift_curve)

    likelihoods = depth_sigmoid * curve_sigmoid

    # Plot the sigmoid curves (does not include value distributions)
    plot_result = False
    if plot_result:
        from matplotlib import pyplot
        xdepth = np.sort(depths_norm)
        xcurve = np.sort(curvatures)
        depth_sigmoid_sort = sigmoid(xdepth, gain_depth, shift_depth)
        curve_sigmoid_sort = sigmoid(xcurve, gain_curve, shift_curve)
        sigmoids = depth_sigmoid_sort * curve_sigmoid_sort
        pyplot.plot(xdepth, depth_sigmoid_sort, 'k')
        pyplot.plot(xcurve, curve_sigmoid_sort, 'b')
        pyplot.plot(xdepth, sigmoids, 'r')
        pyplot.title('Depths, curves: (gains={0:.2f},{1:.2f}; shifts={2:.2f},{3:.2f})'.
               format(gain_depth, gain_curve, shift_depth, shift_curve))
        pyplot.show()

    return likelihoods
