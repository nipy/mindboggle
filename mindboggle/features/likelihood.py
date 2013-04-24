#!/usr/bin/env python
"""
Compute fundus likelihood values.

Compute likelihood values for any surface mesh:

  compute_likelihood()
    Compute likelihood values for a given subject after training on
    distributions of depth and curvature values across VTK surface mesh files.

Learn distributions from training data (different surface meshes):

  estimate_depth_curvature_distributions()
    Estimate distribution means, sigmas (standard deviations), and weights for
    VTK surface mesh depth and curvature scalars along and outside sulcus label
    borders within folds.

  concatenate_sulcus_scalars()
    Prepare data for estimating scalar distributions along and outside fundi.
    Extract (e.g., depth, curvature) scalar values in folds, along sulcus
    label boundaries as well as outside the sulcus label boundaries.
    Concatenate these scalar values across multiple files.

  fit_normals_to_histogram()
    This Estimation-Maximization method returns estimated means, sigmas
    (standard deviations) and weights from a distribution of values.


Authors:
Yrjo Hame, 2012-2013  .  yrjo.hame@gmail.com
Arno Klein, 2012-2013  .  arno@mindboggle.info  .  www.binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#-------------------------------------------------------------------------------
# Compute likelihood values for any surface mesh based on learned parameters.
#-------------------------------------------------------------------------------

def compute_likelihood(trained_file, depth_file, curvature_file, sulci):
    """
    Compute likelihoods based on input values, sulci, and estimated parameters.

    Compute likelihood values for a given VTK surface mesh file, after training
    on distributions of depth and curvature values from multiple files.

    Parameters
    ----------
    trained_file : pickle compressed file
        contains the following dictionaries containing lists of floats
        (estimates of depth or curvature means, sigmas, and weights
         trained on sulcus vertices either on or off sulcus label borders)
        depth_border, curv_border, depth_nonborder, curv_nonborder
    depth_file : string
        name of VTK surface mesh file with depth values in [0,1] for all vertices
    curvature_file : string
        name of VTK surface mesh file with curvature values in [-1,1] for all vertices
    sulci : list of integers
        sulcus number for all vertices (-1 for non-sulcus vertices)

    Returns
    -------
    L : list of floats
        likelihood values for all vertices (0 for non-sulcus vertices)

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.features.likelihood import compute_likelihood
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    >>> curvature_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.avg.vtk')
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> trained_file = os.path.join(path, 'depth_curv_border_nonborder_parameters.pkl')
    >>> sulci, name = read_scalars(sulci_file)
    >>> #
    >>> L = compute_likelihood(trained_file, depth_file, curvature_file, sulci)
    >>> # View:
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> rewrite_scalars(depth_file, 'compute_likelihood.vtk', L, 'likelihoods', folds)
    >>> from mindboggle.utils.plots import plot_vtk
    >>> plot_vtk('compute_likelihood.vtk')

    """
    import numpy as np
    from math import pi
    import cPickle as pickle

    from mindboggle.utils.io_vtk import read_scalars


    # Initialize variables:
    tiny = 0.000000001
    L = np.zeros(len(sulci))
    probs_border = np.zeros(len(sulci))
    probs_nonborder = np.zeros(len(sulci))

    # Load estimated depth and curvature distribution parameters:
    depth_border, curv_border, depth_nonborder, curv_nonborder = pickle.load(
        open(trained_file, "r" ) )

    # Load depths, curvatures:
    depths, name = read_scalars(depth_file, True, True)
    curvatures, name = read_scalars(curvature_file, True, True)

    # Prep for below:
    k = 2
    twopiexp = (2*pi)**(k/2)
    norm_border = 1 / (twopiexp * depth_border['sigmas'] * curv_border['sigmas'])
    norm_nonborder = 1 / (twopiexp * depth_nonborder['sigmas'] * curv_nonborder['sigmas'])
    I = [i for i,x in enumerate(sulci) if x != -1]

    N = depth_border['sigmas'].shape[0]
    for j in range(N):

        # Depth:
        expB = depth_border['weights'][j] * \
               ((depths[I]-depth_border['means'][j])**2) / \
               depth_border['sigmas'][j]**2
        expB = expB + curv_border['weights'][j] * \
               ((curvatures[I]-curv_border['means'][j])**2) / \
               curv_border['sigmas'][j]**2
        expB = -expB / 2
        probs_border[I] = probs_border[I] + norm_border[j] * np.exp(expB)

        # Curvature:
        expNB = depth_nonborder['weights'][j] * \
               ((depths[I]-depth_nonborder['means'][j])**2) / \
               depth_nonborder['sigmas'][j]**2
        expNB = expNB + curv_nonborder['weights'][j] * \
                      ((curvatures[I]-curv_nonborder['means'][j])**2) / \
                      curv_nonborder['sigmas'][j]**2
        expNB = -expNB / 2
        probs_nonborder[I] = probs_nonborder[I] + norm_nonborder[j] * np.exp(expNB)

    L = probs_border / (probs_nonborder + probs_border + tiny)

    return L

#-------------------------------------------------------------------------------
# Learn distributions from training data (different surface meshes).
#-------------------------------------------------------------------------------

def estimate_distribution(scalar_files, scalar_range, sulcus_files, label_files):
    """
    Estimate sulcus label border scalar distributions from VTK files.

    Estimate distribution means, sigmas (standard deviations), and weights
    for VTK surface mesh scalars (e.g., depth, curvature) along and outside
    sulcus label borders within folds.

    Note : The number of classes, k, is currently hard-coded.

    Parameters
    ----------
    scalar_files : list of strings
        names of VTK files with scalar values for all surface vertices
    scalar_range : list of floats
        range of values to estimate distribution
    sulcus_files : list of strings
        names of VTK files with sulcus numbers for scalar values
    label_files : list of strings
        names of VTK files with label numbers for scalar values

    Returns
    -------
    border_parameters : dictionary containing lists of floats
        means, sigmas, weights
    nonborder_parameters : dictionary containing lists of floats
        means, sigmas, weights

    Examples
    --------
    >>> # Train on a single surface mesh:
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.features.likelihood import estimate_distribution
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    >>> curv_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.avg.vtk')
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> depth_files = [depth_file]
    >>> curv_files = [curv_file]
    >>> sulcus_files = [sulci_file]
    >>> label_files = [labels_file]
    >>> scalar_range1 = np.linspace(0, 1, 51, endpoint=True) # (0 to 1 by 0.02)
    >>> scalar_range2 = np.linspace(-1, 1, 101, endpoint=True) # (-1 to 1 by 0.02)
    >>> #
    >>> depth_border, depth_nonborder = estimate_distribution(depth_files,
    >>>     scalar_range1, sulcus_files, label_files)
    >>> #
    >>> curv_border, curv_nonborder = estimate_distribution(curv_files,
    >>>     scalar_range2, sulcus_files, label_files)
    >>> #
    >>> import cPickle as pickle
    >>> pickle.dump( [depth_border, curv_border, depth_nonborder, curv_nonborder],
    >>>     open("depth_curv_border_nonborder_parameters.pkl", "wb" ) )

    """
    from mindboggle.features.likelihood import concatenate_sulcus_scalars, \
        fit_normals_to_histogram


    # Concatenate scalars across multiple training files:
    border_scalars, nonborder_scalars = concatenate_sulcus_scalars(scalar_files,
        sulcus_files, label_files)


    if min(scalar_range) < -0.9:
        import numpy as np
        print('RANDOMIZE ZERO CURVATURE')
        border_scalars = np.array(border_scalars)
        nonborder_scalars = np.array(nonborder_scalars)
        inds = np.where(border_scalars == 0)[0]
        border_scalars[inds] = 0.05 * np.random.randn(len(inds))
        inds = np.where(nonborder_scalars == 0)[0]
        nonborder_scalars[inds] = 0.05 * np.random.randn(len(inds))
        border_scalars.tolist()
        nonborder_scalars.tolist()


    # Estimate distribution parameters:
    border_means, border_sigmas, \
        border_weights = fit_normals_to_histogram(border_scalars, scalar_range)
    nonborder_means, nonborder_sigmas, \
        nonborder_weights = fit_normals_to_histogram(nonborder_scalars, scalar_range)

    # Store outputs in dictionaries:
    border_parameters = {
        'means': border_means,
        'sigmas': border_sigmas,
        'weights': border_weights
    }
    nonborder_parameters = {
        'means': nonborder_means,
        'sigmas': nonborder_sigmas,
        'weights': nonborder_weights
    }

    return border_parameters, nonborder_parameters

def concatenate_sulcus_scalars(scalar_files, sulcus_files, label_files):
    """
    Prepare data for estimating scalar distributions along and outside fundi.

    Extract (e.g., depth, curvature) scalar values in folds, along sulcus
    label boundaries as well as outside the sulcus label boundaries.
    Concatenate these scalar values across multiple files.

    Parameters
    ----------
    scalar_files : list of strings
        names of surface mesh VTK files with scalar values to concatenate
    sulcus_files : list of strings (corr. to each list in scalar_file_lists)
        VTK files with sulcus numbers as scalars (-1 for non-sulcus vertices)
    label_files : list of strings (corr. to sulcus_files)
        VTK files with label numbers (-1 for unlabeled vertices)

    Returns
    -------
    border_scalars : list of floats
        concatenated scalar values within folds along sulcus label boundaries
    nonborder_scalars : list of floats
        concatenated scalar values within folds outside sulcus label boundaries

    Examples
    --------
    >>> # Concatenate (duplicate) depth scalars:
    >>> import os
    >>> from mindboggle.features.likelihood import concatenate_sulcus_scalars
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> scalar_files = [depth_file, depth_file]
    >>> sulcus_files = [sulci_file, sulci_file]
    >>> label_files = [labels_file, labels_file]
    >>> #
    >>> S = concatenate_sulcus_scalars(scalar_files, sulcus_files, label_files)

    """
    import numpy as np

    from mindboggle.utils.io_vtk import read_scalars
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.labels.label import extract_borders
    from mindboggle.labels.protocol.sulci_labelpairs_DKT import sulcus_boundaries


    # Prepare (non-unique) list of sulcus label pairs:
    protocol_label_pair_lists = sulcus_boundaries()
    protocol_label_pairs = [x for lst in protocol_label_pair_lists for x in lst]

    border_scalars = []
    nonborder_scalars = []

    # Loop through files of each type of scalar value:
    for ifile, scalar_file in enumerate(scalar_files):

        # Load scalars, sulci, and labels:
        sulci_file = sulcus_files[ifile]
        labels_file = label_files[ifile]
        scalars, name = read_scalars(scalar_file, True, True)
        if scalars.shape:
            sulci, name = read_scalars(sulci_file)
            labels, name = read_scalars(labels_file)
            indices_sulci = [i for i,x in enumerate(sulci) if x != -1]
            neighbor_lists = find_neighbors_from_file(labels_file)

            # Find all label border pairs within the sulcus folds:
            indices_label_pairs, label_pairs, unique_pairs = extract_borders(
                indices_sulci, labels, neighbor_lists, ignore_values=[-1],
                return_label_pairs=True)
            indices_label_pairs = np.array(indices_label_pairs)

            # Find vertices with label pairs in the sulcus labeling protocol:
            Ipairs_in_protocol = [i for i,x in enumerate(label_pairs)
                                  if x in protocol_label_pairs]
            indices_label_pairs = indices_label_pairs[Ipairs_in_protocol]
            #indices_outside_pairs = [x for x in indices_sulci
            #                         if x not in indices_label_pairs]
            indices_outside_pairs = list(frozenset(indices_sulci).difference(
                indices_label_pairs))

            # Store scalar values in sulci along label border pairs:
            border_scalars.extend(scalars[indices_label_pairs].tolist())

            # Store scalar values in sulci outside label border pairs:
            nonborder_scalars.extend(scalars[indices_outside_pairs].tolist())

    return border_scalars, nonborder_scalars

def fit_normals_to_histogram(data, x):
    """
    This Estimation-Maximization method returns estimated means, sigmas
    (standard deviations) and weights, each of length k (number of classes).

    Parameters
    ----------
    data : list of floats
        data to estimate distribution means, sigmas, and weights
    x : list of floats
        range of values used to initialize distribution means and sigmas

    Returns
    -------
    means : list of floats
        estimated mean for each class
    sigmas : list of floats
        estimated standard deviation for each class
    weights : list of floats
        weight for each class

    Examples
    --------
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.features.likelihood import fit_normals_to_histogram
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    >>> scalars, name = read_scalars(depth_file)
    >>> x = np.linspace(0, 1, 51, endpoint=True)
    >>> #
    >>> means, sigmas, weights = fit_normals_to_histogram(scalars, x)

    """
    import numpy as np
    from math import pi

    # Initialize variables:
    k = 3
    probs = np.zeros((len(data), k))
    W = np.zeros((len(data), k))
    means = np.zeros(k)
    sigmas = np.zeros(k)
    weights = np.zeros(k)

    # Initialize distribution means and sigmas:
    rangex = max(x) - min(x)
    for i in range(1, k + 1):
        means[i-1] = max(x) - rangex/2 - 0.2 * rangex * (i - k/2)
        sigmas[i-1] = 0.2

    print('Fitting normals to histograms...')

    # Iteratively compute probabilities, weights, means and sigmas:
    iter = 0
    while iter < 25:
        iter += 1

        for i in range(k):
            probs[:,i] = 1 / (sigmas[i]*np.sqrt(2*pi)) * \
                         np.exp((-1/(2*(sigmas[i]**2))) * (data-means[i])**2)

        for i in range(k):
            W[:,i] = probs[:,i] / np.sum(probs, axis=1)

        for i in range(k):
            sigmas[i] = np.sqrt(sum(W[:,i]*(data - means[i])**2) / sum(W[:,i]))
            means[i] = sum(W[:,i] * data) / sum(W[:,i])

        print('    means: {0}; sigmas: {1}'.format(means, sigmas))

    for i in range(k):
        weights[i] = sum(W[:,i]) / np.sum(W)

    print('    weights: {0}'.format(weights))

    return means, sigmas, weights
