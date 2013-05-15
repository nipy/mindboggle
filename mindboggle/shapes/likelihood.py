#!/usr/bin/env python
"""
Compute (fundus) likelihood values that highlight deep and highly curved
portions of a surface mesh.

Compute likelihood values for a VTK surface mesh:

  compute_likelihood()
    Compute likelihood values for a given VTK surface mesh.
    This is run after training on the distributions of depth and curvature
    values across multiple VTK surface mesh files in the functions below.

Learn distributions from training data (see Examples below):

  estimate_depth_curvature_distributions()
    Estimate distribution means, sigmas (standard deviations), and weights
    for VTK surface mesh depth and curvature scalars along and outside
    sulcus label borders within folds (as defined by a labeling protocol).

  which calls the functions:

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

def compute_likelihood(trained_file, depth_file, curvature_file, folds,
                       save_file=False):
    """
    Compute likelihoods based on input values, folds, and estimated parameters.

    Compute likelihood values for a given VTK surface mesh file, after training
    on distributions of depth and curvature values from multiple files.

    Parameters
    ----------
    trained_file : pickle compressed file
        contains the following dictionaries containing lists of floats
        (estimates of depth or curvature means, sigmas, and weights
         trained on fold vertices either on or off sulcus label borders)
        depth_border, curv_border, depth_nonborder, curv_nonborder
    depth_file : string
        VTK surface mesh file with depth values in [0,1] for all vertices
    curvature_file : string
        VTK surface mesh file with curvature values in [-1,1] for all vertices
    folds : list of integers
        fold number for all vertices (-1 for non-fold vertices)
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    likelihoods : list of floats
        likelihood values for all vertices (0 for non-fold vertices)
    likelihoods_file : string (if save_file)
        name of output VTK file with likelihood scalars
        (-1 for non-fold vertices)

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.shapes.likelihood import compute_likelihood
    >>> from mindboggle.utils.plots import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> trained_file = os.path.join(path, 'atlases', 'depth_curv_border_nonborder_parameters.pkl')
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    >>> curvature_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> save_file = True
    >>> #
    >>> compute_likelihood(trained_file, depth_file, curvature_file, folds, save_file)
    >>> # View:
    >>> plot_vtk('likelihoods.vtk', folds_file)

    """
    import os
    import numpy as np
    from math import pi
    import cPickle as pickle

    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars


    # Initialize variables:
    tiny = 0.000000001
    L = np.zeros(len(folds))
    probs_border = np.zeros(len(folds))
    probs_nonborder = np.zeros(len(folds))

    # Load estimated depth and curvature distribution parameters:
    depth_border, curv_border, depth_nonborder, curv_nonborder = pickle.load(
        open(trained_file, "r"))

    # Load depths, curvatures:
    depths, name = read_scalars(depth_file, True, True)
    curvatures, name = read_scalars(curvature_file, True, True)

    # Prep for below:
    n = 2
    twopiexp = (2*pi)**(n/2)
    border_sigmas = depth_border['sigmas'] * curv_border['sigmas']
    nonborder_sigmas = depth_nonborder['sigmas'] * curv_nonborder['sigmas']
    norm_border = 1 / (twopiexp * border_sigmas + tiny)
    norm_nonborder = 1 / (twopiexp * nonborder_sigmas + tiny)
    I = [i for i,x in enumerate(folds) if x != -1]

    N = depth_border['sigmas'].shape[0]
    for j in range(N):

        # Depth:
        expB = depth_border['weights'][j] * \
            ((depths[I]-depth_border['means'][j])**2) / \
            depth_border['sigmas'][j]**2
        expB += curv_border['weights'][j] * \
            ((curvatures[I]-curv_border['means'][j])**2) / \
            curv_border['sigmas'][j]**2
        expB = -expB / 2
        probs_border[I] = probs_border[I] + norm_border[j] * np.exp(expB)

        # Curvature:
        expNB = depth_nonborder['weights'][j] * \
            ((depths[I]-depth_nonborder['means'][j])**2) / \
            depth_nonborder['sigmas'][j]**2
        expNB += curv_nonborder['weights'][j] * \
            ((curvatures[I]-curv_nonborder['means'][j])**2) / \
            curv_nonborder['sigmas'][j]**2
        expNB = -expNB / 2
        probs_nonborder[I] = probs_nonborder[I] + norm_nonborder[j] * np.exp(expNB)

    likelihoods = probs_border / (probs_nonborder + probs_border + tiny)
    likelihoods.tolist()

    #-------------------------------------------------------------------------
    # Return likelihoods and output file name
    #-------------------------------------------------------------------------
    if save_file:

        likelihoods_file = os.path.join(os.getcwd(), 'likelihoods.vtk')
        rewrite_scalars(depth_file, likelihoods_file, likelihoods,
                        'likelihoods', likelihoods)
    else:
        likelihoods_file = None

    return likelihoods, likelihoods_file

#-------------------------------------------------------------------------------
# Learn distributions from training data (different surface meshes).
#-------------------------------------------------------------------------------

def estimate_distribution(scalar_files, scalar_range, fold_files, label_files):
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
    fold_files : list of strings
        names of VTK files with fold numbers for scalar values
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
    >>> import os
    >>> import numpy as np
    >>> from mindboggle.shapes.likelihood import estimate_distribution
    >>> from mindboggle.utils.io_file import read_columns
    >>> do_test = False
    >>> # Train on a single surface mesh:
    >>> if do_test:
    >>>     path = os.environ['MINDBOGGLE_DATA']
    >>>     depth_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    >>>     curv_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>>     folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>>     labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>>     depth_files = [depth_file]
    >>>     curv_files = [curv_file]
    >>>     fold_files = [folds_file]
    >>>     label_files = [labels_file]
    >>> # Train on many Mindboggle-101 surface meshes:
    >>> else:
    >>>     mindboggle_path = '/Users/arno/Data/Mindboggle-101/results'
    >>>     label_path = os.environ['SUBJECTS_DIR']
    >>>     x_path = os.path.join(os.environ['MINDBOGGLE'], 'x')
    >>>     atlas_list_file = os.path.join(x_path, 'mindboggle101_atlases.txt')
    >>>     atlas_list = read_columns(atlas_list_file, 1)[0]
    >>>     depth_files = []
    >>>     curv_files = []
    >>>     fold_files = []
    >>>     label_files = []
    >>>     for atlas in atlas_list:
    >>>      if 'OASIS' in atlas or 'NKI' in atlas or 'MMRR-21' in atlas:
    >>>       print(atlas)
    >>>       for h in ['lh','rh']:
    >>>         depth_file = os.path.join(mindboggle_path, 'shapes',
    >>>             '_hemi_'+h+'_subject_'+atlas, 'depth_rescaled.vtk')
    >>>         curv_file = os.path.join(mindboggle_path, 'shapes',
    >>>             '_hemi_'+h+'_subject_'+atlas, h+'.pial.mean_curvature.vtk')
    >>>         folds_file = os.path.join(mindboggle_path, 'features',
    >>>             '_hemi_'+h+'_subject_'+atlas, 'folds.vtk')
    >>>         labels_file = os.path.join(label_path, atlas, 'label',
    >>>             h+'.labels.DKT25.manual.vtk')
    >>>         depth_files.append(depth_file)
    >>>         curv_files.append(curv_file)
    >>>         fold_files.append(folds_file)
    >>>         label_files.append(labels_file)
    >>> scalar_range1 = np.linspace(0, 1, 51, endpoint=True) # (0 to 1 by 0.02)
    >>> scalar_range2 = np.linspace(-1, 1, 101, endpoint=True) # (-1 to 1 by 0.02)
    >>> #
    >>> depth_border, depth_nonborder = estimate_distribution(depth_files,
    >>>     scalar_range1, fold_files, label_files)
    >>> #
    >>> curv_border, curv_nonborder = estimate_distribution(curv_files,
    >>>     scalar_range2, fold_files, label_files)
    >>> #
    >>> import cPickle as pickle
    >>> pickle.dump( [depth_border, curv_border, depth_nonborder, curv_nonborder],
    >>>     open("depth_curv_border_nonborder_parameters.pkl", "wb"))

    """
    from mindboggle.shapes.likelihood import concatenate_sulcus_scalars, \
        fit_normals_to_histogram

    if not scalar_files or not fold_files or not label_files:
        import sys
        sys.exit("Input file lists cannot be empty.")

    # Concatenate scalars across multiple training files:
    border_scalars, nonborder_scalars = concatenate_sulcus_scalars(scalar_files,
        fold_files, label_files)

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

def concatenate_sulcus_scalars(scalar_files, fold_files, label_files):
    """
    Prepare data for estimating scalar distributions along and outside fundi.

    Extract (e.g., depth, curvature) scalar values in folds, along sulcus
    label boundaries as well as outside the sulcus label boundaries.
    Concatenate these scalar values across multiple files.

    Parameters
    ----------
    scalar_files : list of strings
        names of surface mesh VTK files with scalar values to concatenate
    fold_files : list of strings (corr. to each list in scalar_files)
        VTK files with fold numbers as scalars (-1 for non-fold vertices)
    label_files : list of strings (corr. to fold_files)
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
    >>> from mindboggle.shapes.likelihood import concatenate_sulcus_scalars
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> scalar_files = [depth_file, depth_file]
    >>> fold_files = [folds_file, folds_file]
    >>> label_files = [labels_file, labels_file]
    >>> #
    >>> S = concatenate_sulcus_scalars(scalar_files, fold_files, label_files)

    """
    import numpy as np

    from mindboggle.utils.io_vtk import read_scalars
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.labels.labels import extract_borders
    from mindboggle.labels.protocol.sulci_labelpairs_DKT import sulcus_boundaries


    # Prepare (non-unique) list of sulcus label pairs:
    protocol_label_pair_lists = sulcus_boundaries()
    protocol_label_pairs = [x for lst in protocol_label_pair_lists for x in lst]

    border_scalars = []
    nonborder_scalars = []

    # Loop through files with the scalar values:
    for ifile, scalar_file in enumerate(scalar_files):

        # Load scalars, folds, and labels:
        folds_file = fold_files[ifile]
        labels_file = label_files[ifile]
        scalars, name = read_scalars(scalar_file, True, True)
        if scalars.shape:
            folds, name = read_scalars(folds_file)
            labels, name = read_scalars(labels_file)
            indices_folds = [i for i,x in enumerate(folds) if x != -1]
            neighbor_lists = find_neighbors_from_file(labels_file)

            # Find all label border pairs within the folds:
            indices_label_pairs, label_pairs, unique_pairs = extract_borders(
                indices_folds, labels, neighbor_lists, ignore_values=[-1],
                return_label_pairs=True)
            indices_label_pairs = np.array(indices_label_pairs)

            # Find vertices with label pairs in the sulcus labeling protocol:
            Ipairs_in_protocol = [i for i,x in enumerate(label_pairs)
                                  if x in protocol_label_pairs]
            indices_label_pairs = indices_label_pairs[Ipairs_in_protocol]
            indices_outside_pairs = list(frozenset(indices_folds).difference(
                indices_label_pairs))

            # Store scalar values in folds along label border pairs:
            border_scalars.extend(scalars[indices_label_pairs].tolist())

            # Store scalar values in folds outside label border pairs:
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
    >>> from mindboggle.shapes.likelihood import fit_normals_to_histogram
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
    tiny = 0.000000001
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
            m1 = 1 / (sigmas[i] * np.sqrt(2*pi) + tiny)
            m2 = -((data-means[i])**2) / (2 *(sigmas[i]**2) + tiny)
            probs[:,i] = m1 * np.exp(m2)

        for i in range(k):
            W[:,i] = probs[:,i] / (np.sum(probs, axis=1) + tiny)

        for i in range(k):
            n1 = sum(W[:,i] * (data - means[i])**2)
            d1 = sum(W[:,i]) + tiny
            sigmas[i] =  np.sqrt(n1 / d1)
            means[i] = sum(W[:,i] * data) / d1

        print('    means: {0}; sigmas: {1}'.format(means, sigmas))

    for i in range(k):
        weights[i] = sum(W[:,i]) / (np.sum(W) + tiny)

    print('    weights: {0}'.format(weights))

    return means, sigmas, weights
