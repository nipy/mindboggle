#!/usr/bin/python
"""
Extract fundus curves from surface mesh patches (folds).

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

from compute_likelihood import compute_likelihood
from find_points import find_anchors
from connect_points import connect_points
from time import time

import sys
sys.path.append('/projects/Mindboggle/mindboggle/mindboggle/utils/')
import io_vtk


#==================
# Extract all fundi
#==================
def extract_fundi(index_lists_folds, n_folds, neighbor_lists,
                  vertices, faces, depths, mean_curvatures, min_directions,
                  min_fold_size=50, thr=0.5, min_distance=5):
    """
    Extract all fundi.

    A fundus is a connected set of high-likelihood vertices in a surface mesh.

    Inputs:
    ------
    vertices:  [#vertices x 3] numpy array
    faces:  vertices for polygons [#faces x 3] numpy array
    depths:  depth values [#vertices x 1] numpy array
    mean_curvatures:  mean curvature values [#vertices x 1] numpy array
    min_directions:  directions of minimum curvature [3 x #vertices] numpy array
    fraction_folds:  fraction of vertices considered to be folds
    min_fold_size:  minimum fold size from which to find a fundus
    thr:  likelihood threshold
    min_distance:  minimum distance

    Output:
    ------
    fundi:  numpy array of fundi
    fundus_lists:  list of #folds lists of vertex indices.
    likelihoods:  numpy array of likelihood values

    Calls:
    -----
    compute_likelihood()
    find_anchors()
    connect_points()

    """

    # For each fold...
    print("Extract a fundus from each of {} folds...".format(n_folds))
    t1 = time()
    fundus_lists = []
    n_vertices = len(depths)
    Z = np.zeros(n_vertices)
    likelihoods = Z.copy()

    for i_fold, indices_fold in enumerate(index_lists_folds):

        print('  Fold {} of {}:'.format(i_fold + 1, n_folds))

        # Compute fundus likelihood values
        fold_likelihoods = compute_likelihood(depths[indices_fold],
                                              mean_curvatures[indices_fold])
        likelihoods[indices_fold] = fold_likelihoods

        # If the fold has enough high-likelihood vertices, continue
        likelihoods_thr = sum(fold_likelihoods > thr)
        print('    Computed fundus likelihood values: {} > {} (minimum: {})'.
              format(likelihoods_thr, thr, min_fold_size))
        if likelihoods_thr > min_fold_size:

            # Find fundus points
            fold_indices_anchors = find_anchors(vertices[indices_fold, :],
                                                fold_likelihoods,
                                                min_directions[indices_fold],
                                                min_distance, thr)
            indices_anchors = [indices_fold[x] for x in fold_indices_anchors]
            n_anchors = len(indices_anchors)
            if n_anchors > 1:

                # Connect fundus points and extract fundus
                print('    Connect {} fundus points...'.format(n_anchors))
                t2 = time()
                likelihoods_fold = Z.copy()
                likelihoods_fold[indices_fold] = fold_likelihoods

                H = connect_points(indices_anchors, faces, indices_fold,
                                   likelihoods_fold, thr, neighbor_lists)
                fundus_lists.append(H.tolist())
                print('      ...Connected {} fundus points ({:.2f} seconds)'.
                      format(n_anchors, time() - t2))
            else:
                fundus_lists.append([])
        else:
            fundus_lists.append([])

    fundi = np.ones(n_vertices)
    for fundus in fundus_lists:
        if len(fundus) > 0:
            fundi += fundus

    print('  ...Extracted fundi ({:.2f} seconds)'.format(time() - t1))

    return fundi, fundus_lists, likelihoods
