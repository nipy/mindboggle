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
    vertices: [#vertices x 3] numpy array
    faces: vertices for polygons [#faces x 3] numpy array
    depths: depth values [#vertices x 1] numpy array
    mean_curvatures: mean curvature values [#vertices x 1] numpy array
    min_directions: directions of minimum curvature [3 x #vertices] numpy array
    fraction_folds: fraction of vertices considered to be folds
    min_fold_size: minimum fold size from which to find a fundus
    thr: likelihood threshold
    min_distance: minimum distance

    Output:
    ------
    fundi: list of #folds lists of vertex indices.

    Calls:
    -----
    compute_likelihood()
    find_anchors()
    connect_points()

    """

    # For each fold...
    print("Extract a fundus from each of {} folds...".format(n_folds))
    t1 = time()
    fundi = []
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
                fundi.append(H.tolist())
                print('      ...Connected {} fundus points ({:.2f} seconds)'.
                      format(n_anchors, time() - t2))
            else:
                fundi.append([])
        else:
            fundi.append([])

    print('  ...Extracted fundi ({:.2f} seconds)'.format(time() - t1))

"""
    # Remove faces that do not contain three fold vertices
    indices_folds = [x for lst in index_lists_folds for x in lst]
    fs = frozenset(indices_folds)
    faces_folds = [lst for lst in faces if len(fs.intersection(lst)) == 3]
    faces_folds = np.reshape(np.ravel(faces_folds), (-1, 3))

    # Save fold likelihoods
    if save_likelihoods:
        io_vtk.writeSulci(load_path + 'likelihoods.vtk', vertices,
                          indices_folds, faces_folds,
                          LUTs=[likelihoods],
                          LUTNames=['likelihoods'])

    # Save fundi
    if save_fundi:
        fundi_for_vtk = np.ones(n_vertices)
        for fundus in fundi:
            if len(fundus) > 0:
                fundi_for_vtk += fundus
        io_vtk.writeSulci(load_path + 'fundi.vtk', vertices,
            indices_folds, faces_folds,
            LUTs=[fundi_for_vtk], LUTNames=['fundi'])
"""
    return fundi
