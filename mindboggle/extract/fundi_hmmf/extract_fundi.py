#!/usr/bin/python
"""
Extract fundus curves from surface mesh patches (folds).

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

from extract_folds import extract_folds
from compute_likelihood import compute_likelihood
from find_points import find_anchors
from connect_points import connect_points
from time import time

import sys
sys.path.append('/projects/Mindboggle/mindboggle/mindboggle/utils/')
import io_vtk
from percentile import percentile

save_fundi = 1
save_pickles = 0

#==================
# Extract all fundi
#==================
def extract_fundi(vertices, faces, depths, mean_curvatures, min_directions,
                  fraction_folds=0.5, min_fold_size=50, thr=0.5, min_distance=5):
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
    extract_folds()
    compute_likelihood()
    find_anchors()
    connect_points()

    """

    import pickle
    load_path = "/drop/input/"
    load_folds = True
    save_folds = True
    save_likelihoods = True #False
    save_fundi = True

    # Extract folds (vertex indices for each fold)
    if load_folds:
        n_folds = int(pickle.load(open(load_path + "n_folds.p","rb")))
        index_lists_folds = pickle.load(open(load_path + "index_lists_folds.p","rb"))
        neighbor_lists = pickle.load(open(load_path + "neighbor_lists.p","rb"))
    else:
        print("Extract folds from surface mesh...")
        t0 = time()

        # Compute the minimum depth threshold for defining folds by determining the
        # percentile of depth values for the fraction of vertices that are not folds.
        # For example, if we consider the shallowest one-third of vertices not to be
        # folds, we compute the depth percentile, and two-thirds of vertices would
        # have at least this depth value and would be considered folds.
        min_depth = percentile(np.sort(depths), 1 - fraction_folds,
                               key=lambda x:x)

        folds, n_folds, index_lists_folds, neighbor_lists = extract_folds(
            faces, depths, min_depth, min_fold_size)
        print('  ...Extracted folds greater than {:.2f} depth in {:.2f} seconds'.
              format(min_depth, time() - t0))

        if save_folds:
            pickle.dump(folds, open(load_path + "folds.p","wb"))
            pickle.dump(n_folds, open(load_path + "n_folds.p","wb"))
            pickle.dump(index_lists_folds, open(load_path + "index_lists_folds.p","wb"))
            pickle.dump(neighbor_lists, open(load_path + "neighbor_lists.p","wb"))

            indices_folds = [x for lst in index_lists_folds for x in lst]
            # Remove faces that do not contain three fold vertices
            fs = frozenset(indices_folds)
            faces_folds = [lst for lst in faces if len(fs.intersection(lst)) == 3]
            faces_folds = np.reshape(np.ravel(faces_folds), (-1, 3))
            print('  Reduced {} to {} faces.'.format(len(faces),
                                                     len(faces_folds)))
            # Save vtk file
            folds_for_vtk = folds.copy()
            folds_for_vtk[folds == 0] = -1
            LUTs = [[int(x) for x in folds_for_vtk]]
            LUT_names = ['fold'+str(i+1) for i in range(n_folds)]
            io_vtk.writeSulci(load_path + 'folds.vtk', vertices, indices_folds,
                              faces_folds, LUTs=LUTs, LUTNames=LUT_names)

    # For each fold...
    print("Extract a fundus from each of {} folds...".format(n_folds))
    t1 = time()
    fundi = []
    n_vertices = len(depths)
    Z = np.zeros(n_vertices)
    likelihoods = Z.copy()

    for i_fold, indices_fold in enumerate(index_lists_folds):
#      print('Only computing fold 17')
#      if i_fold == 17:
        print('  Fold {} of {}:'.format(i_fold + 1, n_folds))

        # Compute fundus likelihood values
        fold_likelihoods = compute_likelihood(depths[indices_fold],
                                              mean_curvatures[indices_fold])
        likelihoods[indices_fold] = fold_likelihoods

        if save_fundi:

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

    if save_pickles:
        # Save pickled data
        pickle.dump(likelihoods, open(load_path + "likelihoods.p","wb"))
        pickle.dump(fundi, open(load_path + "fundi.p","wb"))

    if save_likelihoods or save_fundi:

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

    return fundi

import pickle
load_path = "/drop/input/"
load_em = 1
save_em = 1
if load_em:
    vertices = pickle.load(open(load_path + "vertices.p","rb"))
    faces = pickle.load(open(load_path + "faces.p","rb"))
    depths = pickle.load(open(load_path + "depths.p","rb"))
    mean_curvatures = pickle.load(open(load_path + "mean_curvatures.p","rb"))
    min_directions = pickle.load(open(load_path + "min_directions.p","rb"))
else:
    depth_file = load_path + 'lh.pial.depth.vtk'
    curv_file = load_path + 'lh.pial.curv.avg.vtk'
    dir_file = load_path + 'lh.pial.curv.min.dir.csv'
    vertices, faces, depths = io_vtk.load_VTK_Map(depth_file)
    vertices, faces, mean_curvatures = io_vtk.load_VTK_Map(curv_file)
    vertices = np.array(vertices)
    faces = np.array(faces)
    depths = np.array(depths)
    mean_curvatures = np.array(mean_curvatures)
    min_directions = np.loadtxt(dir_file)
    if save_em:
        pickle.dump(vertices, open(load_path + "vertices.p","wb"))
        pickle.dump(faces, open(load_path + "faces.p","wb"))
        pickle.dump(depths, open(load_path + "depths.p","wb"))
        pickle.dump(mean_curvatures, open(load_path + "mean_curvatures.p","wb"))
        pickle.dump(min_directions, open(load_path + "min_directions.p","wb"))

fundi = extract_fundi(vertices, faces, depths, mean_curvatures, min_directions,
                      fraction_folds=0.5, min_fold_size=50, thr=0.5, min_distance=5)

"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x, y, z = np.transpose(np.reshape([x for lst in p for x in lst], (-1, 3)))
ax.scatter(x, y, z)

v = np.transpose(np.reshape([x for lst in vertices for x in lst], (-1, 3)))
ax.scatter(v[0], v[1], v[2])

plt.show()
"""
