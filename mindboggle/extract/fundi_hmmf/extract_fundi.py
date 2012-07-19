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
from find_anchors import find_anchors
from connect_anchors import connect_anchors
from time import time

import sys
sys.path.append('/projects/mindboggle/mindboggle/mindboggle/utils/')
import io_vtk

#==================
# Extract all fundi
#==================
def extract_fundi(vertices, faces, depths, mean_curvatures, min_directions,
                  depth_threshold=0.2, thr=0.5, min_fold_size=50,
                  min_distance=5, max_distance=8):
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
    depth_threshold: depth threshold for defining folds
    thr: likelihood threshold
    min_fold_size: minimum fold size from which to find a fundus
    min_distance: minimum distance
    max_distance: maximum distance

    Output:
    ------
    fundi: list of #folds lists of vertex indices.

    Calls:
    -----
    extract_folds()
    compute_likelihood()
    find_anchors()
    connect_anchors()

    """

    import pickle
    load_path = "/drop/input/"
    load_em = 1
    save_em = 1

    # Extract folds (vertex indices for each fold)
    if load_em:
        folds = pickle.load(open(load_path + "folds.p","rb"))
        n_folds = int(pickle.load(open(load_path + "n_folds.p","rb")))
        index_lists_folds = pickle.load(open(load_path + "index_lists_folds.p","rb"))
    else:
        print("Extract folds from surface mesh...")
        t0 = time()
        folds, n_folds, index_lists_folds = extract_folds(faces, 
               depths, depth_threshold, min_fold_size)
        print('  ...Extracted folds in {:.2f} seconds'.
              format(time() - t0))
        if save_em:
            pickle.dump(folds, open(load_path + "folds.p","wb"))
            pickle.dump(n_folds, open(load_path + "n_folds.p","wb"))
            pickle.dump(index_lists_folds, open(load_path + "index_lists_folds.p","wb"))

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
    fundi = []
    n_vertices = len(depths)
    Z = np.zeros(n_vertices)
    likelihoods = Z.copy()
    for i_fold, indices_fold in enumerate(index_lists_folds):

        # Compute fundus likelihood values
        t0 = time()
        fold_likelihoods = compute_likelihood(depths[indices_fold],
                                              mean_curvatures[indices_fold])
        likelihoods[indices_fold] = fold_likelihoods
        print('  Computed fundus likelihood values for fold {} ({:.2f} seconds)'.
              format(i_fold + 1, time() - t0))

#       # If the fold has enough high-likelihood vertices, continue
        print('  {} likelihood values > threshold'.format(sum(fold_likelihoods > thr)))
#       if sum(fold_likelihoods > thr) > min_fold_size: ... else: fundi.append([])

        # Find fundus points
        t0 = time()
        fold_indices_anchors = find_anchors(vertices[indices_fold, :],
                                            fold_likelihoods,
                                            min_directions[indices_fold],
                                            thr, min_distance, max_distance)
        indices_anchors = [indices_fold[x] for x in fold_indices_anchors]
        print('  Found fundus points for fold {} ({:.2f} seconds)'.
              format(i_fold + 1, time() - t0))
        if len(indices_anchors) > 0:

            # Connect fundus points and extract fundus
            t0 = time()
            likelihoods_fold = Z.copy()
            likelihoods_fold[indices_fold] = fold_likelihoods
            fundi.append(
                  connect_anchors(indices_anchors, faces, indices_fold,
                                  likelihoods_fold, thr))
            print('  Connected fundus points for fold {} ({:.2f} seconds)'.
                  format(i_fold + 1, time() - t0))
        else:
            fundi.append([])

    if save_em:
        pickle.dump(likelihoods, open(load_path + "likelihoods.p","wb"))
        pickle.dump(fundi, open(load_path + "fundi.p","wb"))

        # Remove faces that do not contain three fold vertices
        indices_folds = [x for lst in index_lists_folds for x in lst]
        fs = frozenset(indices_folds)
        faces_folds = [lst for lst in faces if len(fs.intersection(lst)) == 3]
        faces_folds = np.reshape(np.ravel(faces_folds), (-1, 3))

        # Save fold likelihoods
        io_vtk.writeSulci(load_path + 'likelihoods.vtk', vertices,
                          indices_folds, faces_folds,
                          LUTs=[likelihoods],
                          LUTNames=['likelihoods'])

        # Save fundi
        fundi_for_vtk = -np.ones(n_vertices)
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
    depths = np.array([x/max(depths) for x in depths])
    mean_curvatures = np.array(mean_curvatures)
    min_directions = np.loadtxt(dir_file)
    if save_em:
        pickle.dump(vertices, open(load_path + "vertices.p","wb"))
        pickle.dump(faces, open(load_path + "faces.p","wb"))
        pickle.dump(depths, open(load_path + "depths.p","wb"))
        pickle.dump(mean_curvatures, open(load_path + "mean_curvatures.p","wb"))
        pickle.dump(min_directions, open(load_path + "min_directions.p","wb"))

fundi = extract_fundi(vertices, faces, depths, mean_curvatures, min_directions,
    depth_threshold=0.2, thr=0.5, min_fold_size=50,
    min_distance=5, max_distance=8)

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
