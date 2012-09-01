#!/usr/bin/python
"""
Extract fundus curves from surface mesh patches (folds).

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#==================
# Extract all fundi
#==================
def extract_fundi(folds, n_folds, neighbor_lists,
                  depth_file, mean_curvature_file, min_curvature_vector_file,
                  min_fold_size=50, min_distance=5, thr=0.5):
    """
    Extract all fundi.

    A fundus is a connected set of high-likelihood vertices in a surface mesh.

    Parameters
    ----------
    folds : label indices for folds: [#vertices x 1] numpy array
    n_folds :  #folds [int]
    neighbor_lists : list of lists of neighboring vertex indices (see return_arrays)
    depth_file : surface mesh file in VTK format with scalar values
    mean_curvature_file : surface mesh file in VTK format with scalar values
    min_curvature_vector_file : surface mesh file in VTK format with scalar values
    min_fold_size : minimum fold size
    min_distance :  minimum distance
    thr :  likelihood threshold

    Returns
    -------
    fundi :  numpy array of fundi

    """

    import numpy as np
    from time import time

    from extract.fundi_hmmf.compute_likelihood import compute_likelihood
    from extract.fundi_hmmf.find_points import find_anchors
    from extract.fundi_hmmf.connect_points import connect_points

    from utils.io_vtk import load_scalar, inside_faces

    # Convert folds array to a list of lists of vertex indices
    index_lists_folds = [np.where(folds == i)[0].tolist()
                         for i in range(1, n_folds+1)]

    # Load depth and curvature values from VTK and text files
    vertices, Faces, depths = load_scalar(depth_file, return_arrays=1)
    vertices, Faces, mean_curvatures = load_scalar(mean_curvature_file,
                                                   return_arrays=1)
    min_directions = np.loadtxt(min_curvature_vector_file)

    # For each fold...
    print("Extract a fundus from each of {0} folds...".format(n_folds))
    t1 = time()
    fundus_lists = []
    n_vertices = len(depths)
    Z = np.zeros(n_vertices)
    likelihoods = Z.copy()

    for i_fold, indices_fold in enumerate(index_lists_folds):

        print('  Fold {0} of {1}:'.format(i_fold + 1, n_folds))

        # Compute fundus likelihood values
        fold_likelihoods = compute_likelihood(depths[indices_fold],
                                              mean_curvatures[indices_fold])
        likelihoods[indices_fold] = fold_likelihoods

        # If the fold has enough high-likelihood vertices, continue
        likelihoods_thr = sum(fold_likelihoods > thr)
        print('    Computed fundus likelihood values: {0} > {1} (minimum: {2})'.
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
                print('    Connect {0} fundus points...'.format(n_anchors))
                t2 = time()
                likelihoods_fold = Z.copy()
                likelihoods_fold[indices_fold] = fold_likelihoods

                # Remove surface mesh faces whose three vertices
                # are not all in "indices_fold"
                faces_folds = inside_faces(Faces, indices_fold)

                H = connect_points(indices_anchors, faces_folds, indices_fold,
                                   likelihoods_fold, thr, neighbor_lists)
                fundus_lists.append(H.tolist())
                print('      ...Connected {0} fundus points ({1:.2f} seconds)'.
                      format(n_anchors, time() - t2))
            else:
                fundus_lists.append([])
        else:
            fundus_lists.append([])

    fundi = np.ones(n_vertices)
    for fundus in fundus_lists:
        if len(fundus) > 0:
            fundi += fundus

    print('  ...Extracted fundi ({0:.2f} seconds)'.format(time() - t1))

    return fundi #np.array(fundus_lists), likelihoods
                 #fundus_lists, likelihoods
