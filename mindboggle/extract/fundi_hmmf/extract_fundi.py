#!/usr/bin/python
"""
Extract fundus curves from surface mesh patches (sulci).

Authors:
Yrjo Hame  .  yrjo.hame@gmail.com
Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import numpy as np

from extract_sulci import extract_sulci
from compute_likelihood import compute_likelihood
from find_anchors import find_anchors
from connect_anchors import connect_anchors
from time import time

#==================
# Extract all fundi
#==================
def extract_fundi(vertices, faces, depths, mean_curvatures, min_directions,
                  depth_threshold=0.2, thr=0.5, min_sulcus_size=50,
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
    depth_threshold: depth threshold for defining sulci
    thr: likelihood threshold
    min_sulcus_size: minimum sulcus size from which to find a fundus
    min_distance: minimum distance
    max_distance: maximum distance

    Output:
    ------
    fundi: list of #sulci lists of vertex indices.

    Calls:
    -----
    extract_sulci()
    compute_likelihood()
    find_anchors()
    connect_anchors()

    """

    load_em = 0
    save_em = 1
    if load_em:
        import pickle
        load_path = "/drop/yrjo_code_io_data/"

    # Extract sulci (vertex indices for each sulcus)
    if load_em:
        sulci = pickle.load(open(load_path + "sulci.p","rb"))
        n_sulci = int(pickle.load(open(load_path + "n_sulci.p","rb")))
    else:
        t0 = time()
        sulci, n_sulci = extract_sulci(faces, depths, depth_threshold,
                                       min_sulcus_size)
        print(str(time() - t0) + 'seconds')
        if save_em:
            pickle.dump(sulci, open(load_path + "sulci.p","wb"))
            pickle.dump(n_sulci, open(load_path + "n_sulci.p","wb"))

    # For each sulcus...
    print('Extract fundi from ' + str(n_sulci) + ' sulci...')
    fundi = []
    for i_sulcus, sulcus in enumerate(sulci):

        # Compute fundus likelihood values
        print('Compute fundus likelihood for sulcus ' +
              str(i_sulcus + 1) + '...')
        t0 = time()
        sulcus_likelihoods = compute_likelihood(depths[sulcus],
                                                mean_curvatures[sulcus])
        print(str(time() - t0) + 'seconds: likelihood')
        if save_em:
            pickle.dump(sulcus_likelihoods, open(load_path + "sulcus_likelihoods"+str(i_sulcus)+".p","wb"))

        # If the sulcus has enough high-likelihood vertices, continue
        if sum(sulcus_likelihoods > thr) > min_sulcus_size:

            # Find fundus points
            print('Find fundus points for sulcus ' +
                  str(i_sulcus + 1) + '...')
            t0 = time()
            anchors = find_anchors(vertices[sulcus, :], sulcus_likelihoods,
                                   min_directions[sulcus],
                                   thr, min_distance, max_distance)
            print(str(time() - t0) + 'seconds: anchors')
            if len(anchors) > 0:

                if save_em:
                    pickle.dump(anchors, open(load_path + "anchors"+str(i_sulcus)+".p","wb"))

                # Remove faces that have a non-sulcus vertex (output: 1-D array)
                fs = frozenset(sulcus)
                faces_sulcus1 = np.ravel([lst for lst in faces
                                          if len(fs.intersection(lst))==3])
                if save_em:
                    pickle.dump(faces_sulcus1, open(load_path + "faces_sulcus1_"+str(i_sulcus)+".p","wb"))

                # Replace mesh indices with sulcus indices
                print('Replace mesh indices with sulcus indices...')
                t0 = time()
                # anchor vertices
                V = np.unique(faces_sulcus1)
                anchors = np.array(anchors)
                anchors2 = anchors.copy()
                for index_new, index_old in enumerate(V):
                    anchors2[anchors == index_old] = index_new
                # faces vertices
                faces_table = np.transpose([[i,x] for i,x in enumerate(V)])
                faces_sulcus2 = [faces_table[0, faces_table[1, :] == x]
                                 for x in faces_sulcus1]
                faces_sulcus2 = np.reshape(np.ravel(faces_sulcus2), (-1, 3))
                print(str(time() - t0) + 'seconds: reindexing')
                if save_em:
                    pickle.dump(faces_sulcus2, open(load_path + "faces_sulcus2_"+str(i_sulcus)+".p","wb"))

                # Connect fundus points and extract fundus
                print('Connect fundus points for sulcus ' +
                      str(i_sulcus + 1) + '...')
                t0 = time()
                fundi.append(
                      connect_anchors(anchors2, faces_sulcus2,
                                      sulcus_likelihoods, thr))
                print(str(time() - t0) + 'seconds: connect')
            else:
                fundi.append([])
        else:
            fundi.append([])

    if save_em:
        pickle.dump(fundi, open(load_path + "fundi.p","wb"))
    return fundi

load_em = 1
if load_em:
    import pickle
    load_path = "/drop/yrjo_code_io_data/"
    min_directions = pickle.load(open(load_path + "min_directions.p","rb"))
    mean_curvatures = pickle.load(open(load_path + "mean_curvatures.p","rb"))
    depths = pickle.load(open(load_path + "depths.p","rb"))
    vertices = pickle.load(open(load_path + "vertices.p","rb"))
    faces = pickle.load(open(load_path + "faces.p","rb"))
else:
    from test import test_fundi_hmmf
    mean_curvatures, depths, vertices, faces, min_directions = test_fundi_hmmf()
    #mean_curvatures, depths, vertices, faces, min_directions, output_sulci, \
    #output_anchor_points, output_L, output_fundi = test_fundi_hmmf()

#sulci, n_sulci = extract_sulci(faces, depths, depth_threshold=0.2, min_sulcus_size=50)

fundi = extract_fundi(vertices, faces, depths, mean_curvatures, min_directions,
    depth_threshold=0.2, thr=0.5, min_sulcus_size=50,
    min_distance=5, max_distance=8)

# import time
# t0 = time.time()
# print time.time() - t0, "seconds"
