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
from find_neighbors import find_neighbors
from find_anchors import find_anchors
from connect_anchors import connect_anchors
from test import test_fundi_hmmf

#==================
# Extract all fundi
#==================
def extract_fundi(vertices, faces, depths, mean_curvatures, min_directions,
                  depth_threshold=0.2, min_sulcus_size=50, thr=0.5,
                  min_distance=5, max_distance=8):
#max_neighbors=15
    """
    Extract all fundi.

    A fundus is a connected set of high-likelihood vertices in a surface mesh.

    Inputs:
    ------
    L: fundus likelihood values [#vertices x 1] numpy array
    sulci: sulcus values [#vertices x 1] numpy array
    sulcus_index: sulcus number [int]
    vertices: [#vertices x 3] numpy array
    faces: vertices for polygons [#faces x 3] numpy array
    min_distances: [#vertices x 1] numpy array
    thr: likelihood threshold
    min_sulcus_size: minimum sulcus size from which to find a fundus
    min_distance: minimum distance
    max_distance: maximum distance

#    max_neighbors: maximum number of neighbors per sulcus vertex

    vertices: [#vertices x 3] numpy array
    faces: vertices for polygons [#faces x 3] numpy array
    depths: depth values [#vertices x 1] numpy array
    mean_curvatures: mean curvature values [#vertices x 1] numpy array
    min_directions: directions of minimum curvature [3 x #vertices] numpy array
    depth_threshold: depth threshold for defining sulci
    thr: likelihood threshold

    Output:
    ------
    fundi: list of #sulci lists of vertex indices.

    Calls:
    -----
    extract_sulci()
    compute_likelihood()
    find_neighbors()
    find_anchors()
    connect_anchors()

    """

    load_em = 1
    if load_em:
        import pickle
        load_path = "/drop/yrjo_code_io_data/"

    # Extract sulci (vertex indices for each sulcus)
    if load_em:
        sulci = pickle.load(open(load_path + "sulci.p","rb"))
        n_sulci = int(pickle.load(open(load_path + "n_sulci.p","rb")))
    else:
        sulci, n_sulci = extract_sulci(faces, depths, depth_threshold,
                                       min_sulcus_size)

    # For each sulcus...
    print('Extract fundi from ' + str(n_sulci) + ' sulci...')
    fundi = []
    for i_sulcus, sulcus in enumerate(sulci):

        # Compute fundus likelihood values
        print('Compute fundus likelihood for sulcus ' + str(i_sulcus + 1) + '...')
        """
        if load_em:
            sulcus_likelihoods = pickle.load(open(load_path + "L.p","rb"))
        else:
        """
        sulcus_likelihoods = compute_likelihood(depths[sulcus],
                                                mean_curvatures[sulcus])

        # If the size of the sulcus (according to likelihood values)
        # is sufficiently large, continue
        if sum(sulcus_likelihoods > thr) > min_sulcus_size:
            """
            if load_em:

                L = pickle.load(open(load_path + "L.p","rb"))
                I = pickle.load(open(load_path + "I.p","rb"))
            else:
            """
            I = [i for i,x in enumerate(sulcus) if sulcus_likelihoods[i] > 0]
            L = sulcus_likelihoods[sulcus_likelihoods > 0]
                #np.array([x for x in sulcus_likelihoods if x > 0])

            # Find fundus points
            print('Find fundus points for sulcus ' + str(i_sulcus + 1) + '...')
            anchors = find_anchors(vertices[I], L, min_directions[I],
                                   thr, min_distance, max_distance)
            if len(anchors) > 0:

                # Remove faces that contain a non-sulcus vertex
                print('Remove faces that contain a non-sulcus vertex...')
                fs = frozenset(sulcus)
                faces_sulcus = [lst for lst in faces
                                if len(fs.intersection(lst))==3]
                faces_sulcus = np.reshape(np.ravel(faces_sulcus), (-1, 3))

                # Replace mesh indices with sulcus indices
                print('Replace mesh indices with sulcus indices')
                V = np.unique(faces_sulcus)
                anchors2 = anchors.copy()
                faces_sulcus2 = faces_sulcus.copy()
                for index_new, index_old in enumerate(V):
                    anchors2[anchors == index_old] = index_new
                    faces_sulcus2[faces_sulcus == index_old] = index_new
                print(faces_sulcus)
                print(faces_sulcus2)
                # Connect fundus points and extract fundus
                print('Connect fundus points for sulcus ' + str(i_sulcus + 1) + '...')
                fundi.append(
                      connect_anchors(anchors2, faces_sulcus2, L, thr))
            else:
                fundi.append([])
        else:
            fundi.append([])

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
    mean_curvatures, depths, vertices, faces, min_directions = test_fundi_hmmf()
    #mean_curvatures, depths, vertices, faces, min_directions, output_sulci, \
    #output_anchor_points, output_L, output_fundi = test_fundi_hmmf()

#sulci, n_sulci = extract_sulci(faces, depths, depth_threshold=0.2, min_sulcus_size=50)

fundi = extract_fundi(vertices, faces, depths, mean_curvatures, min_directions, depth_threshold=0.2)
