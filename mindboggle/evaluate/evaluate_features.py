#!/usr/bin/env python
"""
Evaluate deep surface features by computing the minimum distance from each
label boundary vertex to all of the feature vertices in the same sulcus.
The label boundaries run along the deepest parts of many sulci and
correspond to fundi in the DKT cortical labeling protocol.

Examples
--------
$ python evaluate_fundi.py


Authors:
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)
    - Arno Klein, 2012-2014  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2014,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def evaluate_deep_features(features_file, labels_file, sulci_file='',
                           output_vtk='', excludeIDs=[-1]):
    """
    Evaluate deep surface features by computing the minimum distance from each
    label boundary vertex to all of the feature vertices in the same sulcus.
    The label boundaries run along the deepest parts of sulci and
    correspond to fundi in the DKT cortical labeling protocol.

    Parameters
    ----------
    features_file : string
        VTK surface file with feature numbers for vertex scalars
    labels_file : string
        VTK surface file with label numbers for vertex scalars
    sulci_file : string
        VTK surface file with sulcus numbers for vertex scalars
    output_vtk : string
        name of output VTK file containing surface with distances as scalars
        (not written if empty string)
    excludeIDs : list of integers
        feature/sulcus/label IDs to exclude (background set to -1)

    Returns
    -------
    feature_mean_distances: numpy array [number of features x 1]
        mean distances from source to target features
    feature_std_distances: numpy array [number of features x 1]
        standard deviations of distances from source to target features

    """
    import sys
    import numpy as np
    from mindboggle.utils.io_vtk import read_vtk, read_scalars, write_vtk
    from mindboggle.utils.mesh import find_neighbors
    from mindboggle.utils.segment import extract_borders
    from mindboggle.utils.compute import source_to_target_distances
    from mindboggle.LABELS import DKTprotocol

    dkt = DKTprotocol()

    # Load labels, features, and sulci:
    faces, lines, indices, points, npoints, labels, scalar_names, \
        input_vtk = read_vtk(labels_file, True, True)
    features, name = read_scalars(features_file, True, True)
    if sulci_file:
        sulci, name = read_scalars(sulci_file, True, True)
        # List of indices to sulcus vertices:
        sulcus_indices = [i for i,x in enumerate(sulci) if x != -1]
        segmentIDs = sulci
    else:
        sulcus_indices = range(len(labels))
        segmentIDs = []

    # Calculate neighbor lists for all points:
    print('Find neighbors to all vertices...')
    neighbor_lists = find_neighbors(faces, npoints)

    # Prepare list of all unique sorted label pairs in the labeling protocol:
    print('Prepare a list of unique, sorted label pairs in the protocol...')
    nsulcus_lists = len(dkt.sulcus_label_pair_lists)

    # Find label boundary points in any of the sulci:
    print('Find label boundary points in any of the sulci...')
    border_indices, border_label_tuples, unique_border_label_tuples = \
        extract_borders(sulcus_indices, labels, neighbor_lists,
                        ignore_values=[], return_label_pairs=True)
    if not len(border_indices):
        sys.exit('There are no label boundary points!')

    # Initialize an array of label boundaries fundus IDs
    # (label boundary vertices that define sulci in the labeling protocol):
    print('Build an array of label boundary fundus IDs...')
    label_boundary_fundi = -1 * np.ones(npoints)

    # For each list of sorted label pairs (corresponding to a sulcus):
    for isulcus, label_pairs in enumerate(dkt.sulcus_label_pair_lists):
        #print('  Sulcus ' + str(isulcus + 1))

        # Keep the boundary points with label pair labels:
        fundus_indices = [x for i,x in enumerate(border_indices)
                          if np.unique(border_label_tuples[i]).tolist()
                              in label_pairs]

        # Store the points as sulcus IDs in the fundus IDs array:
        if fundus_indices:
            label_boundary_fundi[fundus_indices] = isulcus + 1

    feature_mean_distances = -1 * np.ones(nsulcus_lists)
    feature_std_distances = -1 * np.ones(nsulcus_lists)

    if len(np.unique(label_boundary_fundi)) > 1:

        # Construct a distance matrix:
        print('Construct a boundary-to-feature distance matrix...')
        sourceIDs = features
        targetIDs = label_boundary_fundi
        ntargets = nsulcus_lists
        distances, distance_matrix = source_to_target_distances(
            sourceIDs, targetIDs, points, ntargets, segmentIDs,
            excludeIDs)

        # Compute mean distances for each feature:
        for ifeature in range(nsulcus_lists):
            feature_distances = [x for x in distance_matrix[:, ifeature]
                                 if x != -1]
            feature_mean_distances[ifeature] = np.mean(feature_distances)
            feature_std_distances[ifeature] = np.std(feature_distances)

        print('Feature mean distances')
        print(feature_mean_distances)
        print('Feature standard deviations of distances')
        print(feature_std_distances)

        # Write resulting feature-label boundary distances to VTK file:
        if output_vtk:
            print('Write distances to a VTK file for visualization...')
            write_vtk(output_vtk, points, [], [], faces,
                      [distances], ['feature-label_boundary_distances'],
                      'float')

    return feature_mean_distances, feature_std_distances, output_vtk


if __name__ == "__main__":

    import os
    from mindboggle.evaluate.evaluate_features import evaluate_deep_features

    mindboggled = '/homedir/mindboggled'
    fundi_dir = '/home/arno/Projects/BrainImaging/fundus_evaluation_2014'
    output_dir = '/desk'

    # Features: fundus method:
    # 0 = mindboggle
    # 1 = forrest scalars
    # 2 = forrest lines
    # 3 = gang li
    # 4 = olivier coulon
    nmethod = 0
    fmethods = [mindboggled,
                os.path.join(fundi_dir, 'Forrest_Fundi.scalars'),
                os.path.join(fundi_dir, 'Forrest_Fundi.lines'),
                os.path.join(fundi_dir, 'GangLiFundi'),
                os.path.join(fundi_dir, 'Olivier_fundi')]
    fmethod = fmethods[nmethod]

    # Loop through subjects and hemispheres:
    names = ['OASIS-TRT-20', 'MMRR-21', 'NKI-RS-22', 'NKI-TRT-20']
            #['Afterthought', 'Colin27', 'Twins-2', 'MMRR-3T7T-2', 'HLN-12']
    numbers = [20,21,22,20]  #[1,1,2,2,12]
    surfs = ['left_surface','right_surface']
    hemis = ['lh','rh']
    for iname, name in enumerate(names):
        number = numbers[iname]
        for n in range(1,number+1):
            subject = name+'-'+str(n)
            for isurf, surf in enumerate(surfs):
                hemi = hemis[isurf]

                # Identify surface files with labels and sulci:
                mdir = os.path.join(fmethod, subject)
                labels_file = os.path.join(mdir, 'labels', surf,
                                           'relabeled_labels.DKT31.manual.vtk')
                sulci_file = os.path.join(mdir, 'features', surf, 'sulci.vtk')
                output_vtk = os.path.join(output_dir,
                                          subject + '_fundus-label_distances.vtk')

                # Identify features file:
                if nmethod == 0:
                    features_file = os.path.join(mdir, 'features', surf,
                                                 'fundus_per_sulcus.vtk')
                else:
                    features_file = os.path.join(fmethod,
                        '_hemi_' + hemi + '_subject_' + name + '-' + str(n),
                        hemi + '.pial.fundi.vtk')

                # Compute distances from features to label boundaries in sulci
                # corresponding to fundi:
                feature_mean_distances, feature_std_distances, \
                output_vtk = evaluate_deep_features(features_file, labels_file,
                                                    sulci_file, output_vtk,
                                                    excludeIDs=[-1])
                print('***')
                print('Input features:' + features_file)
                print('Input sulci:' + sulci_file)
                print('Input labels:' + labels_file)
                print('Output VTK file:' + output_vtk)
                print('feature_mean_distances:')
                print(feature_mean_distances)
                print('***')
