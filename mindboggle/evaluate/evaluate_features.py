#!/usr/bin/env python
"""
Evaluate deep surface features by computing the minimum distance from each
label boundary vertex to all of the feature vertices in the same fold.
The label boundaries run along the deepest parts of many folds and
correspond to fundi in the DKT cortical labeling protocol.

Examples
--------
$ python evaluate_fundi.py fundi.vtk folds.vtk lh.labels.DKT25.manual.vtk


Authors:
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)
    - Arno Klein, 2012-2014  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2014,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def evaluate_deep_features(features_file, labels_file, folds_file=[],
                           output_vtk=''):
    """
    Evaluate deep surface features by computing the minimum distance from each
    label boundary vertex to all of the feature vertices in the same fold.
    The label boundaries run along the deepest parts of many folds and
    correspond to fundi in the DKT cortical labeling protocol.

    Parameters
    ----------
    features_file : string
        VTK surface file with feature numbers for vertex scalars
    folds_file : string
        VTK surface file with fold numbers for vertex scalars
    labels_file : string
        VTK surface file with label numbers for vertex scalars
    output_vtk : string
        name of output VTK file containing surface with distances as scalars

    Returns
    -------
    feature_mean_distances: numpy array [number of features x 1]
        mean distances from source to target features
    feature_std_distances: numpy array [number of features x 1]
        standard deviations of distances from source to target features

    """
    import sys
    import numpy as np
    from mindboggle.utils.io_vtk import load_vtk, read_scalars, write_scalars
    from mindboggle.utils.mesh import find_neighbors
    from mindboggle.utils.segment import extract_borders
    from mindboggle.utils.compute import source_to_target_distances
    from mindboggle.LABELS import DKTprotocol

    dkt = DKTprotocol()

    # Load features, folds, labels
    features, name = read_scalars(features_file, return_arrays=True)
    folds, name = read_scalars(folds_file, return_arrays=True)
    faces, lines, indices, points, npoints, labels, \
        scalar_names = load_vtk(labels_file, return_arrays=True)

    # List of indices to fold vertices
    fold_indices = [i for i,x in enumerate(folds) if x > 0]

    # Calculate neighbor lists for all points
    print('Find neighbors to all vertices...')
    neighbor_lists = find_neighbors(faces, npoints)

    # Prepare list of all unique sorted label pairs in the labeling protocol
    print('Prepare a list of unique, sorted label pairs in the protocol...')
    n_features = len(dkt.sulcus_label_pair_lists)

    # Find label boundary points in any of the folds
    print('Find label boundary points in any of the folds...')
    border_indices, border_label_tuples, unique_border_label_tuples = \
        extract_borders(fold_indices, labels, neighbor_lists)
    if not len(border_indices):
        sys.exit('There are no label boundary points!')

    # Initialize an array of label boundaries fundus IDs
    # (label boundary vertices that define sulci in the labeling protocol)
    print('Build an array of label boundary fundus IDs...')
    label_boundary_fundi = np.zeros(npoints)

    # For each list of sorted label pairs (corresponding to a sulcus)
    for isulcus, label_pairs in enumerate(dkt.sulcus_label_pair_lists):
        print('  Sulcus ' + str(isulcus + 1))

        # Keep the boundary points with label pair labels
        fundus_indices = [x for i,x in enumerate(border_indices)
                          if np.unique(border_label_tuples[i]).tolist()
                              in label_pairs]

        # Store the points as sulcus IDs in the fundus IDs array
        label_boundary_fundi[fundus_indices] = isulcus + 1

    # Construct a distance matrix
    print('Construct a boundary-to-feature distance matrix...')
    distances, distance_matrix = source_to_target_distances(
        sourceIDs=label_boundary_fundi,
        targetIDs=features, points=points, segmentIDs=[])

    # Compute mean distances for each feature
    feature_mean_distances = np.zeros(n_features)
    feature_std_distances = np.zeros(n_features)
    for ifeature in range(n_features):
        feature_distances = distance_matrix[:, ifeature]
        feature_distances = [x for x in feature_distances if x != -1]
        feature_mean_distances[ifeature] = np.mean(feature_distances)
        feature_std_distances[ifeature] = np.std(feature_distances)

    print('Feature mean distances')
    print(feature_mean_distances)
    print('Feature standard deviations of distances')
    print(feature_std_distances)

    # Write resulting feature-label boundary distances to VTK file:
    if output_vtk:
        print('Write distances to a VTK file for visualization...')
        write_scalars(output_vtk, points, range(npoints), faces,
                      [distances], ['feature-label_boundary_distances'])

    return feature_mean_distances, feature_std_distances, output_vtk


if __name__ == "__main__":

    import os
    from mindboggle.evaluate import evaluate_deep_features

    subjects_dir = '/data/Brains/Mindboggle101/subjects'
    mindboggle_dir = '/homedir/mindboggled'
    output_dir = '/desk/fundi/'

    names = ['OASIS-TRT-20', 'MMRR-21', 'NKI-RS-22', 'NKI-TRT-20',
             'Afterthought', 'Colin27', 'Twins-2', 'MMRR-3T7T-2', 'HLN-12']
    numbers = [20,21,22,20,
               1,1,2,2,12]

    for i,name in enumerate(names):
        number = numbers[i]
        for n in range(1,number+1):
            subject = name+'-'+str(n)
            sdir = os.path.join(subjects_dir, subject, 'mri')
            mdir = os.path.join(mindboggle_dir, subject)

            labels_file = os.path.join(sdir, 'labels.nii.gz')
            features_file = os.path.join(mdir, 'features', 'fundi.vtk')
            folds_file = os.path.join(mdir, 'features', 'folds.vtk')
            output_vtk = os.path.join(output_dir,
                                      subject + '_fundus-label_distances.vtk')

            # Compute distances from features to label boundaries in folds
            # corresponding to fundi:
            feature_mean_distances, feature_std_distances, \
            output_vtk = evaluate_deep_features(features_file, labels_file,
                                                folds_file, output_vtk)

            print('***')
            print('Input features:' + features_file)
            print('Input folds:' + folds_file)
            print('Input labels:' + labels_file)
            print('Output VTK file:' + output_vtk)
            print('feature_mean_distances:')
            print(feature_mean_distances)
            print('***')
