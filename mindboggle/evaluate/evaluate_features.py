#!/usr/bin/env python
"""
Evaluate deep surface features by computing the minimum distance from each
label border vertex to all of the feature vertices in the same sulcus.
The label boundaries run along the deepest parts of many sulci and
correspond to fundi in the DKT cortical labeling protocol.

Examples
--------
$ python evaluate_features.py


Authors:
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)
    - Arno Klein, 2012-2014  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2014,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def evaluate_deep_features(features_file, labels_file, sulci_file='', hemi='',
                           excludeIDs=[-1], output_vtk_name='', verbose=True):
    """
    Evaluate deep surface features by computing the minimum distance from each
    label border vertex to all of the feature vertices in the same sulcus,
    and from each feature vertex to all of the label border vertices in the
    same sulcus.  The label boundaries run along the deepest parts of sulci
    and correspond to fundi in the DKT cortical labeling protocol.

    Parameters
    ----------
    features_file : string
        VTK surface file with feature numbers for vertex scalars
    labels_file : string
        VTK surface file with label numbers for vertex scalars
    sulci_file : string
        VTK surface file with sulcus numbers for vertex scalars
    excludeIDs : list of integers
        feature/sulcus/label IDs to exclude (background set to -1)
    output_vtk_name : Boolean
        if not empty, output a VTK file beginning with output_vtk_name that
        contains a surface with mean distances as scalars
    verbose : Boolean
        print mean distances to standard output?

    Returns
    -------
    feature_to_border_mean_distances : numpy array [number of features x 1]
        mean distance from each feature to sulcus label border
    feature_to_border_sd_distances : numpy array [number of features x 1]
        standard deviations of feature-to-border distances
    feature_to_border_distances_vtk : string
        VTK surface file containing feature-to-border distances
    border_to_feature_mean_distances : numpy array [number of features x 1]
        mean distances from each sulcus label border to feature
    border_to_feature_sd_distances : numpy array [number of features x 1]
        standard deviations of border-to-feature distances
    border_to_feature_distances_vtk : string
        VTK surface file containing border-to-feature distances

    """
    import os
    import sys
    import numpy as np
    from mindboggle.utils.io_vtk import read_vtk, read_scalars, write_vtk
    from mindboggle.utils.mesh import find_neighbors, remove_faces
    from mindboggle.utils.segment import extract_borders
    from mindboggle.utils.compute import source_to_target_distances
    from mindboggle.LABELS import DKTprotocol

    dkt = DKTprotocol()
    #-------------------------------------------------------------------------
    # Load labels, features, and sulci:
    #-------------------------------------------------------------------------
    faces, lines, indices, points, npoints, labels, scalar_names, \
        input_vtk = read_vtk(labels_file, True, True)
    features, name = read_scalars(features_file, True, True)
    if sulci_file:
        sulci, name = read_scalars(sulci_file, True, True)
        # List of indices to sulcus vertices:
        sulcus_indices = [i for i,x in enumerate(sulci) if x != -1]
        segmentIDs = sulci
        sulcus_faces = remove_faces(faces, sulcus_indices)
    else:
        sulcus_indices = range(len(labels))
        segmentIDs = []
        sulcus_faces = faces

    #-------------------------------------------------------------------------
    # Prepare neighbors, label pairs, border IDs, and outputs:
    #-------------------------------------------------------------------------
    # Calculate neighbor lists for all points:
    print('Find neighbors to all vertices...')
    neighbor_lists = find_neighbors(faces, npoints)

    # Find label border points in any of the sulci:
    print('Find label border points in any of the sulci...')
    border_indices, border_label_tuples, unique_border_label_tuples = \
        extract_borders(sulcus_indices, labels, neighbor_lists,
                        ignore_values=[], return_label_pairs=True)
    if not len(border_indices):
        sys.exit('There are no label border points!')

    # Initialize an array of label border IDs
    # (label border vertices that define sulci in the labeling protocol):
    print('Build an array of label border IDs...')
    label_borders = -1 * np.ones(npoints)

    if hemi == 'lh':
        nsulcus_lists = len(dkt.left_sulcus_label_pair_lists)
    else:
        nsulcus_lists = len(dkt.right_sulcus_label_pair_lists)
    feature_to_border_mean_distances = -1 * np.ones(nsulcus_lists)
    feature_to_border_sd_distances = -1 * np.ones(nsulcus_lists)
    border_to_feature_mean_distances = -1 * np.ones(nsulcus_lists)
    border_to_feature_sd_distances = -1 * np.ones(nsulcus_lists)
    feature_to_border_distances_vtk = ''
    border_to_feature_distances_vtk = ''

    #-------------------------------------------------------------------------
    # Loop through sulci:
    #-------------------------------------------------------------------------
    # For each list of sorted label pairs (corresponding to a sulcus):
    for isulcus, label_pairs in enumerate(dkt.sulcus_label_pair_lists):

        # Keep the border points with label pair labels:
        border_indices = [x for i,x in enumerate(border_indices)
                          if np.unique(border_label_tuples[i]).tolist()
                          in label_pairs]

        # Store the points as sulcus IDs in the border IDs array:
        if border_indices:
            label_borders[border_indices] = isulcus

    if len(np.unique(label_borders)) > 1:

        #---------------------------------------------------------------------
        # Construct a feature-to-border distance matrix and VTK file:
        #---------------------------------------------------------------------
        # Construct a distance matrix:
        print('Construct a feature-to-border distance matrix...')
        sourceIDs = features
        targetIDs = label_borders
        distances, distance_matrix = source_to_target_distances(
            sourceIDs, targetIDs, points, segmentIDs, excludeIDs)

        # Compute mean distances for each feature:
        nfeatures = min(np.shape(distance_matrix)[1], nsulcus_lists)
        for ifeature in range(nfeatures):
            feature_distances = [x for x in distance_matrix[:, ifeature]
                                 if x != -1]
            feature_to_border_mean_distances[ifeature] = \
                np.mean(feature_distances)
            feature_to_border_sd_distances[ifeature] = \
                np.std(feature_distances)

        if verbose:
            print('Feature-to-border mean distances:')
            print(feature_to_border_mean_distances)
            print('Feature-to-border standard deviations of distances:')
            print(feature_to_border_sd_distances)

        # Write resulting feature-label border distances to VTK file:
        if output_vtk_name:
            feature_to_border_distances_vtk = os.path.join(os.getcwd(),
                output_vtk_name + '_feature_to_border_mean_distances.vtk')
            print('Write feature-to-border distances to {0}...'.
                  format(feature_to_border_distances_vtk))
            write_vtk(feature_to_border_distances_vtk, points,
                      [], [], sulcus_faces, [distances],
                      ['feature-to-border_distances'], 'float')

        #---------------------------------------------------------------------
        # Construct a border-to-feature distance matrix and VTK file:
        #---------------------------------------------------------------------
        # Construct a distance matrix:
        print('Construct a border-to-feature distance matrix...')
        sourceIDs = label_borders
        targetIDs = features
        distances, distance_matrix = source_to_target_distances(
            sourceIDs, targetIDs, points, segmentIDs, excludeIDs)

        # Compute mean distances for each feature:
        nfeatures = min(np.shape(distance_matrix)[1], nsulcus_lists)
        for ifeature in range(nfeatures):
            border_distances = [x for x in distance_matrix[:, ifeature]
                                if x != -1]
            border_to_feature_mean_distances[ifeature] = \
                np.mean(border_distances)
            border_to_feature_sd_distances[ifeature] = \
                np.std(border_distances)

        if verbose:
            print('border-to-feature mean distances:')
            print(border_to_feature_mean_distances)
            print('border-to-feature standard deviations of distances:')
            print(border_to_feature_sd_distances)

        # Write resulting feature-label border distances to VTK file:
        if output_vtk_name:
            border_to_feature_distances_vtk = os.path.join(os.getcwd(),
                output_vtk_name + '_border_to_feature_mean_distances.vtk')
            print('Write border-to-feature distances to {0}...'.
                  format(border_to_feature_distances_vtk))
            write_vtk(border_to_feature_distances_vtk, points,
                      [], [], sulcus_faces, [distances],
                      ['border-to-feature_distances'], 'float')

    #-------------------------------------------------------------------------
    # Return outputs:
    #-------------------------------------------------------------------------
    return feature_to_border_mean_distances, feature_to_border_sd_distances,\
           feature_to_border_distances_vtk,\
           border_to_feature_mean_distances, border_to_feature_sd_distances,\
           border_to_feature_distances_vtk



if __name__ == "__main__":

    import os
    import numpy as np

    from mindboggle.evaluate.evaluate_features import evaluate_deep_features

    feature_type = 'fundi'  # Otherwise measure distances to all vertices

    mindboggled = '/homedir/mindboggled'
    output_dir = '/desk'

    if feature_type == 'fundi':
        feature_dir = '/homedir/Projects/BrainImaging/fundus_evaluation_2014'
    else:
        feature_dir = ''

    # Features: fundus method:
    # 0 = mindboggle
    # 1 = forrest scalars
    # 2 = forrest lines
    # 3 = gang li
    # 4 = olivier coulon
    nmethod = 0
    fmethods = [mindboggled,
                os.path.join(feature_dir, 'Forrest_Fundi.scalars'),
                os.path.join(feature_dir, 'Forrest_Fundi.lines'),
                os.path.join(feature_dir, 'GangLiFundi'),
                os.path.join(feature_dir, 'Olivier_fundi')]
    fmethod = fmethods[nmethod]

    # Loop through subjects and hemispheres:
    names = ['OASIS-TRT-20'] #, 'MMRR-21', 'NKI-RS-22', 'NKI-TRT-20']
            #['Afterthought', 'Colin27', 'Twins-2', 'MMRR-3T7T-2', 'HLN-12']
    numbers = [1] #[20,21,22,20]  #[1,1,2,2,12]
    surfs = ['left_surface','right_surface']
    hemis = ['lh','rh']
    nsubjects = sum(numbers)
    feature_to_border_mean_distances_left = -1 * np.ones((nsubjects, 25))
    feature_to_border_sd_distances_left = -1 * np.ones((nsubjects, 25))
    border_to_feature_mean_distances_left = -1 * np.ones((nsubjects, 25))
    border_to_feature_sd_distances_left = -1 * np.ones((nsubjects, 25))
    feature_to_border_mean_distances_right = -1 * np.ones((nsubjects, 25))
    feature_to_border_sd_distances_right = -1 * np.ones((nsubjects, 25))
    border_to_feature_mean_distances_right = -1 * np.ones((nsubjects, 25))
    border_to_feature_sd_distances_right = -1 * np.ones((nsubjects, 25))
    isubject = 0
    for iname, name in enumerate(names):
        number = numbers[iname]
        for n in range(1,number+1):
            subject = name+'-'+str(n)
            for isurf, surf in enumerate(surfs):
                hemi = hemis[isurf]
                print('{0}: {1}'.format(subject, hemi))

                # Identify surface files with labels and sulci:
                mdir = os.path.join(fmethod, subject)
                sulci_file = os.path.join(mdir, 'features', surf, 'sulci.vtk')



                #labels_file = os.path.join(mdir, 'labels', surf,
                #                           'relabeled_labels.DKT31.manual.vtk')
                labels_file = os.path.join(mdir, 'labels', surf,
                                           'relabeled_classifier.vtk')


                # Identify features file:
                if nmethod == 0:
                    features_file = os.path.join(mdir, 'features', surf,
                                                 'border_per_sulcus.vtk')
                else:
                    features_file = os.path.join(fmethod,
                        '_hemi_' + hemi + '_subject_' + name + '-' + str(n),
                        hemi + '.pial.fundi.vtk')

                # Compute distances between features and label boundaries
                # in sulci corresponding to fundi:
                feature_to_border_mean_distances, \
                feature_to_border_sd_distances,\
                feature_to_border_distances_vtk,\
                border_to_feature_mean_distances, \
                border_to_feature_sd_distances,\
                border_to_feature_distances_vtk = \
                    evaluate_deep_features(features_file, labels_file,
                                           sulci_file, hemi, excludeIDs=[-1],
                                           output_vtk_name=subject+'_'+hemi,
                                           verbose=True)
                print('*' * 79)

                if isurf == 0:
                    feature_to_border_mean_distances_left[isubject, :] = \
                        feature_to_border_mean_distances
                    feature_to_border_sd_distances_left[isubject, :] = \
                       feature_to_border_sd_distances
                    border_to_feature_mean_distances_left[isubject, :] = \
                        border_to_feature_mean_distances
                    border_to_feature_sd_distances_left[isubject, :] = \
                        border_to_feature_sd_distances
                else:
                    feature_to_border_mean_distances_right[isubject, :] = \
                        feature_to_border_mean_distances
                    feature_to_border_sd_distances_right[isubject, :] = \
                       feature_to_border_sd_distances
                    border_to_feature_mean_distances_right[isubject, :] = \
                        border_to_feature_mean_distances
                    border_to_feature_sd_distances_right[isubject, :] = \
                        border_to_feature_sd_distances

            isubject += 1

    # Save tables of mean distances:
    np.savetxt('feature_to_border_mean_distances_left.csv',
               feature_to_border_mean_distances_left)
    np.savetxt('feature_to_border_sd_distances_left.csv',
               feature_to_border_sd_distances_left)
    np.savetxt('border_to_feature_mean_distances_left.csv',
               border_to_feature_mean_distances_left)
    np.savetxt('border_to_feature_sd_distances_left.csv',
               border_to_feature_sd_distances_left)

    np.savetxt('feature_to_border_mean_distances_right.csv',
               feature_to_border_mean_distances_right)
    np.savetxt('feature_to_border_sd_distances_right.csv',
               feature_to_border_sd_distances_right)
    np.savetxt('border_to_feature_mean_distances_right.csv',
               border_to_feature_mean_distances_right)
    np.savetxt('border_to_feature_sd_distances_right.csv',
               border_to_feature_sd_distances_right)

    # Average of the mean distances across all subjects:
    mean_feature_to_border_mean_distances_left = \
        np.mean(feature_to_border_mean_distances_left, axis=0)
    mean_feature_to_border_sd_distances_left = \
        np.mean(feature_to_border_mean_distances_left, axis=0)
    mean_border_to_feature_mean_distances_left = \
        np.mean(border_to_feature_mean_distances_left, axis=0)
    mean_border_to_feature_sd_distances_left = \
        np.mean(border_to_feature_mean_distances_left, axis=0)

    mean_feature_to_border_mean_distances_right = \
        np.mean(feature_to_border_mean_distances_right, axis=0)
    mean_feature_to_border_sd_distances_right = \
        np.mean(feature_to_border_mean_distances_right, axis=0)
    mean_border_to_feature_mean_distances_right = \
        np.mean(border_to_feature_mean_distances_right, axis=0)
    mean_border_to_feature_sd_distances_right = \
        np.mean(border_to_feature_mean_distances_right, axis=0)

    # Save tables of mean distances averaged across all subjects:
    np.savetxt('mean_feature_to_border_mean_distances_left.csv',
               mean_feature_to_border_mean_distances_left)
    np.savetxt('mean_feature_to_border_sd_distances_left.csv',
               mean_feature_to_border_sd_distances_left)
    np.savetxt('mean_border_to_feature_mean_distances_left.csv',
               mean_border_to_feature_mean_distances_left)
    np.savetxt('mean_border_to_feature_sd_distances_left.csv',
               mean_border_to_feature_sd_distances_left)

    np.savetxt('mean_feature_to_border_mean_distances_right.csv',
               mean_feature_to_border_mean_distances_right)
    np.savetxt('mean_feature_to_border_sd_distances_right.csv',
               mean_feature_to_border_sd_distances_right)
    np.savetxt('mean_border_to_feature_mean_distances_right.csv',
               mean_border_to_feature_mean_distances_right)
    np.savetxt('mean_border_to_feature_sd_distances_right.csv',
               mean_border_to_feature_sd_distances_right)
