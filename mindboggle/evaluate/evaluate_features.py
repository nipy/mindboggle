#!/usr/bin/env python
"""
Evaluate deep surface features by computing the minimum distance from each
label border vertex to all of the feature vertices in the same sulcus.
The label borders run along the deepest parts of many sulci and
correspond to fundi in the DKT cortical labeling protocol.

Examples
--------
$ python evaluate_features.py


Authors:
    - Yrjo Hame, 2012  (yrjo.hame@gmail.com)
    - Arno Klein, 2012-2015  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2015,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def evaluate_deep_features(features_file, labels_file, sulci_file='', hemi='',
                           excludeIDs=[-1], output_vtk_name='', verbose=True):
    """
    Evaluate deep surface features by computing the minimum distance from each
    label border vertex to all of the feature vertices in the same sulcus,
    and from each feature vertex to all of the label border vertices in the
    same sulcus.  The label borders run along the deepest parts of sulci
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
    from mindboggle.io.vtks import read_vtk, read_scalars, write_vtk
    from mindboggle.guts.mesh import find_neighbors, remove_faces
    from mindboggle.guts.segment import extract_borders
    from mindboggle.guts.compute import source_to_target_distances
    from mindboggle.io.labels import DKTprotocol

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
    print('Find neighbors for all vertices...')
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
        label_pair_border_indices = [x for i,x in enumerate(border_indices)
                          if np.unique(border_label_tuples[i]).tolist()
                          in label_pairs]

        # Store the points as sulcus IDs in the border IDs array:
        if label_pair_border_indices:
            label_borders[label_pair_border_indices] = isulcus

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

    #-------------------------------------------------------------------------
    # Set feature type ('fundi' or '' for every sulcus vertex), subjects:
    #-------------------------------------------------------------------------
    feature_type = 'fundi' #'sulci'  # If 'fundi', select 'nmethod' below.
    names = ['OASIS-TRT-20', 'MMRR-21', 'NKI-RS-22', 'NKI-TRT-20',
             'Afterthought', 'Colin27', 'Twins-2', 'MMRR-3T7T-2', 'HLN-12']
    numbers = [20,21,22,20, 1,1,2,2,12]
    mindboggled = '/mnt/nfs-share/Mindboggle101/mindboggled/manual'
    labels_dir = '/mnt/nfs-share/Mindboggle101/mindboggled/manual' 

    #-------------------------------------------------------------------------
    # Feature-specific settings:
    #-------------------------------------------------------------------------
    if feature_type == 'fundi':
        # Features: fundus method:
        # 0 = mindboggle
        # 1 = forrest scalars
        # 2 = forrest lines
        # 3 = gang li
        # 4 = olivier coulon
        nmethod = 0
        feature_dir = '/homedir/fundus_evaluation_2014/fundi_vtk'
        fmethods = ['mindboggle_fundi',
                    'ForrestBao_scalar_fundi',
                    'ForrestBao_line_fundi',
                    'GangLi_fundi',
                    'OlivierCoulon_fundi']
        fmethod_dirs = [mindboggled,
                        os.path.join(feature_dir, fmethods[1]),
                        os.path.join(feature_dir, fmethods[2]),
                        os.path.join(feature_dir, fmethods[3]),
                        os.path.join(feature_dir, fmethods[4])]
        fmethod_dir = fmethod_dirs[nmethod]
        fmethod = fmethods[nmethod]
    else:
        fmethod = 'all'

    #-------------------------------------------------------------------------
    # Miscellaneous defaults:
    #-------------------------------------------------------------------------
    surfs = ['left_cortical_surface', 'right_cortical_surface']
    hemis = ['lh', 'rh']
    nsulci = 25

    #-------------------------------------------------------------------------
    # Loop through subjects and hemispheres:
    #-------------------------------------------------------------------------
    nsubjects = sum(numbers)
    feature_to_border_mean_distances_left = -1 * np.ones((nsubjects, nsulci))
    feature_to_border_sd_distances_left = -1 * np.ones((nsubjects, nsulci))
    border_to_feature_mean_distances_left = -1 * np.ones((nsubjects, nsulci))
    border_to_feature_sd_distances_left = -1 * np.ones((nsubjects, nsulci))
    feature_to_border_mean_distances_right = -1 * np.ones((nsubjects, nsulci))
    feature_to_border_sd_distances_right = -1 * np.ones((nsubjects, nsulci))
    border_to_feature_mean_distances_right = -1 * np.ones((nsubjects, nsulci))
    border_to_feature_sd_distances_right = -1 * np.ones((nsubjects, nsulci))
    isubject = 0
    for iname, name in enumerate(names):
        number = numbers[iname]
        for n in range(1, number+1):
            subject = name+'-'+str(n)
            for isurf, surf in enumerate(surfs):
                hemi = hemis[isurf]
                #print('{0}: {1}'.format(subject, hemi))
                #-------------------------------------------------------------
                # Identify surface files with labels and with sulci:
                #-------------------------------------------------------------
                mdir = os.path.join(mindboggled, subject)
                ldir = os.path.join(labels_dir, subject)
                sulci_file = os.path.join(mdir, 'features', surf, 'sulci.vtk')
                labels_file = os.path.join(ldir, 'labels', surf,
                                           'relabeled_labels.DKT31.manual.vtk')
                #-------------------------------------------------------------
                # Identify features file:
                #-------------------------------------------------------------
                if feature_type == 'fundi':
                    if nmethod == 0:
                        features_file = os.path.join(mdir, 'features', surf,
                                                     'fundus_per_sulcus.vtk')
                    else:
                        features_file = os.path.join(fmethod_dir,
                            '_hemi_' + hemi + '_subject_' + name + '-' + str(n),
                            hemi + '.pial.fundi.vtk')
                else:
                    features_file = sulci_file
                #if not os.path.exists(features_file):
                #    print(features_file)

                #-------------------------------------------------------------
                # Compute distances between features and label borders
                # in sulci corresponding to fundi:
                #-------------------------------------------------------------
                if os.path.exists(features_file) \
                   and os.path.exists(labels_file) \
                   and os.path.exists(sulci_file):
                    feature_to_border_mean_distances, \
                    feature_to_border_sd_distances,\
                    feature_to_border_distances_vtk,\
                    border_to_feature_mean_distances, \
                    border_to_feature_sd_distances,\
                    border_to_feature_distances_vtk = \
                        evaluate_deep_features(features_file, labels_file,
                            sulci_file, hemi, excludeIDs=[-1],
                            output_vtk_name=subject+'_'+hemi+'_'+fmethod,
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

    #-------------------------------------------------------------------------
    # Save tables of mean distances:
    #-------------------------------------------------------------------------
    np.savetxt(fmethod + '_mean_distances_to_border_left.csv',
               feature_to_border_mean_distances_left)
    np.savetxt(fmethod + '_sd_distances_to_border_left.csv',
               feature_to_border_sd_distances_left)
    np.savetxt(fmethod + '_mean_distances_from_border_left.csv',
               border_to_feature_mean_distances_left)
    np.savetxt(fmethod + '_sd_distances_from_border_left.csv',
               border_to_feature_sd_distances_left)

    np.savetxt(fmethod + '_mean_distances_to_border_right.csv',
               feature_to_border_mean_distances_right)
    np.savetxt(fmethod + '_sd_distances_to_border_right.csv',
               feature_to_border_sd_distances_right)
    np.savetxt(fmethod + '_mean_distances_from_border_right.csv',
               border_to_feature_mean_distances_right)
    np.savetxt(fmethod + '_sd_distances_from_border_right.csv',
               border_to_feature_sd_distances_right)

    # #-------------------------------------------------------------------------
    # # Save tables of mean distances averaged across all subjects:
    # # NOTE: np.mean() results in nan's if any element has a nan.
    # #-------------------------------------------------------------------------
    # mean_feature_to_border_mean_distances_left = \
    #     np.mean(feature_to_border_mean_distances_left, axis=0)
    # mean_feature_to_border_sd_distances_left = \
    #     np.mean(feature_to_border_mean_distances_left, axis=0)
    # mean_border_to_feature_mean_distances_left = \
    #     np.mean(border_to_feature_mean_distances_left, axis=0)
    # mean_border_to_feature_sd_distances_left = \
    #     np.mean(border_to_feature_mean_distances_left, axis=0)
    #
    # mean_feature_to_border_mean_distances_right = \
    #     np.mean(feature_to_border_mean_distances_right, axis=0)
    # mean_feature_to_border_sd_distances_right = \
    #     np.mean(feature_to_border_mean_distances_right, axis=0)
    # mean_border_to_feature_mean_distances_right = \
    #     np.mean(border_to_feature_mean_distances_right, axis=0)
    # mean_border_to_feature_sd_distances_right = \
    #     np.mean(border_to_feature_mean_distances_right, axis=0)
    #
    # np.savetxt('avg_per_fundus_' + fmethod + '_mean_distances_to_border_left.csv',
    #            mean_feature_to_border_mean_distances_left)
    # np.savetxt('avg_per_fundus_' + fmethod + '_sd_distances_to_border_left.csv',
    #            mean_feature_to_border_sd_distances_left)
    # np.savetxt('avg_per_fundus_' + fmethod + '_mean_distances_from_border_left.csv',
    #            mean_border_to_feature_mean_distances_left)
    # np.savetxt('avg_per_fundus_' + fmethod + '_sd_distances_from_border_left.csv',
    #            mean_border_to_feature_sd_distances_left)
    #
    # np.savetxt('avg_per_fundus_' + fmethod + '_mean_distances_to_border_right.csv',
    #            mean_feature_to_border_mean_distances_right)
    # np.savetxt('avg_per_fundus_' + fmethod + '_sd_distances_to_border_right.csv',
    #            mean_feature_to_border_sd_distances_right)
    # np.savetxt('avg_per_fundus_' + fmethod + '_mean_distances_from_border_right.csv',
    #            mean_border_to_feature_mean_distances_right)
    # np.savetxt('avg_per_fundus_' + fmethod + '_sd_distances_from_border_right.csv',
    #            mean_border_to_feature_sd_distances_right)
