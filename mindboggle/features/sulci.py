#!/usr/bin/env python
"""
Functions to extract sulci from folds.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Forrest Sheng Bao, 2012  (forrest.bao@gmail.com)  http://fsbao.net

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#===============================================================================
# Extract sulci
#===============================================================================
def extract_sulci(labels_file, folds, label_pair_lists,
                  min_boundary=1, sulcus_names=[], save_file=False):
    """
    Identify sulci from folds in a brain surface according to a labeling
    protocol that includes a list of label pairs defining each sulcus.

    Definitions ::
        ``fold``: a group of connected, deep vertices
        ``label list``: list of labels used to define a single sulcus
        ``label pair list``: list of pairs of labels, where each pair
                             defines a boundary between two labeled regions;
                             no two label pair lists share a label pair
        ``sulcus ID``: ID to uniquely identify a sulcus,
                       as index to a sulcus label pair list

    Steps for each fold::
        1. Remove fold if it has fewer than two labels.
        2. Remove fold if its labels do not contain a sulcus label pair.
        3. Find vertices with labels that are in only one of the fold's
           label boundary pairs. Assign the vertices the sulcus with the
           label pair if they are connected to the label boundary for that pair.
        4. If there are remaining vertices, segment into sets of vertices
           connected to label boundaries, and assign a sulcus ID to each segment.

    Parameters
    ----------
    labels_file : string
        file name for surface mesh vtk containing labels for all vertices
    folds : list or string
        fold ID for each vertex or name of folds file containing folds scalars
    label_pair_lists : list of sublists of subsublists of integers
        each subsublist contains a pair of labels, and the sublist of these
        label pairs represents the label boundaries defining a sulcus
    min_boundary : integer
        minimum number of vertices for a sulcus label boundary segment
    sulcus_names : list of strings (optional)
        names of sulci
    save_file : Boolean
        save output VTK file?

    Returns
    -------
    sulci : array of integers
        sulcus numbers for all vertices (-1 for non-sulcus vertices)
    n_sulci : integers
        number of sulci
    sulci_file : string (if save_file)
        name of output VTK file with sulcus numbers (-1 for non-sulcus vertices)

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    >>> from mindboggle.labels.protocol.sulci_labelpairs_DKT import sulcus_boundaries
    >>> from mindboggle.utils.mesh import find_neighbors_from_file
    >>> from mindboggle.features.sulci import extract_sulci
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>>
    >>> # Load labels, folds, neighbor lists, and sulcus names and label pairs
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> folds, name = read_scalars(folds_file)
    >>> label_pair_lists = sulcus_boundaries()
    >>> min_boundary = 10
    >>> sulcus_names_file = os.path.join(path, 'info', 'sulcus_names.txt')
    >>> fid = open(sulcus_names_file, 'r')
    >>> sulcus_names = fid.readlines()
    >>> sulcus_names = [x.strip('\n') for x in sulcus_names]
    >>>
    >>> # Extract sulci
    >>> sulci, n_sulci = extract_sulci(labels_file, folds, label_pair_lists,
    >>>                                min_boundary, sulcus_names)
    >>>
    >>> # Write sulci to a new VTK file and view:
    >>> rewrite_scalars(labels_file, 'test_extract_sulci.vtk', sulci, 'sulci', sulci)
    >>> from mindboggle.utils.mesh import plot_vtk
    >>> plot_vtk('test_extract_sulci.vtk')

    """
    import os
    from time import time
    import numpy as np
    from mindboggle.utils.io_vtk import read_scalars, read_vtk, rewrite_scalars
    from mindboggle.utils.mesh import find_neighbors
    from mindboggle.labels.label import extract_borders
    from mindboggle.labels.segment import propagate, segment


    if isinstance(folds, str):
        folds, name = read_scalars(folds)

    # Prepare list of sulcus label lists (labels within a sulcus)
    label_lists = []
    for row in label_pair_lists:
        label_lists.append(list(np.unique(np.asarray(row))))

    # Prepare list of all unique sorted label pairs in the labeling protocol
    protocol_pairs = []
    [protocol_pairs.append(np.unique(x).tolist()) for lst in label_pair_lists
     for x in lst if np.unique(x).tolist() not in protocol_pairs]

    # Load points, faces, and neighbors
    faces, lines, indices, points, npoints, labels, name = read_vtk(labels_file)
    neighbor_lists = find_neighbors(faces, npoints)

    # Array of sulcus IDs for fold vertices, initialized as -1.
    # Since we do not touch gyral vertices and vertices whose labels
    # are not in the label list, or vertices having only one label,
    # their sulcus IDs will remain -1.
    sulci = -1 * np.ones(npoints)

    #---------------------------------------------------------------------------
    # Loop through folds
    #---------------------------------------------------------------------------
    fold_numbers = [int(x) for x in np.unique(folds) if x > -1]
    n_folds = len(fold_numbers)
    print("Extract sulci from {0} folds...".format(n_folds))
    t0 = time()
    for n_fold in fold_numbers:
        fold = [i for i,x in enumerate(folds) if x == n_fold]
        len_fold = len(fold)
        # List the labels in this fold (greater than zero)
        fold_labels = [labels[x] for x in fold]
        unique_fold_labels = [int(x) for x in np.unique(fold_labels) if x > 0]

        #-----------------------------------------------------------------------
        # NO MATCH -- fold has fewer than two labels
        #-----------------------------------------------------------------------
        if len(unique_fold_labels) < 2:
            # Ignore: sulci already initialized with -1 values
            if not unique_fold_labels:
                print("  Fold {0} ({1} vertices): NO MATCH -- fold has no labels".
                      format(n_fold, len_fold))
            else:
                print("  Fold {0} ({1} vertices): "
                  "NO MATCH -- fold has only one label ({2})".
                  format(n_fold, len_fold, unique_fold_labels[0]))
            # Ignore: sulci already initialized with -1 values

        else:
            # Find all label boundary pairs within the fold
            indices_fold_pairs, fold_pairs, unique_fold_pairs = extract_borders(
                fold, labels, neighbor_lists)

            # Find fold label pairs in the protocol (pairs are already sorted)
            fold_pairs_in_protocol = [x for x in unique_fold_pairs
                                      if x in protocol_pairs]

            print("  Fold {0} labels: {1} ({2} vertices)".format(n_fold,
                  ', '.join([str(x) for x in unique_fold_labels]), len_fold))
            print("  Fold {0} label pairs in protocol: {1}".format(n_fold,
                  ', '.join([str(x) for x in fold_pairs_in_protocol])))

            #-------------------------------------------------------------------
            # NO MATCH -- fold has no sulcus label pair
            #-------------------------------------------------------------------
            if not len(fold_pairs_in_protocol):
                print("  Fold {0}: NO MATCH -- fold has no sulcus label pair".
                      format(n_fold, len_fold))

            #-------------------------------------------------------------------
            # Possible matches
            #-------------------------------------------------------------------
            else:
                # Labels in the protocol (includes repeats across label pairs)
                labels_in_pairs = [x for lst in fold_pairs_in_protocol for x in lst]

                # Labels that appear in one or more than one sulcus label boundary
                unique_labels = []
                nonunique_labels = []
                for label in np.unique(labels_in_pairs):
                    if len([x for x in labels_in_pairs if x == label]) == 1:
                        unique_labels.append(label)
                    else:
                        nonunique_labels.append(label)

                #---------------------------------------------------------------
                # Vertices whose labels are in only one sulcus label pair
                #---------------------------------------------------------------
                # Find vertices with a label that is in only one of the fold's
                # label pairs (the other label in the pair can exist
                # in other pairs). Assign the vertices the sulcus with the label
                # pair if they are connected to the label boundary for that pair.
                #---------------------------------------------------------------
                if len(unique_labels):

                    for pair in fold_pairs_in_protocol:
                        # If one or both labels in label pair is/are unique
                        unique_labels_in_pair = [x for x in pair if x in unique_labels]
                        n_unique = len(unique_labels_in_pair)
                        if n_unique:

                            ID = [i for i,x in enumerate(label_pair_lists)
                                  if pair in x][0]

                            # Construct seeds from label boundary vertices
                            # (fold_pairs and pair already sorted)
                            indices_pair = [x for i,x in enumerate(indices_fold_pairs)
                                            if fold_pairs[i] == pair]

                            # Identify vertices with unique label(s) in pair
                            indices_unique_labels = [fold[i]
                                                     for i,x in enumerate(fold_labels)
                                                     if x in unique_labels_in_pair]

                            # Propagate from seeds to labels in label pair
                            sulci2 = segment(indices_unique_labels, neighbor_lists,
                                             min_region_size=1,
                                             seed_lists=[indices_pair],
                                             keep_seeding=False,
                                             spread_within_labels=True,
                                             labels=labels)
                            sulci[sulci2 > -1] = ID

                            # Print statement
                            if n_unique == 1:
                                ps1 = '1 label'
                            else:
                                ps1 = 'Both labels'
                            if len(sulcus_names):
                                ps2 = sulcus_names[ID]
                            else:
                                ps2 = ''
                            print("    {0} unique to one fold pair: {1} {2}".
                                  format(ps1, ps2, unique_labels_in_pair))

                #---------------------------------------------------------------
                # Vertex labels shared by multiple label pairs
                #---------------------------------------------------------------
                # Propagate labels from label boundaries to vertices with labels
                # that are shared by multiple label pairs in the fold.
                #---------------------------------------------------------------
                if len(nonunique_labels):
                    # For each label shared by different label pairs
                    for label in nonunique_labels:
                        # Print statement
                        print("    Propagate sulcus label boundaries with label {0}".
                              format(int(label)))

                        # Construct seeds from label boundary vertices
                        seeds = -1 * np.ones(len(points))
                        for ID, label_pair_list in enumerate(label_pair_lists):
                            label_pairs = [x for x in label_pair_list if label in x]
                            for label_pair in label_pairs:
                                indices_pair = [x for i,x in enumerate(indices_fold_pairs)
                                    if np.sort(fold_pairs[i]).tolist() == label_pair]
                                if len(indices_pair):

                                    # Do not include short boundary segments
                                    if min_boundary > 1:
                                        indices_pair2 = []
                                        seeds2 = segment(indices_pair, neighbor_lists)
                                        for seed2 in range(int(max(seeds2))+1):
                                            iseed2 = [i for i,x in enumerate(seeds2)
                                                      if x == seed2]
                                            if len(iseed2) >= min_boundary:
                                                indices_pair2.extend(iseed2)
                                            else:
                                                if len(iseed2) == 1:
                                                    print("    Remove assignment "
                                                          "of ID {0} from 1 vertex".
                                                          format(seed2))
                                                else:
                                                    print("    Remove assignment "
                                                          "of ID {0} from {1} vertices".
                                                          format(seed2, len(iseed2)))
                                        indices_pair = indices_pair2

                                    # Assign sulcus IDs to seeds
                                    seeds[indices_pair] = ID

                        # Identify vertices with the label
                        label_array = -1 * np.ones(len(points))
                        indices_label = [fold[i] for i,x in enumerate(fold_labels)
                                         if x == label]
                        if len(indices_label):
                            label_array[indices_label] = 1

                            # Propagate from seeds to vertices with label
                            #indices_seeds = []
                            #for seed in range(int(max(seeds))+1):
                            #    indices_seeds.append([i for i,x in enumerate(seeds)
                            #                          if x == seed])
                            #sulci2 = segment(indices_label, neighbor_lists,
                            #                 50, indices_seeds, False, True, labels)
                            sulci2 = propagate(points, faces, label_array,
                                               seeds, sulci, max_iters=10000,
                                               tol=0.001, sigma=5)
                            sulci[sulci2 > -1] = sulci2[sulci2 > -1]

    #---------------------------------------------------------------------------
    # Print out assigned sulci
    #---------------------------------------------------------------------------
    sulcus_numbers = [int(x) for x in np.unique(sulci) if x > -1]
    n_sulci = len(sulcus_numbers)
    print("Extracted {0} sulci from {1} folds ({2:.1f}s):".
          format(n_sulci, n_folds, time()-t0))
    if len(sulcus_names):
        for sulcus_number in sulcus_numbers:
            print("  {0}: {1}".format(sulcus_number, sulcus_names[sulcus_number]))
    else:
        print("  " + ", ".join([str(x) for x in sulcus_numbers]))

    #---------------------------------------------------------------------------
    # Print out unresolved sulci
    #---------------------------------------------------------------------------
    unresolved = [i for i in range(len(label_pair_lists))
                  if i not in sulcus_numbers]
    if len(unresolved) == 1:
        print("The following sulcus is unaccounted for:")
    else:
        print("The following {0} sulci are unaccounted for:".format(len(unresolved)))
    if len(sulcus_names):
        for sulcus_number in unresolved:
            print("  {0}: {1}".format(sulcus_number, sulcus_names[sulcus_number]))
    else:
        print("  " + ", ".join([str(x) for x in unresolved]))

    #---------------------------------------------------------------------------
    # Return sulci, number of sulci, and file name
    #---------------------------------------------------------------------------
    if save_file:

        sulci_file = os.path.join(os.getcwd(), 'sulci.vtk')
        rewrite_scalars(labels_file, sulci_file, sulci, 'sulci', sulci)

    else:
        sulci_file = None

    return sulci, n_sulci, sulci_file

#===============================================================================
# Example: extract_sulci()
#===============================================================================
if __name__ == "__main__":

    import os
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    from mindboggle.labels.protocol.sulci_labelpairs_DKT import sulcus_boundaries
    from mindboggle.utils.mesh import find_neighbors_from_file
    from mindboggle.features.sulci import extract_sulci
    from mindboggle.utils.mesh import plot_vtk
    path = os.environ['MINDBOGGLE_DATA']

    # Load labels, folds, neighbor lists, and sulcus names and label pairs
    labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    folds_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    folds, name = read_scalars(folds_file)
    neighbor_lists = find_neighbors(labels_file)
    label_pair_lists = sulcus_boundaries()
    min_boundary = 10
    sulcus_names_file = os.path.join(path, 'info', 'sulcus_names.txt')
    fid = open(sulcus_names_file, 'r')
    sulcus_names = fid.readlines()
    sulcus_names = [x.strip('\n') for x in sulcus_names]

    # Extract sulci
    sulci, n_sulci = extract_sulci(labels_file, folds, neighbor_lists,
                                   label_pair_lists, min_boundary, sulcus_names)

    # Finally, write points, faces and sulci to a new vtk file and view:
    rewrite_scalars(labels_file, 'test_extract_sulci.vtk', sulci, 'sulci', sulci)
    plot_vtk('test_extract_sulci.vtk')
