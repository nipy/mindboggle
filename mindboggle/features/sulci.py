#!/usr/bin/env python
"""
Functions to extract sulci from folds.

Authors:
    - Arno Klein, 2012-2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def extract_sulci(labels_file, folds_or_file, hemi, min_boundary=1,
                  sulcus_names=[], background_value=-1, verbose=False):
    """
    Identify sulci from folds in a brain surface according to a labeling
    protocol that includes a list of label pairs defining each sulcus.

    A fold is a group of connected, deep vertices.

    Steps for each fold ::

        1. Remove fold if it has fewer than two labels.
        2. Remove fold if its labels do not contain a sulcus label pair.
        3. Find vertices with labels that are in only one of the fold's
           label boundary pairs. Assign the vertices the sulcus with the label
           pair if they are connected to the label boundary for that pair.
        4. If there are remaining vertices, segment into sets of vertices
           connected to label boundaries, and assign a unique ID to each set.

    Parameters
    ----------
    labels_file : string
        file name for surface mesh VTK containing labels for all vertices
    folds_or_file : numpy array, list or string
        fold number for each vertex / name of VTK file containing fold scalars
    hemi : string
        hemisphere abbreviation in {'lh', 'rh'} for sulcus labels
    min_boundary : integer
        minimum number of vertices for a sulcus label boundary segment
    sulcus_names : list of strings
        names of sulci
    background_value : integer or float
        background value
    verbose : bool
        print statements?

    Returns
    -------
    sulci : list of integers
        sulcus numbers for all vertices (-1 for non-sulcus vertices)
    n_sulci : integers
        number of sulci
    sulci_file : string
        output VTK file with sulcus numbers (-1 for non-sulcus vertices)

    Examples
    --------
    >>> import numpy as np
    >>> from mindboggle.features.sulci import extract_sulci
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> # Load labels, folds, neighbor lists, and sulcus names and label pairs
    >>> labels_file = fetch_data(urls['left_freesurfer_labels'])
    >>> folds_file = fetch_data(urls['left_folds'])
    >>> folds_or_file, name = read_scalars(folds_file, True, True)
    >>> # Limit number of folds to speed up the test:
    >>> background_value = -1
    >>> limit_folds = True
    >>> if limit_folds:
    ...     fold_numbers = [4, 7] #[4, 6]
    ...     i0 = [i for i,x in enumerate(folds_or_file) if x not in fold_numbers]
    ...     folds_or_file[i0] = background_value
    >>> hemi = 'lh'
    >>> min_boundary = 10
    >>> sulcus_names = []
    >>> verbose = True #False
    >>> sulci, n_sulci, sulci_file = extract_sulci(labels_file, folds_or_file,
    ...     hemi, min_boundary, sulcus_names, background_value, verbose)
    >>> n_sulci  # 23
    3
    >>> lens = [len([x for x in sulci if x==y])
    ...         for y in np.unique(sulci) if y != -1]
    >>> lens[0:10]  # [6358, 3288, 7612, 5205, 4414, 6251, 3493, 2566, 4436, 739]
    [1151]

    View result (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces # doctest: +SKIP
    >>> from mindboggle.mio.vtks import rewrite_scalars # doctest: +SKIP
    >>> rewrite_scalars(sulci_file, 'just_sulci.vtk', sulci, 'sulci', sulci) # doctest: +SKIP
    >>> plot_surfaces('just_sulci.vtk') # doctest: +SKIP

    """
    import os
    from time import time
    import numpy as np

    from mindboggle.mio.vtks import read_scalars, read_vtk, rewrite_scalars
    from mindboggle.guts.mesh import find_neighbors
    from mindboggle.guts.segment import extract_borders, propagate, segment_regions
    from mindboggle.mio.labels import DKTprotocol

    # Load fold numbers if folds_or_file is a string:
    if isinstance(folds_or_file, str):
        folds, name = read_scalars(folds_or_file)
    elif isinstance(folds_or_file, list):
        folds = folds_or_file
    elif isinstance(folds_or_file, np.ndarray):
        folds = folds_or_file.tolist()

    dkt = DKTprotocol()

    if hemi == 'lh':
        pair_lists = dkt.left_sulcus_label_pair_lists
    elif hemi == 'rh':
        pair_lists = dkt.right_sulcus_label_pair_lists
    else:
        raise IOError("Warning: hemisphere not properly specified ('lh' or 'rh').")

    # Load points, faces, and neighbors:
    points, indices, lines, faces, labels, scalar_names, npoints, \
            input_vtk = read_vtk(labels_file)
    neighbor_lists = find_neighbors(faces, npoints)

    # Array of sulcus IDs for fold vertices, initialized as -1.
    # Since we do not touch gyral vertices and vertices whose labels
    # are not in the label list, or vertices having only one label,
    # their sulcus IDs will remain -1:
    sulci = background_value * np.ones(npoints)

    #-------------------------------------------------------------------------
    # Loop through folds
    #-------------------------------------------------------------------------
    fold_numbers = [int(x) for x in np.unique(folds) if x != background_value]
    n_folds = len(fold_numbers)
    if verbose:
        print("Extract sulci from {0} folds...".format(n_folds))
    t0 = time()
    for n_fold in fold_numbers:
        fold_indices = [i for i,x in enumerate(folds) if x == n_fold]
        len_fold = len(fold_indices)

        # List the labels in this fold:
        fold_labels = [labels[x] for x in fold_indices]
        unique_fold_labels = [int(x) for x in np.unique(fold_labels)
                              if x != background_value]

        #---------------------------------------------------------------------
        # NO MATCH -- fold has fewer than two labels
        #---------------------------------------------------------------------
        if verbose and len(unique_fold_labels) < 2:
            # Ignore: sulci already initialized with -1 values:
            if not unique_fold_labels:
                print("  Fold {0} ({1} vertices): "
                      "NO MATCH -- fold has no labels".
                      format(n_fold, len_fold))
            else:
                print("  Fold {0} ({1} vertices): "
                  "NO MATCH -- fold has only one label ({2})".
                  format(n_fold, len_fold, unique_fold_labels[0]))
            # Ignore: sulci already initialized with -1 values

        else:
            # Find all label boundary pairs within the fold:
            indices_fold_pairs, fold_pairs, unique_fold_pairs = \
                extract_borders(fold_indices, labels, neighbor_lists,
                                ignore_values=[], return_label_pairs=True)

            # Find fold label pairs in the protocol (pairs are already sorted):
            fold_pairs_in_protocol = [x for x in unique_fold_pairs
                                      if x in dkt.unique_sulcus_label_pairs]

            if verbose and unique_fold_labels:
                print("  Fold {0} labels: {1} ({2} vertices)".format(n_fold,
                      ', '.join([str(x) for x in unique_fold_labels]),
                      len_fold))
            #-----------------------------------------------------------------
            # NO MATCH -- fold has no sulcus label pair
            #-----------------------------------------------------------------
            if verbose and not fold_pairs_in_protocol:
                print("  Fold {0}: NO MATCH -- fold has no sulcus label pair".
                      format(n_fold, len_fold))

            #-----------------------------------------------------------------
            # Possible matches
            #-----------------------------------------------------------------
            else:
                if verbose:
                    print("  Fold {0} label pairs in protocol: {1}".format(n_fold,
                          ', '.join([str(x) for x in fold_pairs_in_protocol])))

                # Labels in the protocol (includes repeats across label pairs):
                labels_in_pairs = [x for lst in fold_pairs_in_protocol
                                   for x in lst]

                # Labels that appear in one or more sulcus label boundary:
                unique_labels = []
                nonunique_labels = []
                for label in np.unique(labels_in_pairs):
                    if len([x for x in labels_in_pairs if x == label]) == 1:
                        unique_labels.append(label)
                    else:
                        nonunique_labels.append(label)

                #-------------------------------------------------------------
                # Vertices whose labels are in only one sulcus label pair
                #-------------------------------------------------------------
                # Find vertices with a label that is in only one of the fold's
                # label pairs (the other label in the pair can exist in other
                # pairs). Assign the vertices the sulcus with the label pair
                # if they are connected to the label boundary for that pair.
                #-------------------------------------------------------------
                if unique_labels:

                    for pair in fold_pairs_in_protocol:

                        # If one or both labels in label pair is/are unique:
                        unique_labels_in_pair = [x for x in pair
                                                 if x in unique_labels]
                        n_unique = len(unique_labels_in_pair)
                        if n_unique:

                            ID = None
                            for i, pair_list in enumerate(pair_lists):
                                if not isinstance(pair_list, list):
                                    pair_list = [pair_list]
                                if pair in pair_list:
                                    ID = i
                                    break
                            if ID:
                                # Seeds from label boundary vertices
                                # (fold_pairs and pair already sorted):
                                indices_pair = [x for i,x
                                    in enumerate(indices_fold_pairs)
                                    if fold_pairs[i] == pair]

                                # Vertices with unique label(s) in pair:
                                indices_unique_labels = [fold_indices[i]
                                     for i,x in enumerate(fold_labels)
                                     if x in unique_labels_in_pair]
                                             #dkt.unique_sulcus_label_pairs]

                                # Propagate sulcus ID from seeds to vertices
                                # with "unique" labels (only exist in one
                                # label pair in a fold); propagation ensures
                                # that sulci consist of contiguous vertices
                                # for each label boundary:
                                sulci2 = segment_regions(indices_unique_labels,
                                         neighbor_lists,
                                         min_region_size=1,
                                         seed_lists=[indices_pair],
                                         keep_seeding=False,
                                         spread_within_labels=True,
                                         labels=labels,
                                         label_lists=[],
                                         values=[], max_steps='',
                                         background_value=background_value,
                                         verbose=False)

                                sulci[sulci2 != background_value] = ID

                                # Print statement:
                                if verbose:
                                    if n_unique == 1:
                                        ps1 = 'One label'
                                    else:
                                        ps1 = 'Both labels'
                                    if len(sulcus_names):
                                        ps2 = sulcus_names[ID]
                                    else:
                                        ps2 = ''
                                    print("    {0} unique to one fold pair: "
                                          "{1} {2}".
                                          format(ps1, ps2,
                                                 unique_labels_in_pair))

                #-------------------------------------------------------------
                # Vertex labels shared by multiple label pairs
                #-------------------------------------------------------------
                # Propagate labels from label borders to vertices with labels
                # that are shared by multiple label pairs in the fold.
                #-------------------------------------------------------------
                if len(nonunique_labels):
                    # For each label shared by different label pairs:
                    for label in nonunique_labels:
                        # Print statement:
                        if verbose:
                            print("    Propagate sulcus borders with label {0}".
                                  format(int(label)))

                        # Construct seeds from label boundary vertices:
                        seeds = background_value * np.ones(npoints)

                        for ID, pair_list in enumerate(pair_lists):
                            if not isinstance(pair_list, list):
                                pair_list = [pair_list]
                            label_pairs = [x for x in pair_list if label in x]
                            for label_pair in label_pairs:
                                indices_pair = [x for i,x
                                    in enumerate(indices_fold_pairs)
                                    if np.sort(fold_pairs[i]).
                                    tolist() == label_pair]
                                if indices_pair:

                                    # Do not include short boundary segments:
                                    if min_boundary > 1:
                                        indices_pair2 = []
                                        seeds2 = segment_regions(indices_pair,
                                                    neighbor_lists, 1, [],
                                                    False, False, [], [],
                                                    [], '', background_value,
                                                    verbose)

                                        useeds2 = [x for x in
                                                   np.unique(seeds2)
                                                   if x != background_value]
                                        for seed2 in useeds2:
                                            iseed2 = [i for i,x
                                                      in enumerate(seeds2)
                                                      if x == seed2]
                                            if len(iseed2) >= min_boundary:
                                                indices_pair2.extend(iseed2)
                                            elif verbose:
                                                if len(iseed2) == 1:
                                                    print("    Remove "
                                                          "assignment "
                                                          "of ID {0} from "
                                                          "1 vertex".
                                                          format(seed2))
                                                else:
                                                    print("    Remove "
                                                          "assignment "
                                                          "of ID {0} from "
                                                          "{1} vertices".
                                                          format(seed2,
                                                                 len(iseed2)))
                                        indices_pair = indices_pair2

                                    # Assign sulcus IDs to seeds:
                                    seeds[indices_pair] = ID

                        # Identify vertices with the label:
                        indices_label = [fold_indices[i] for i,x
                                         in enumerate(fold_labels)
                                         if x == label]
                        if len(indices_label):

                            # Propagate sulcus ID from seeds to vertices
                            # with a given shared label:
                            seg_vs_prop = True
                            if seg_vs_prop:
                                indices_seeds = []
                                for seed in np.unique(seeds):
                                   indices_seeds.append([i for i,x
                                                         in enumerate(seeds)
                                                         if x == seed])
                                sulci2 = segment_regions(indices_label,
                                            neighbor_lists, 50, indices_seeds,
                                            False, True, labels, [], [], '',
                                            background_value, verbose)
                            else:
                                label_array = background_value * \
                                              np.ones(npoints)
                                label_array[indices_label] = 1
                                sulci2 = propagate(points, faces,
                                            label_array, seeds, sulci,
                                            max_iters=10000,
                                            tol=0.001, sigma=5,
                                            background_value=background_value,
                                            verbose=verbose)
                            sulci[sulci2 != background_value] = \
                                sulci2[sulci2 != background_value]

    sulcus_numbers = [int(x) for x in np.unique(sulci)
                      if x != background_value]
    n_sulci = len(sulcus_numbers)

    #-------------------------------------------------------------------------
    # Print statements
    #-------------------------------------------------------------------------
    if verbose:
        if n_sulci == 1:
            sulcus_str = 'sulcus'
        else:
            sulcus_str = 'sulci'
        if n_folds == 1:
            folds_str = 'fold'
        else:
            folds_str = 'folds'
        print("Extracted {0} {1} from {2} {3} ({4:.1f}s):".
                  format(n_sulci, sulcus_str, n_folds, folds_str, time()-t0))
        if sulcus_names:
            for sulcus_number in sulcus_numbers:
                print("  {0}: {1}".format(sulcus_number,
                                          sulcus_names[sulcus_number]))
        elif sulcus_numbers:
            print("  " + ", ".join([str(x) for x in sulcus_numbers]))

        unresolved = [i for i in range(len(pair_lists))
                      if i not in sulcus_numbers]
        if len(unresolved) == 1:
            print("The following sulcus is unaccounted for:")
        else:
            print("The following {0} sulci are unaccounted for:".
                  format(len(unresolved)))
        if sulcus_names:
            for sulcus_number in unresolved:
                print("  {0}: {1}".format(sulcus_number,
                                          sulcus_names[sulcus_number]))
        else:
            print("  " + ", ".join([str(x) for x in unresolved]))

    #-------------------------------------------------------------------------
    # Return sulci, number of sulci, and file name
    #-------------------------------------------------------------------------
    sulci = [int(x) for x in sulci]
    sulci_file = os.path.join(os.getcwd(), 'sulci.vtk')
    rewrite_scalars(labels_file, sulci_file, sulci, 'sulci', [],
                    background_value)

    if not os.path.exists(sulci_file):
        raise IOError(sulci_file + " not found")

    return sulci, n_sulci, sulci_file


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)  # py.test --doctest-modules