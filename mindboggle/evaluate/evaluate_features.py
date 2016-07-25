#!/usr/bin/env python
"""
Evaluate deep surface features by computing the minimum distance from each
label border vertex to all of the feature vertices in the same sulcus.
The label borders run along the deepest parts of many sulci and
correspond to fundi in the DKT cortical labeling protocol.

Results have been saved to https://osf.io/r95wb/

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
    output_vtk_name : bool
        if not empty, output a VTK file beginning with output_vtk_name that
        contains a surface with mean distances as scalars
    verbose : bool
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
    from mindboggle.mio.vtks import read_vtk, read_scalars, write_vtk
    from mindboggle.guts.mesh import find_neighbors, keep_faces
    from mindboggle.guts.segment import extract_borders
    from mindboggle.guts.compute import source_to_target_distances
    from mindboggle.mio.labels import DKTprotocol

    dkt = DKTprotocol()
    # ------------------------------------------------------------------------
    # Load labels, features, and sulci:
    # ------------------------------------------------------------------------
    points, indices, lines, faces, labels, scalar_names, npoints, \
    input_vtk = read_vtk(labels_file, True, True)
    features, name = read_scalars(features_file, True, True)
    if sulci_file:
        sulci, name = read_scalars(sulci_file, True, True)
        # List of indices to sulcus vertices:
        sulcus_indices = [i for i,x in enumerate(sulci)
                          if x not in excludeIDs]
        segmentIDs = sulci
        sulcus_faces = keep_faces(faces, sulcus_indices)
    else:
        sulcus_indices = list(range(len(labels)))
        segmentIDs = []
        sulcus_faces = faces

    # ------------------------------------------------------------------------
    # Prepare neighbors, label pairs, border IDs, and outputs:
    # ------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------
    # Loop through sulci:
    # ------------------------------------------------------------------------
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

        # --------------------------------------------------------------------
        # Construct a feature-to-border distance matrix and VTK file:
        # --------------------------------------------------------------------
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

        # --------------------------------------------------------------------
        # Construct a border-to-feature distance matrix and VTK file:
        # --------------------------------------------------------------------
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

    # ------------------------------------------------------------------------
    # Return outputs:
    # ------------------------------------------------------------------------
    return feature_to_border_mean_distances, feature_to_border_sd_distances,\
           feature_to_border_distances_vtk,\
           border_to_feature_mean_distances, border_to_feature_sd_distances,\
           border_to_feature_distances_vtk


# ----------------------------------------------------------------------------
# Run evaluate_deep_features() on fundi extracted from Mindboggle-101 data
# by Mindboggle, and Forrest Bao's, Gang Li's, and Olivier Coulon's methods.
# ----------------------------------------------------------------------------
if __name__ == "__main__":

    import os
    import numpy as np
    import pandas as pd

    from mindboggle.mio.labels import DKTprotocol
    from mindboggle.evaluate.evaluate_features import evaluate_deep_features

    # For plotting:
    from math import pi
    from bokeh.models import HoverTool
    from bokeh.plotting import ColumnDataSource, figure, show, save, output_file
    from mindboggle.mio.colors import viridis_colormap

    measure_feature_distances = False
    maxd = 53
    feature_dir = '/homedir/fundus_evaluation_2014/fundi_vtk'

    plot_and_describe_fundus_distances = False
    use_all_subjects = True
    exclude_sulci = [20] # Sulcus 20 removed from protocol since initial run
    colors = viridis_colormap()
    #from matplotlib import cm as cmaps
    #import matplotlib.pyplot as plt
    #plt.register_cmap(name='viridis', cmap=cmaps.viridis)
    #plt.set_cmap(cmaps.viridis)

    compute_results = True

    mindboggled = '/mnt/nfs-share/Mindboggle101/mindboggled/manual'
    labels_dir = '/mnt/nfs-share/Mindboggle101/mindboggled/manual'

    # --------------------------------------------------------------------
    # Subjects:
    # --------------------------------------------------------------------
    # Select only those subjects for which all methods were applied:
    names = ['OASIS-TRT-20', 'MMRR-21', 'NKI-RS-22', 'NKI-TRT-20',
             'Afterthought', 'Colin27', 'Twins-2', 'MMRR-3T7T-2', 'HLN-12']
    numbers = [20,21,22,20, 1,1,2,2,12]
    scale_rect = 20
    # if not use_all_subjects:
    #     names = ['OASIS-TRT-20', 'MMRR-21', 'NKI-TRT-20']
    #     numbers = [20,21,20]
    #     scale_rect = 40

    nsubjects = sum(numbers)
    subjects = []
    for iname, name in enumerate(names):
        for n in range(1, numbers[iname]+1):
            subjects.append(name+'-'+str(n))

    # ------------------------------------------------------------------------
    # Alternating left, right cortex label numbers:
    # ------------------------------------------------------------------------
    dkt = DKTprotocol()
    labels = [str(x) for x in range(len(dkt.sulcus_names))]
    label_names = dkt.sulcus_names
    nsulci = len(label_names)
    #if exclude_sulci:
    #    nsulci = nsulci - len(exclude_sulci)

    # ------------------------------------------------------------------------
    # Measure feature distances:
    # ------------------------------------------------------------------------
    if measure_feature_distances:

        # --------------------------------------------------------------------
        # Set feature type ('fundi' or '' for every sulcus vertex), subjects:
        # --------------------------------------------------------------------
        feature_type = 'fundi' #'sulci'  # If 'fundi', select 'nmethod' below.
        if feature_type == 'fundi':
            # Features: fundus method:
            # 0 = mindboggle
            # 1 = forrest scalars
            # 2 = forrest lines
            # 3 = gang li
            # 4 = olivier coulon
            nmethod = 0
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

        # --------------------------------------------------------------------
        # Miscellaneous defaults:
        # --------------------------------------------------------------------
        surfs = ['left_cortical_surface', 'right_cortical_surface']
        hemis = ['lh', 'rh']

        # --------------------------------------------------------------------
        # Loop through subjects and hemispheres:
        # --------------------------------------------------------------------
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
                    # --------------------------------------------------------
                    # Identify surface files with labels and with sulci:
                    # --------------------------------------------------------
                    mdir = os.path.join(mindboggled, subject)
                    ldir = os.path.join(labels_dir, subject)
                    sulci_file = os.path.join(mdir, 'features', surf, 'sulci.vtk')
                    labels_file = os.path.join(ldir, 'labels', surf,
                                               'relabeled_labels.DKT31.manual.vtk')
                    # --------------------------------------------------------
                    # Identify features file:
                    # --------------------------------------------------------
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

                    # --------------------------------------------------------
                    # Compute distances between features and label borders
                    # in sulci corresponding to fundi:
                    # --------------------------------------------------------
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

        # --------------------------------------------------------------------
        # Save tables of mean distances:
        # --------------------------------------------------------------------
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

        # # --------------------------------------------------------------------
        # # Save tables of mean distances averaged across all subjects:
        # # NOTE: np.mean() results in nan's if any element has a nan.
        # # --------------------------------------------------------------------
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


    # ------------------------------------------------------------------------
    # Plot and describe fundus distances:
    # ------------------------------------------------------------------------
    if plot_and_describe_fundus_distances:

        # --------------------------------------------------------------------
        # Features:
        # --------------------------------------------------------------------
        # Features: fundus method:
        # 0 = mindboggle
        # 1 = forrest scalars
        # 2 = forrest lines
        # 3 = gang li
        # 4 = olivier coulon
        fmethods = ['mindboggle_fundi',
                    'ForrestBao_fundi',
                    'GangLi_fundi',
                    'OlivierCoulon_fundi']

        # --------------------------------------------------------------------
        # File names, paths:
        # --------------------------------------------------------------------
        mindboggled = '/mnt/nfs-share/Mindboggle101/mindboggled/auto'
        fmethod_dirs = [mindboggled,
                        os.path.join(feature_dir, fmethods[0]),
                        os.path.join(feature_dir, fmethods[1]),
                        os.path.join(feature_dir, fmethods[2]),
                        os.path.join(feature_dir, fmethods[3])]
        fmethod_names = ["Mindboggle's fundi",
                         "Forrest Bao's fundi",
                         "Gang Li's fundi",
                         "Olivier Coulon's_fundi"]

        feature_tables_dir = '/Users/arno/Data/tables_fundus_label_distances'

        # --------------------------------------------------------------------
        # Loop through methods:
        # --------------------------------------------------------------------
        for imethod, fmethod in enumerate(fmethods):
            fmethod_name = fmethod_names[imethod]

            # ----------------------------------------------------------------
            # Tables of mean distances:
            # ----------------------------------------------------------------
            desc1 = 'Mean distances from left ' + fmethod_name + ' to label borders'
            desc2 = 'Mean distances from left label borders to ' + fmethod_name
            desc3 = 'Mean distances from right ' + fmethod_name + ' to label borders'
            desc4 = 'Mean distances from right label borders to ' + fmethod_name
            desc5 = 'Standard deviation of distances from left ' + fmethod_name + ' to label borders'
            desc6 = 'Standard deviation of distances from left label borders to ' + fmethod_name
            desc7 = 'Standard deviation of distances from right ' + fmethod_name + ' to label borders'
            desc8 = 'Standard deviation of distances from right label borders to ' + fmethod_name
            descriptions = [desc1, desc2, desc3, desc4, desc5, desc6, desc7, desc8]
            table1 = fmethod + '_mean_distances_to_border_left'
            table2 = fmethod + '_mean_distances_from_border_left'
            table3 = fmethod + '_mean_distances_to_border_right'
            table4 = fmethod + '_mean_distances_from_border_right'
            table5 = fmethod + '_sd_distances_to_border_left'
            table6 = fmethod + '_sd_distances_from_border_left'
            table7 = fmethod + '_sd_distances_to_border_right'
            table8 = fmethod + '_sd_distances_from_border_right'
            tables = [table1, table2, table3, table4, table5, table6, table7, table8]
            feature_tables = [os.path.join(feature_tables_dir, x + '.csv')
                              for x in tables]

            # ----------------------------------------------------------------
            # Plot heatmap for labels X subjects array for each table:
            # ----------------------------------------------------------------
            for itable, feature_table in enumerate(feature_tables):
                title = descriptions[itable]
                table = tables[itable]

                summary_file = table + '_summary.csv'
                html_file = table + '_summary.html'

                data_file = feature_table
                data = pd.read_csv(data_file, sep=" ",
                                   index_col=False, header=None)
                if exclude_sulci:
                    data = data.iloc[:,
                           [x for x in range(len(labels) + len(exclude_sulci))
                            if x not in exclude_sulci]]
                print(title)
                print('Maximum distance = ' + str(data.max().max()))
                #maxd = data.max().max()
                data_summary = data.describe()
                data_summary.to_csv(summary_file)

                # Set up the data for plotting. We will need to have values for every
                # pair of subject/label names. Map the value to a color.
                subjectx = []
                labelx = []
                label_namex = []
                valuex = []
                colorx = []
                for ilabel, label in enumerate(labels):
                    for isubject, subject in enumerate(subjects):
                        labelx.append(label)
                        label_namex.append(label_names[ilabel])
                        subjectx.append(subject)
                        value = data.iloc[isubject, ilabel]
                        valuex.append(value)
                        if np.isnan(value):
                            rgb = [0, 0, 0]
                        else:
                            rgb = [np.int(255 * x) for x in
                                   colors[np.int(255 * value / maxd)]]
                        hex = "#%02x%02x%02x" % tuple(rgb)
                        colorx.append(hex)

                output_file(html_file, title="Distance between fundi and label borders")
                source = ColumnDataSource(dict(subject=subjectx, label=labelx,
                                               label_name=label_namex,
                                               color=colorx, value=valuex))
                TOOLS = "hover,save,pan,box_zoom,wheel_zoom"

                plot_width = len(subjects) * scale_rect
                plot_height = nsulci * scale_rect
                p = figure(title=title, x_range=subjects, y_range=list(reversed(label_names)),
                           plot_width=plot_width, plot_height=plot_height,
                           x_axis_location="above", tools=TOOLS)


                p.grid.grid_line_color = None
                p.axis.axis_line_color = None
                p.axis.major_tick_line_color = None
                p.axis.major_label_text_font_size = "10pt"
                p.axis.major_label_standoff = 0
                p.xaxis.major_label_orientation = pi/3
                p.rect(x="subject", y="label_name", width=1, height=1, source=source,
                       color="color", line_color=None)

                p.select_one(HoverTool).tooltips = [
                    ('subject', '@subject'),
                    ('label', '@label'),
                    ('label name', '@label_name'),
                    ('value', '@value'),
                ]

                #show(p)      # show the plot
                #import sys; sys.exit()
                save(p)      # save the plot

    # ------------------------------------------------------------------------
    # Compute results from descriptions of fundus distances:
    # ------------------------------------------------------------------------
    if compute_results:

        fmethods = ['mindboggle_fundi',
                    'ForrestBao_fundi',
                    'GangLi_fundi']
                    #'OlivierCoulon_fundi']

        # --------------------------------------------------------------------
        # Loop through methods:
        # --------------------------------------------------------------------
        mean_to_border_left = []
        mean_to_border_right = []
        mean_from_border_left = []
        mean_from_border_right = []
        sd_to_border_left = []
        sd_to_border_right = []
        sd_from_border_left = []
        sd_from_border_right = []
        # --------------------------------------------------------------------
        # Tables of mean distances:
        # --------------------------------------------------------------------
        table_lists = []
        for imethod, fmethod in enumerate(fmethods):
            mean_to_border_left.append(fmethod + '_mean_distances_to_border_left')
            mean_from_border_left.append(fmethod + '_mean_distances_from_border_left')
            mean_to_border_right.append(fmethod + '_mean_distances_to_border_right')
            mean_from_border_right.append(fmethod + '_mean_distances_from_border_right')
            sd_to_border_left.append(fmethod + '_sd_distances_to_border_left')
            sd_from_border_left.append(fmethod + '_sd_distances_from_border_left')
            sd_to_border_right.append(fmethod + '_sd_distances_to_border_right')
            sd_from_border_right.append(fmethod + '_sd_distances_from_border_right')
        table_lists = [mean_to_border_left, mean_from_border_left,
                       mean_to_border_right, mean_from_border_right,
                       sd_to_border_left, sd_from_border_left,
                       sd_to_border_right, sd_from_border_right]

        # --------------------------------------------------------------------
        # Sort methods by a given statistic (mean):
        # sort_values.shape = (8, 24, 3):
        #     8 tables (mean_distances_to_border_left, etc.),
        #     24 sulci (frontomarginal, etc.),
        #     3 methods (save the ascending sort order)
        # --------------------------------------------------------------------
        counts = np.zeros(len(fmethods))
        sort_values = np.zeros((len(table_lists), nsulci, 3))
        min_values = np.zeros((len(table_lists), nsulci))
        mean_values= np.zeros((len(table_lists), len(fmethods)))
        mean_max_values= np.zeros((len(table_lists), len(fmethods)))
        for itable_list, table_list in enumerate(table_lists):
            mean_tables = np.zeros((len(table_list), nsulci))
            max_tables = np.zeros((len(table_list), nsulci))
            for itable, table in enumerate(table_list):
                summary_file = table + '_summary.csv'
                data = pd.read_csv(summary_file, sep=",",
                                   index_col='Unnamed: 0')
                mean_tables[itable, :] = data.iloc[[1], :]
                max_tables[itable, :] = data.iloc[[7], :]
                if itable_list == 0:
                    counts[itable] = np.sum(nsubjects - data.values[0])
            sort_values[itable_list, :, :] = mean_tables.argsort(axis=0).transpose()
            min_values[itable_list, :] = mean_tables.argmin(axis=0)
            mean_values[itable_list, :] = mean_tables.mean(axis=1)
            mean_max_values[itable_list, :] = max_tables.mean(axis=1)

        # --------------------------------------------------------------------
        # Save table with sum of sort indices. Minimum values are best:
        # --------------------------------------------------------------------
        np.savetxt('sort_8tables_3methods-mb-bao-li.csv',
                   sort_values.sum(axis=1), delimiter=',')

        # --------------------------------------------------------------------
        # Save table with mean values. Minimum values are best:
        # --------------------------------------------------------------------
        np.savetxt('means_8tables_3methods-mb-bao-li.csv',
                   mean_values, delimiter=',')
        # Save results as separate tables:
        np.savetxt('mean-to-border_mb-bao-li.csv',
                   (mean_values[0] + mean_values[2]) / 2, delimiter=',')
        np.savetxt('mean-from-border_mb-bao-li.csv',
                   (mean_values[1] + mean_values[3]) / 2, delimiter=',')
        np.savetxt('sd-to-border_mb-bao-li.csv',
                   (mean_values[4] + mean_values[6]) / 2, delimiter=',')
        np.savetxt('sd-from-border_mb-bao-li.csv',
                   (mean_values[5] + mean_values[7]) / 2, delimiter=',')

        # --------------------------------------------------------------------
        # Save table with mean of maximum values. Minimum values are best:
        # --------------------------------------------------------------------
        np.savetxt('mean-of-max_8tables_3methods-mb-bao-li.csv',
                   mean_max_values, delimiter=',')
        # Save results as separate tables:
        np.savetxt('mean-of-max-to-border_mb-bao-li.csv',
                   (mean_max_values[0] + mean_max_values[2]) / 2, delimiter=',')
        np.savetxt('mean-of-max-from-border_mb-bao-li.csv',
                   (mean_max_values[1] + mean_max_values[3]) / 2, delimiter=',')
        np.savetxt('sd-mean-of-max-to-border_mb-bao-li.csv',
                   (mean_max_values[4] + mean_max_values[6]) / 2, delimiter=',')
        np.savetxt('sd-mean-of-max-from-border_mb-bao-li.csv',
                   (mean_max_values[5] + mean_max_values[7]) / 2, delimiter=',')

        # --------------------------------------------------------------------
        # Count how many sulci have minimum values for each method:
        # number_of_sort_values.shape = (8, 4, 3):
        #     8 tables (mean_distances_to_border_left, etc.),
        #     4 measures (mean, sd, min, max)
        #     3 methods (mindboggle, etc.)
        # --------------------------------------------------------------------
        number_of_min_values = 10 * np.ones((len(table_lists), len(fmethods)))
        for itable, table in enumerate(min_values):
            for imethod  in range(len(fmethods)):
                number_of_min_values[itable, imethod] = \
                    np.size(np.where(table == imethod))

        # --------------------------------------------------------------------
        # Save table with number of minimum values. Maximum values are best:
        # --------------------------------------------------------------------
        np.savetxt('number_of_min_values_8tables_3methods-mb-bao-li.csv',
                   number_of_min_values, delimiter=',')

        # # Save results as separate tables:
        # np.savetxt('mean-to-border-left_mu-sd-min-max_mb-bao-li.csv', number_of_min_values[0,:,:], delimiter=',')
        # np.savetxt('mean-from-border-left_mu-sd-min-max_mb-bao-li.csv', number_of_min_values[1,:,:], delimiter=',')
        # np.savetxt('mean-to-border-right_mu-sd-min-max_mb-bao-li.csv', number_of_min_values[2,:,:], delimiter=',')
        # np.savetxt('mean-from-border-right_mu-sd-min-max_mb-bao-li.csv', number_of_min_values[3,:,:], delimiter=',')
        # np.savetxt('sd-to-border-left_mu-sd-min-max_mb-bao-li.csv', number_of_min_values[4,:,:], delimiter=',')
        # np.savetxt('sd-from-border-left_mu-sd-min-max_mb-bao-li.csv', number_of_min_values[5,:,:], delimiter=',')
        # np.savetxt('sd-to-border-right_mu-sd-min-max_mb-bao-li.csv', number_of_min_values[6,:,:], delimiter=',')
        # np.savetxt('sd-from-border-right_mu-sd-min-max_mb-bao-li.csv', number_of_min_values[7,:,:], delimiter=',')

