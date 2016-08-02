#!/usr/bin/env python
"""
Evaluate the consistency of Mindboggle shape measures for rescanned MRI data.

Results have been saved to https://osf.io/mhc37/

Examples
--------
$ python evaluate_shapes.py


Authors:
    - Arno Klein, 2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

# ----------------------------------------------------------------------------
# (1) Compare related surface shape measures across vertices, average over subjects.
# (2) Compare thickness measures.
# (3) Compare shapes between MRI scans.
# (4) Compare shapes between hemispheres.
# ----------------------------------------------------------------------------
if __name__ == "__main__":

    compare_shapes_between_scans = True
    compare_shapes_between_hemispheres = False
    compare_surface_shape_measures_by_vertex = False
    compare_thickness_measures = False

    if compare_surface_shape_measures_by_vertex:
        compare_surface_shape_measures_by_vertex()
    if compare_thickness_measures:
        compare_thickness_measures()
    if compare_shapes_between_scans:
        compare_shapes_between_scans()
    if compare_shapes_between_hemispheres:
        compare_shapes_between_hemispheres()

    # ------------------------------------------------------------------------
    # (1) Compare related surface shape measures across vertices,
    #     average over subjects (NOTE: SLOW!):
    # ------------------------------------------------------------------------
    def compare_surface_shape_measures_by_vertex():

        import os
        import pandas as pd
        import numpy as np

        from mindboggle.guts.compute import distcorr
        from mindboggle.mio.labels import DKTprotocol

        dkt = DKTprotocol()
        label_namesL = dkt.left_cerebrum_cortex_DKT31_names
        label_namesR = dkt.right_cerebrum_cortex_DKT31_names
        labelsL = dkt.left_cerebrum_cortex_DKT31_numbers
        labelsR = dkt.right_cerebrum_cortex_DKT31_numbers
        label_names_bilateral = dkt.DKT31_names

        subject_list = '/Users/arno/Data/subject_list_Mindboggle101.txt'
        fid = open(subject_list, 'r')
        subjects = [x.strip() for x in fid.readlines()]

        table_dir = '/Users/arno/Data/manual_tables'
        table_pathL = 'tables/left_cortical_surface/vertices.csv'
        table_pathR = 'tables/right_cortical_surface/vertices.csv'

        # --------------------------------------------------------------------
        # Loop through subjects and save distance correlations between
        # different curvature and between different depth shape measures:
        # --------------------------------------------------------------------
        dcors = np.zeros((len(subjects), len(labelsL), 4))
        for isubject, subject in enumerate(subjects):

            # Load shape tables:
            tableL = os.path.join(table_dir, subject, table_pathL)
            tableR = os.path.join(table_dir, subject, table_pathR)
            columnsL = pd.read_csv(tableL, sep=",", index_col="label ID")
            columnsR = pd.read_csv(tableR, sep=",", index_col="label ID")

            for ilabel, labelL in enumerate(labelsL):
                print(subject + ', ' + str(labelL))
                labelR = labelsR[ilabel]
                columnc1L = columnsL.loc[[labelL], ['mean curvature']].iloc[:,0].values
                columnc2L = columnsL.loc[[labelL], ['freesurfer curvature']].iloc[:,0].values
                columnd1L = columnsL.loc[[labelL], ['travel depth']].iloc[:,0].values
                columnd2L = columnsL.loc[[labelL], ['geodesic depth']].iloc[:,0].values
                columnc1R = columnsR.loc[[labelR], ['mean curvature']].iloc[:,0].values
                columnc2R = columnsR.loc[[labelR], ['freesurfer curvature']].iloc[:,0].values
                columnd1R = columnsR.loc[[labelR], ['travel depth']].iloc[:,0].values
                columnd2R = columnsR.loc[[labelR], ['geodesic depth']].iloc[:,0].values

                # Compute distance correlations:
                dcors[isubject, ilabel, 0] = distcorr(columnc1L, columnc2L)
                dcors[isubject, ilabel, 1] = distcorr(columnc1R, columnc2R)
                dcors[isubject, ilabel, 2] = distcorr(columnd1L, columnd2L)
                dcors[isubject, ilabel, 3] = distcorr(columnd1R, columnd2R)

        # --------------------------------------------------------------------
        # Save csv files:
        # --------------------------------------------------------------------
        data = pd.DataFrame(dcors[:,:,0].transpose(),
                            index=label_names_bilateral,
                            columns=[x for x in range(101)])
        data.to_csv('mean_and_FS_curvature_distance_correlation_'
                    'per_left_label_vertices_Mindboggle101.csv')

        data = pd.DataFrame(dcors[:,:,1].transpose(),
                            index=label_names_bilateral,
                            columns=[x for x in range(101)])
        data.to_csv('mean_and_FS_curvature_distance_correlation_'
                    'per_right_label_vertices_Mindboggle101.csv')

        data = pd.DataFrame(dcors[:,:,2].transpose(),
                            index=label_names_bilateral,
                            columns=[x for x in range(101)])
        data.to_csv('geodesic_and_travel_depth_distance_correlation_'
                    'per_left_label_vertices_Mindboggle101.csv')

        data = pd.DataFrame(dcors[:,:,3].transpose(),
                            index=label_names_bilateral,
                            columns=[x for x in range(101)])
        data.to_csv('geodesic_and_travel_depth_distance_correlation_'
                    'per_right_label_vertices_Mindboggle101.csv')

        data = dcors.mean(axis=0)
        data = pd.DataFrame(data, index=label_names_bilateral,
                            columns=['mean / freesurfer curvature distance correlation (left)',
                                     'mean / freesurfer curvature distance correlation (right)',
                                     'geodesic / travel depth distance correlation (left)',
                                     'geodesic / travel depth distance correlation (right)'])
        data.to_csv('mean_and_FS_curvature_geodesic_and_travel_depth_distance_correlations_'
                    'per_label_vertices_avg_over_Mindboggle101.csv')


    # ------------------------------------------------------------------------
    # (2) Compare thickness measures across subjects for each label:
    # ------------------------------------------------------------------------
    def compare_thickness_measures():

        import os
        import pandas as pd
        import numpy as np

        from mindboggle.guts.compute import distcorr
        from mindboggle.mio.labels import DKTprotocol

        dkt = DKTprotocol()
        label_names = dkt.cerebrum_cortex_DKT31_names

        subject_list = '/Users/arno/Data/subject_list_Mindboggle101.txt'
        fid = open(subject_list, 'r')
        subjects = [x.strip() for x in fid.readlines()]

        table_dir = '/Users/arno/Data/manual_tables'
        table_path1a = 'tables/left_cortical_surface/label_shapes.csv'
        table_path1b = 'tables/right_cortical_surface/label_shapes.csv'
        table_path2 = 'tables/thickinthehead_per_freesurfer_cortex_label.csv'

        # --------------------------------------------------------------------
        # Loop through subjects and table columns:
        # --------------------------------------------------------------------
        subjects_by_labels1 = np.zeros((len(subjects), len(label_names)))
        subjects_by_labels2 = np.zeros((len(subjects), len(label_names)))
        for isubject, subject in enumerate(subjects):
            # Load shape tables:
            table1a = os.path.join(table_dir, subject, table_path1a)
            table1b = os.path.join(table_dir, subject, table_path1b)
            table2 = os.path.join(table_dir, subject,  table_path2)
            columns1a = pd.read_csv(table1a, sep=",", index_col='name')
            columns1b = pd.read_csv(table1b, sep=",", index_col='name')
            column1 = columns1a['freesurfer thickness: median'] + \
                      columns1b['freesurfer thickness: median']
            column1index = column1.index
            columns2 = pd.read_csv(table2, sep=",", index_col='name')
            column2 = columns2.iloc[:, 1]
            column2_match = []
            for icolumn2, column2_index in enumerate(column2.index):
                if column2_index in column1index:
                    column2_match.append(column2[icolumn2])
            subjects_by_labels1[isubject, :] = column1
            subjects_by_labels2[isubject, :] = column2_match

        dcors = []
        for ilabel in range(len(label_names)):
            dcors.append(distcorr(subjects_by_labels1[:, ilabel],
                                  subjects_by_labels2[:, ilabel]))

        # --------------------------------------------------------------------
        # Save csv files:
        # --------------------------------------------------------------------
        data = pd.DataFrame(dcors, index=label_names, #index=columns1.columns)
                            columns=['freesurfer / thickinthehead cortical thickness distance correlation'])
        data.to_csv('thickinthehead_FSthickness_distance_correlations_'
                    'per_label_Mindboggle101.csv')


    # ------------------------------------------------------------------------
    # (3) Compare shapes between MRI scans.
    #
    # To evaluate Mindboggle shape measure consistency, run Mindboggle on 41
    # Mindboggle-101 subjects for which there is a second MRI scan
    # (OASIS-TRT-20 and MMRR-21 groups). Since cortical surface reconstruction
    # results in a different number of vertices for the two scans, we compared
    # median values for each cortical region across the two scans.
    #
    # The shape differences are computed for each of the 62 (volume) or 31
    # (left surface) cortical regions as the difference between the region's
    # shape values between the two scans divided by the first scan's shape
    # value. For the surface-based shape values, we used the median value for
    # all vertices within each region.
    #
    # Note: some right hemisphere surfaces had NaN values, so we analyzed
    # whole volumes and left hemisphere surfaces.
    # ------------------------------------------------------------------------
    def compare_shapes_between_scans():
        import os
        import numpy as np
        import pandas as pd

        # For plotting:
        from math import pi
        from bokeh.models import HoverTool
        from bokeh.plotting import ColumnDataSource, figure, show, save, output_file
        from mindboggle.mio.colors import viridis_colormap
        from mindboggle.mio.labels import DKTprotocol
        #from mindboggle.mio.plots import histograms_of_lists

        titles = ["Fractional difference between re/scan volumes",
                  "Fractional difference between re/scan thickinthehead cortical thicknesses",
                  "Fractional difference between re/scan left cortical label median areas",
                  "Fractional difference between re/scan left cortical label median travel depths",
                  "Fractional difference between re/scan left cortical label median geodesic depths",
                  "Fractional difference between re/scan left cortical label median mean curvatures",
                  "Fractional difference between re/scan left cortical label median FreeSurfer curvatures",
                  "Fractional difference between re/scan left cortical label median FreeSurfer thicknesses"]
                  # "Fractional difference between re/scan right cortical label median areas",
                  # "Fractional difference between re/scan right cortical label median travel depths",
                  # "Fractional difference between re/scan right cortical label median geodesic depths",
                  # "Fractional difference between re/scan right cortical label median mean curvatures",
                  # "Fractional difference between re/scan right cortical label median FreeSurfer curvatures",
                  # "Fractional difference between re/scan right cortical label median FreeSurfer thicknesses"]
                  # "Fractional difference between re/scan right cortical label median FreeSurfer convexities"]
        names = ["volume_for_each_freesurfer_label",
                 "thickinthehead_per_freesurfer_cortex_label",
                 "median_area_per_freesurfer_left_cortex_label",
                 "median_travel_depth_per_freesurfer_left_cortex_label",
                 "median_geodesic_depth_per_freesurfer_left_cortex_label",
                 "median_mean_curvatures_per_freesurfer_left_cortex_label",
                 "median_freesurfer_curvature_per_freesurfer_left_cortex_label",
                 "median_freesurfer_thickness_per_freesurfer_left_cortex_label"]
                 # "median_area_per_freesurfer_right_cortex_label",
                 # "median_travel_depth_per_freesurfer_right_cortex_label",
                 # "median_geodesic_depth_per_freesurfer_right_cortex_label",
                 # "median_mean_curvatures_per_freesurfer_right_cortex_label",
                 # "median_freesurfer_curvature_per_freesurfer_right_cortex_label",
                 # "median_freesurfer_thickness_per_freesurfer_right_cortex_label"]
                 # "median_freesurfer_convexity_per_freesurfer_right_cortex_label"]
        table_dir = '/Users/arno/Data/shape_tables_for_auto_labels_of_Mindboggle101_rescans'
        tables = [os.path.join('tables', 'volume_for_each_freesurfer_label.csv'),
                  os.path.join('tables', 'thickinthehead_per_freesurfer_cortex_label.csv'),
                  os.path.join('tables', 'left_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'left_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'left_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'left_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'left_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'left_cortical_surface', 'label_shapes.csv')]
                  # os.path.join('tables', 'right_cortical_surface', 'label_shapes.csv'),
                  # os.path.join('tables', 'right_cortical_surface', 'label_shapes.csv'),
                  # os.path.join('tables', 'right_cortical_surface', 'label_shapes.csv'),
                  # os.path.join('tables', 'right_cortical_surface', 'label_shapes.csv'),
                  # os.path.join('tables', 'right_cortical_surface', 'label_shapes.csv'),
                  # os.path.join('tables', 'right_cortical_surface', 'label_shapes.csv')]
        column_indices = [1, 1, 1, 2, 10, 18, 26, 34] #, 1, 2, 10, 18, 26, 34] #, 42]

        # --------------------------------------------------------------------
        # Alternating left, right cortex label numbers (for volume shapes):
        # --------------------------------------------------------------------
        dkt = DKTprotocol()
        labels_left = dkt.left_cerebrum_cortex_DKT31_numbers
        labels_right = dkt.right_cerebrum_cortex_DKT31_numbers
        DKT31_names = dkt.DKT31_names
        # label_list = []
        # label_name_list = []
        # for ilabel, label_left in enumerate(labels_left):
        #     label_list.append(label_left)
        #     label_list.append(labels_right[ilabel])
        #     label_name_list.append(DKT31_names[ilabel] + ' (left)')
        #     label_name_list.append(DKT31_names[ilabel] + ' (right)')
        #label_list = [str(x) for x in label_list]
        ##exclude_sulci = [20] # Sulcus 20 removed from protocol since initial run

        #label_lists = [label_list,
        #               label_list,
        label_lists = [labels_left,
                       labels_left,
                       labels_left, labels_left, labels_left,
                       labels_left, labels_left, labels_left, labels_left,
                       labels_right, labels_right, labels_right,
                       labels_right, labels_right, labels_right, labels_right]
        label_name_lists = [DKT31_names for x in range (len(titles))]
        #label_name_lists[0] = label_name_list
        #label_name_lists[1] = label_name_list

        # --------------------------------------------------------------------
        # Colors:
        # --------------------------------------------------------------------
        colors = viridis_colormap()
        #from matplotlib import cm as cmaps
        #import matplotlib.pyplot as plt
        #plt.register_cmap(name='viridis', cmap=cmaps.viridis)
        #plt.set_cmap(cmaps.viridis)

        scale_rect = 40

        # --------------------------------------------------------------------
        # Subjects with second scans:
        # --------------------------------------------------------------------
        groups = ['OASIS-TRT-20', 'MMRR-21']
        numbers = [20, 21]
        nsubjects = sum(numbers)
        subjects = []
        subjects2 = []
        for igroup, group in enumerate(groups):
            for n in range(1, numbers[igroup]+1):
                subjects.append(group+'-'+str(n))
                subjects2.append(group+'-rescan-'+str(n))

        # --------------------------------------------------------------------
        # Loop through tables:
        # --------------------------------------------------------------------
        data_means = np.zeros((len(labels_left), len(titles)))
        data_summaries = np.zeros((len(titles), 8))
        for ititle, title in enumerate(titles):
            table = tables[ititle]
            name = names[ititle]
            index = column_indices[ititle]
            labels = label_lists[ititle]
            label_names = label_name_lists[ititle]

            # ----------------------------------------------------------------
            # Loop through subjects:
            # ----------------------------------------------------------------
            subject_shapes = np.zeros((nsubjects, len(labels)))
            subject2_shapes = np.zeros((nsubjects, len(labels)))
            for isubject, subject in enumerate(subjects):
                subject2 = subjects2[isubject]
                table_file = os.path.join(table_dir, subject, table)
                table_file2 = os.path.join(table_dir, subject2, table)
                columns = pd.read_csv(table_file, sep=",", index_col='name')
                columns2 = pd.read_csv(table_file2, sep=",", index_col='name')

                # ------------------------------------------------------------
                # Loop through labels:
                # ------------------------------------------------------------
                for ilabel, label in enumerate(labels):
                    for irow in range(columns.shape[0]):
                        if int(columns.iloc[irow][0]) == int(label):
                            value = columns.iloc[irow][index]
                            value2 = columns2.iloc[irow][index]
                            subject_shapes[isubject, ilabel] = value
                            subject2_shapes[isubject, ilabel] = value2

            # ----------------------------------------------------------------
            # Save csv files:
            # ----------------------------------------------------------------
            data = pd.DataFrame(subject_shapes, index=subjects, columns=labels)
            data.to_csv(name + '_scans.csv')
            data_summary = data.describe(include='all')
            data_summary.to_csv(name + '_scans_summary.csv')

            data = pd.DataFrame(subject2_shapes, index=subjects, columns=labels)
            data.to_csv(name + '_rescans.csv')
            data_summary = data.describe(include='all')
            data_summary.to_csv(name + '_rescans_summary.csv')

            subject_shape_diffs = subject2_shapes - subject_shapes
            data = pd.DataFrame(subject_shape_diffs, index=subjects, columns=labels)
            data.to_csv(name + '_differences.csv')
            data_summary = data.describe(include='all')
            data_summary.to_csv(name + '_differences_summary.csv')

            subject_shape_abs_diffs = np.abs(subject_shape_diffs)
            max_diffs = subject_shape_abs_diffs.max(axis=0)
            subject_shape_frac_diffs = subject_shape_diffs / subject_shapes
            data = pd.DataFrame(subject_shape_frac_diffs,
                                index=subjects, columns=labels)
            #iInf, jInf = np.where(data.values == np.inf)
            #data.iloc[iInf, jInf] = 'NaN'
            data.to_csv(name + '_fractional_differences.csv')
            data_summary = data.describe(include='all')
            data_summary.to_csv(name + '_fractional_differences_summary.csv')

            # max_array = np.zeros((nsubjects, len(labels), 2))
            # max_array[:, :, 0] = subject_shapes
            # max_array[:, :, 1] = subject2_shapes
            # max_diff = np.abs(np.max(max_array, axis=2))

            # maxs = np.max(np.max([subject_shapes, subject2_shapes], axis=0), axis=0)
            # mins = np.min(np.min([subject_shapes, subject2_shapes], axis=0), axis=0)
            # max_diff = np.abs(maxs - mins)

            # max_diff = np.max([subject_shapes, subject2_shapes]) -\
            #           np.min([subject_shapes, subject2_shapes])

            subject_shape_frac_abs_diffs = np.abs(subject_shape_abs_diffs / subject_shapes)
            data = pd.DataFrame(subject_shape_frac_abs_diffs,
                                index=subjects, columns=labels)
            n50 = len(np.where(data.values > 0.5)[0])
            n25 = len(np.where(data.values > 0.25)[0])
            n10 = len(np.where(data.values > 0.1)[0])
            print(title)
            print("Fractional absolute differences above "
                  "0.5: {0}; 0.25: {1}; 0.1: {2}".format(n50, n25, n10))
            print("")
            #iInf, jInf = np.where(data.values == np.inf)
            #data.iloc[iInf, jInf] = 'NaN'
            data.to_csv(name + '_fractional_abs_differences.csv')
            data_summary = data.describe(include='all')
            data_summary.to_csv(name + '_fractional_abs_differences_summary.csv')

            data_means[:, ititle] = data_summary.loc['mean'].values
            data_summaries[ititle, :] = data_summary.mean(axis=1)

            # ignore_columns = []
            # nbins = 100
            # axis_limits = []
            # histograms_of_lists(subject_shape_diffs, title, ignore_columns,
            #                     nbins, axis_limits, [title])

            # ----------------------------------------------------------------
            # Plot heatmap for labels X subjects array for each table:
            # ----------------------------------------------------------------
            html_file = name + '_fractional_abs_differences.html'
            print(html_file)

            # Set up the data for plotting. We will need to have values for every
            # pair of subject/label names. Map the value to a color.
            subjectx = []
            labelx = []
            label_namex = []
            value1x = []
            value2x = []
            differencex = []
            fractionx = []
            colorx = []
            for ilabel, label in enumerate(labels):
                for isubject, subject in enumerate(subjects):
                    labelx.append(label)
                    label_namex.append(label_names[ilabel])
                    subjectx.append(subject)
                    value1 = subject_shapes[isubject, ilabel]
                    value2 = subject2_shapes[isubject, ilabel]
                    difference = subject_shape_diffs[isubject, ilabel]
                    value1x.append(value1)
                    value2x.append(value2)
                    differencex.append(difference)
                    fraction = subject_shape_frac_diffs[isubject, ilabel]
                    fractionx.append(fraction)
                    abs_fraction = subject_shape_frac_abs_diffs[isubject, ilabel]
                    if np.isnan(value) or np.isnan(abs_fraction):
                        rgb = [0, 0, 0]
                    elif abs_fraction > 1.0:
                        rgb = [1, 1, 1]
                    else:
                        rgb = [np.int(255 * x) for x in
                               colors[np.int(255 * abs_fraction)]]
                    hex = "#%02x%02x%02x" % tuple(rgb)
                    colorx.append(hex)

            output_file(html_file, title=title)
            source = ColumnDataSource(dict(subject=subjectx,
                                           label=labelx,
                                           label_name=label_namex,
                                           color=colorx,
                                           value1=value1x,
                                           value2=value2x,
                                           difference=differencex,
                                           fraction=fractionx))
            TOOLS = "hover,save,pan,box_zoom,wheel_zoom"

            plot_width = len(subjects) * scale_rect
            plot_height = len(labels) * scale_rect
            p = figure(title=title, x_range=subjects, y_range=list(reversed(label_names)),
                       plot_width=plot_width, plot_height=plot_height,
                       x_axis_location="above", tools=TOOLS)
            p.grid.grid_line_color = None
            p.axis.axis_line_color = None
            p.axis.major_tick_line_color = None
            p.axis.major_label_text_font_size = "10pt"
            p.axis.major_label_standoff = 0
            p.xaxis.major_label_orientation = pi/3
            p.rect(x="subject", y="label_name", width=1, height=1,
                   source=source, color="color", line_color=None)

            p.select_one(HoverTool).tooltips = [
                ('subject', '@subject'),
                ('label', '@label'),
                ('label name', '@label_name'),
                ('value1', '@value1'),
                ('value2', '@value2'),
                ('difference', '@difference'),
                ('fraction', '@fraction'),
            ]

            #show(p)      # show the plot
            #import sys; sys.exit()
            save(p)      # save the plot

        data_means_df = pd.DataFrame(data_means,
                                     index=label_names,
                                     columns=names)
        data_means_df.to_csv('means_of_rescan_fractional_abs_shape_differences.csv')

        data_summaries_df = pd.DataFrame(data_summaries,
                                         index=names,
                                         columns=data_summary.index)
        data_summaries_df.to_csv('summary_of_rescan_fractional_abs_shape_differences.csv')


    # ------------------------------------------------------------------------
    # (4) Compare shapes between hemispheres for each label:
    #
    # The shape differences are computed for each of the 31 cortical regions
    # as the difference between the region's shape values between the two
    # hemispheres divided by the first scan's shape value.
    # For the surface-based shape values, we used the median value for
    # all vertices within each region.
    # ------------------------------------------------------------------------
    def compare_shapes_between_hemispheres():

        import os
        import numpy as np
        import pandas as pd

        # For plotting:
        from math import pi
        from bokeh.models import HoverTool
        from bokeh.plotting import ColumnDataSource, figure, show, save, output_file
        from mindboggle.mio.colors import viridis_colormap
        from mindboggle.mio.labels import DKTprotocol

        titles = ["Fractional difference between interhemispheric volumes",
                  "Fractional difference between interhemispheric thickinthehead cortical thicknesses",
                  "Fractional difference between interhemispheric cortical label median areas",
                  "Fractional difference between interhemispheric cortical label median travel depths",
                  "Fractional difference between interhemispheric cortical label median geodesic depths",
                  "Fractional difference between interhemispheric cortical label median mean curvatures",
                  "Fractional difference between interhemispheric cortical label median FreeSurfer curvatures",
                  "Fractional difference between interhemispheric cortical label median FreeSurfer thicknesses"]
        names = ["volume_for_each_freesurfer_label",
                 "thickinthehead_per_freesurfer_cortex_label",
                 "median_area_per_freesurfer_cortex_label",
                 "median_travel_depth_per_freesurfer_cortex_label",
                 "median_geodesic_depth_per_freesurfer_cortex_label",
                 "median_mean_curvatures_per_freesurfer_cortex_label",
                 "median_freesurfer_curvature_per_freesurfer_cortex_label",
                 "median_freesurfer_thickness_per_freesurfer_cortex_label"]
        table_dir = '/Users/arno/Data/manual_tables'
        tablesL = [os.path.join('tables', 'volume_for_each_freesurfer_label.csv'),
                  os.path.join('tables', 'thickinthehead_per_freesurfer_cortex_label.csv'),
                  os.path.join('tables', 'left_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'left_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'left_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'left_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'left_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'left_cortical_surface', 'label_shapes.csv')]
        tablesR = [os.path.join('tables', 'volume_for_each_freesurfer_label.csv'),
                  os.path.join('tables', 'thickinthehead_per_freesurfer_cortex_label.csv'),
                  os.path.join('tables', 'right_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'right_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'right_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'right_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'right_cortical_surface', 'label_shapes.csv'),
                  os.path.join('tables', 'right_cortical_surface', 'label_shapes.csv')]
        column_indices = [1, 1, 1, 2, 10, 18, 26, 34]

        # --------------------------------------------------------------------
        # Alternating left, right cortex label numbers (for volume shapes):
        # --------------------------------------------------------------------
        dkt = DKTprotocol()
        labels_left = dkt.left_cerebrum_cortex_DKT31_numbers
        labels_right = dkt.right_cerebrum_cortex_DKT31_numbers
        label_names = dkt.DKT31_names
        #exclude_sulci = [20] # Sulcus 20 removed from protocol since initial run

        # --------------------------------------------------------------------
        # Colors:
        # --------------------------------------------------------------------
        colors = viridis_colormap()
        #from matplotlib import cm as cmaps
        #import matplotlib.pyplot as plt
        #plt.register_cmap(name='viridis', cmap=cmaps.viridis)
        #plt.set_cmap(cmaps.viridis)

        scale_rect = 20

        # --------------------------------------------------------------------
        # Subjects:
        # --------------------------------------------------------------------
        subject_list = '/Users/arno/Data/subject_list_Mindboggle101.txt'
        fid = open(subject_list, 'r')
        subjects = [x.strip() for x in fid.readlines()]

        # --------------------------------------------------------------------
        # Loop through tables:
        # --------------------------------------------------------------------
        data_means = np.zeros((len(labels_left), len(titles)))
        data_summaries = np.zeros((len(titles), 8))
        for ititle, title in enumerate(titles):
            tableL_file = tablesL[ititle]
            tableR_file = tablesR[ititle]
            name = names[ititle]
            index = column_indices[ititle]

            # ----------------------------------------------------------------
            # Loop through subjects:
            # ----------------------------------------------------------------
            subject_shapesL = np.zeros((len(subjects), len(labels_left)))
            subject_shapesR = np.zeros((len(subjects), len(labels_right)))
            for isubject, subject in enumerate(subjects):
                tableL = os.path.join(table_dir, subject, tableL_file)
                tableR = os.path.join(table_dir, subject, tableR_file)
                columnsL = pd.read_csv(tableL, sep=",", index_col='name')
                columnsR = pd.read_csv(tableR, sep=",", index_col='name')

                # ------------------------------------------------------------
                # Loop through labels:
                # ------------------------------------------------------------
                for ilabel, labelL in enumerate(labels_left):
                    for irow in range(columnsL.shape[0]):
                        if int(columnsL.iloc[irow][0]) == int(labelL):
                            valueL = columnsL.iloc[irow][index]
                            subject_shapesL[isubject, ilabel] = valueL
                for ilabel, labelR in enumerate(labels_right):
                    for irow in range(columnsR.shape[0]):
                        if int(columnsR.iloc[irow][0]) == int(labelR):
                            valueR = columnsR.iloc[irow][index]
                            subject_shapesR[isubject, ilabel] = valueR

            # ----------------------------------------------------------------
            # Save csv files:
            # ----------------------------------------------------------------
            data = pd.DataFrame(subject_shapesL, index=subjects, columns=labels_left)
            data.to_csv(name + '_left.csv')
            data_summary = data.describe(include='all')
            data_summary.to_csv(name + '_left_summary.csv')

            data = pd.DataFrame(subject_shapesR, index=subjects, columns=labels_right)
            data.to_csv(name + '_right.csv')
            data_summary = data.describe(include='all')
            data_summary.to_csv(name + '_right_summary.csv')

            subject_shape_diffs = subject_shapesL - subject_shapesR
            data = pd.DataFrame(subject_shape_diffs, index=subjects, columns=label_names)
            data.to_csv(name + '_differences.csv')
            data_summary = data.describe(include='all')
            data_summary.to_csv(name + '_differences_summary.csv')

            subject_shape_abs_diffs = np.abs(subject_shape_diffs)
            subject_shape_frac_diffs = subject_shape_diffs / subject_shapesL
            data = pd.DataFrame(subject_shape_frac_diffs,
                                index=subjects, columns=label_names)
            data.to_csv(name + '_fractional_differences.csv')
            data_summary = data.describe(include='all')
            data_summary.to_csv(name + '_fractional_differences_summary.csv')

            subject_shape_frac_abs_diffs = np.abs(subject_shape_abs_diffs / subject_shapesL)
            data = pd.DataFrame(subject_shape_frac_abs_diffs,
                                index=subjects, columns=label_names)
            n50 = len(np.where(data.values > 0.5)[0])
            n25 = len(np.where(data.values > 0.25)[0])
            n10 = len(np.where(data.values > 0.1)[0])
            print(title)
            print("Fractional absolute differences above "
                  "0.5: {0}; 0.25: {1}; 0.1: {2}".format(n50, n25, n10))
            print("")
            data.to_csv(name + '_fractional_abs_differences.csv')
            data_summary = data.describe(include='all')
            data_summary.to_csv(name + '_fractional_abs_differences_summary.csv')

            data_means[:, ititle] = data_summary.loc['mean'].values
            data_summaries[ititle, :] = data_summary.mean(axis=1)

            # ----------------------------------------------------------------
            # Plot heatmap for labels X subjects array for each table:
            # ----------------------------------------------------------------
            html_file = name + '_fractional_abs_differences.html'
            print(html_file)

            # Set up the data for plotting. We will need to have values for every
            # pair of subject/label names. Map the value to a color.
            subjectx = []
            label_namex = []
            value1x = []
            value2x = []
            differencex = []
            fractionx = []
            colorx = []
            for ilabel in range(len(labels_left)):
                for isubject, subject in enumerate(subjects):
                    label_namex.append(label_names[ilabel])
                    subjectx.append(subject)
                    value1 = subject_shapesL[isubject, ilabel]
                    value2 = subject_shapesR[isubject, ilabel]
                    difference = subject_shape_diffs[isubject, ilabel]
                    value1x.append(value1)
                    value2x.append(value2)
                    differencex.append(difference)
                    fraction = subject_shape_frac_diffs[isubject, ilabel]
                    fractionx.append(fraction)
                    abs_fraction = subject_shape_frac_abs_diffs[isubject, ilabel]
                    if np.isnan(value1) or np.isnan(value2) or np.isnan(abs_fraction):
                        rgb = [0, 0, 0]
                    elif abs_fraction > 1.0:
                        rgb = [1, 1, 1]
                    else:
                        rgb = [np.int(255 * x) for x in
                               colors[np.int(255 * abs_fraction)]]
                    hex = "#%02x%02x%02x" % tuple(rgb)
                    colorx.append(hex)

            output_file(html_file, title=title)
            source = ColumnDataSource(dict(subject=subjectx,
                                           label_name=label_namex,
                                           color=colorx,
                                           value1=value1x,
                                           value2=value2x,
                                           difference=differencex,
                                           fraction=fractionx))
            TOOLS = "hover,save,pan,box_zoom,wheel_zoom"

            plot_width = len(subjects) * scale_rect
            plot_height = len(labels_left) * scale_rect
            p = figure(title=title, x_range=subjects, y_range=list(reversed(label_names)),
                       plot_width=plot_width, plot_height=plot_height,
                       x_axis_location="above", tools=TOOLS)
            p.grid.grid_line_color = None
            p.axis.axis_line_color = None
            p.axis.major_tick_line_color = None
            p.axis.major_label_text_font_size = "10pt"
            p.axis.major_label_standoff = 0
            p.xaxis.major_label_orientation = pi/3
            p.rect(x="subject", y="label_name", width=1, height=1,
                   source=source, color="color", line_color=None)

            p.select_one(HoverTool).tooltips = [
                ('subject', '@subject'),
                ('label', '@label_name'),
                ('value1', '@value1'),
                ('value2', '@value2'),
                ('difference', '@difference'),
                ('fraction', '@fraction'),
            ]

            #show(p)      # show the plot
            #import sys; sys.exit()
            save(p)      # save the plot

        data_means_df = pd.DataFrame(data_means,
                                     index=label_names,
                                     columns=names)
        data_means_df.to_csv('means_of_interhemispheric_fractional_abs_shape_differences.csv')

        data_summaries_df = pd.DataFrame(data_summaries,
                                         index=names,
                                         columns=data_summary.index)
        data_summaries_df.to_csv('summary_of_interhemispheric_fractional_abs_shape_differences.csv')


    # def compare_shapes_between_hemispheres():
    #
    #     import os
    #     import pandas as pd
    #     import numpy as np
    #
    #     from mindboggle.guts.compute import distcorr
    #     from mindboggle.mio.labels import DKTprotocol
    #
    #     dkt = DKTprotocol()
    #     label_names = dkt.cerebrum_cortex_DKT31_names
    #     label_names_bilateral = dkt.DKT31_names
    #
    #     subject_list = '/Users/arno/Data/subject_list_Mindboggle101.txt'
    #     fid = open(subject_list, 'r')
    #     subjects = [x.strip() for x in fid.readlines()]
    #
    #     table_dir = '/Users/arno/Data/manual_tables'
    #     table_pathL = 'tables/left_cortical_surface/label_shapes.csv'
    #     table_pathR = 'tables/right_cortical_surface/label_shapes.csv'
    #     table_path2 = 'tables/thickinthehead_per_freesurfer_cortex_label.csv'
    #
    #     # --------------------------------------------------------------------
    #     # Loop through subjects and save distance correlations
    #     # between hemispheres:
    #     # --------------------------------------------------------------------
    #     Z = np.zeros((len(subjects), len(label_names_bilateral), 2))
    #     areas = Z.copy()
    #     thickintheheads = Z.copy()
    #     fs_thicknesses = Z.copy()
    #     mean_curvatures = Z.copy()
    #     fs_curvatures = Z.copy()
    #     travel_depths = Z.copy()
    #     geodesic_depths = Z.copy()
    #     for isubject, subject in enumerate(subjects):
    #         # Load shape tables:
    #         tableL = os.path.join(table_dir, subject, table_pathL)
    #         tableR = os.path.join(table_dir, subject, table_pathR)
    #         table2 = os.path.join(table_dir, subject, table_path2)
    #         columnsL = pd.read_csv(tableL, sep=",", index_col='name')
    #         columnsR = pd.read_csv(tableR, sep=",", index_col='name')
    #         columns2 = pd.read_csv(table2, sep=",", index_col='name')
    #         columnsL = columnsL.iloc[0:len(label_names_bilateral), :]
    #         columnsR = columnsR.iloc[len(label_names_bilateral):len(label_names), :]
    #
    #         # area:
    #         columnL = columnsL['area']
    #         columnR = columnsR['area']
    #         areas[isubject, :, 0] = columnL
    #         areas[isubject, :, 1] = columnR
    #         columnLindex = columnL.index
    #         columnRindex = columnR.index
    #
    #         # thickinthehead:
    #         column2 = columns2.iloc[:, 1]
    #         columnt2L = []
    #         columnt2R = []
    #         for icolumn2, column2_index in enumerate(column2.index):
    #             if column2_index in columnLindex:
    #                 columnt2L.append(column2[icolumn2])
    #             elif column2_index in columnRindex:
    #                 columnt2R.append(column2[icolumn2])
    #         thickintheheads[isubject, :, 0] = columnt2L
    #         thickintheheads[isubject, :, 1] = columnt2R
    #
    #         # FreeSurfer thickness:
    #         columnL = columnsL['freesurfer thickness: median']
    #         columnR = columnsR['freesurfer thickness: median']
    #         fs_thicknesses[isubject, :, 0] = columnL
    #         fs_thicknesses[isubject, :, 1] = columnR
    #
    #         # mean curvature:
    #         columnL = columnsL['mean curvature: median']
    #         columnR = columnsR['mean curvature: median']
    #         mean_curvatures[isubject, :, 0] = columnL
    #         mean_curvatures[isubject, :, 1] = columnR
    #
    #         # FreeSurfer curvature:
    #         columnL = columnsL['freesurfer curvature: median']
    #         columnR = columnsR['freesurfer curvature: median']
    #         fs_curvatures[isubject, :, 0] = columnL
    #         fs_curvatures[isubject, :, 1] = columnR
    #
    #         # travel depth:
    #         columnL = columnsL['travel depth: median']
    #         columnR = columnsR['travel depth: median']
    #         travel_depths[isubject, :, 0] = columnL
    #         travel_depths[isubject, :, 1] = columnR
    #
    #         # geodesic depth:
    #         columnL = columnsL['geodesic depth: median']
    #         columnR = columnsR['geodesic depth: median']
    #         geodesic_depths[isubject, :, 0] = columnL
    #         geodesic_depths[isubject, :, 1] = columnR
    #
    #     dcors = np.zeros((len(label_names_bilateral), 7))
    #     for ilabel in range(len(label_names_bilateral)):
    #         dcors[ilabel, 0] = distcorr(areas[:, ilabel, 0],
    #                                     areas[:, ilabel, 1])
    #         dcors[ilabel, 1] = distcorr(thickintheheads[:, ilabel, 0],
    #                                     thickintheheads[:, ilabel, 1])
    #         dcors[ilabel, 2] = distcorr(fs_thicknesses[:, ilabel, 0],
    #                                     fs_thicknesses[:, ilabel, 1])
    #         dcors[ilabel, 3] = distcorr(mean_curvatures[:, ilabel, 0],
    #                                     mean_curvatures[:, ilabel, 1])
    #         dcors[ilabel, 4] = distcorr(fs_curvatures[:, ilabel, 0],
    #                                     fs_curvatures[:, ilabel, 1])
    #         dcors[ilabel, 5] = distcorr(travel_depths[:, ilabel, 0],
    #                                     travel_depths[:, ilabel, 1])
    #         dcors[ilabel, 6] = distcorr(geodesic_depths[:, ilabel, 0],
    #                                     geodesic_depths[:, ilabel, 1])
    #
    #     # --------------------------------------------------------------------
    #     # Save csv files:
    #     # --------------------------------------------------------------------
    #     data = pd.DataFrame(dcors, index=label_names_bilateral, #index=columns1.columns)
    #                         columns=['area', 'thickinthehead',
    #                                  'freesurfer thickness', 'mean curvature',
    #                                  'freesurfer curvature', 'travel depth',
    #                                  'geodesic depth'])
    #     data.to_csv('distance_correlations_for_shapes_between_hemispheres_'
    #                 'per_label_Mindboggle101.csv')

