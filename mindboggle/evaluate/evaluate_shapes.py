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
# Evaluate consistency of Mindboggle shape measures for rescanned MRI data.
# To evaluate Mindboggle shape measure consistency, run Mindboggle on 41
# Mindboggle-101 subjects for which there is a second MRI scan (OASIS-TRT-20
# and MMRR-21 groups). Since cortical surface reconstruction results in a
# different number of vertices for the two scans, we compared median values
# for each cortical region across the two scans.
#
# Note: some right hemisphere surfaces had NaN values, so we analyzed whole
# volumes and left hemisphere surfaces.
# ----------------------------------------------------------------------------
if __name__ == "__main__":

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

    # ------------------------------------------------------------------------
    # Alternating left, right cortex label numbers (for volume shapes):
    # ------------------------------------------------------------------------
    dkt = DKTprotocol()
    labels_left = dkt.left_cerebrum_cortex_DKT31_numbers
    labels_right = dkt.right_cerebrum_cortex_DKT31_numbers
    DKT31_names = dkt.DKT31_names
    label_list = []
    label_name_list = []
    for ilabel, label_left in enumerate(labels_left):
        label_list.append(label_left)
        label_list.append(labels_right[ilabel])
        label_name_list.append(DKT31_names[ilabel] + ' (left)')
        label_name_list.append(DKT31_names[ilabel] + ' (right)')
    label_list = [str(x) for x in label_list]
    #exclude_sulci = [20] # Sulcus 20 removed from protocol since initial run

    label_lists = [label_list,
                   label_list,
                   labels_left, labels_left, labels_left,
                   labels_left, labels_left, labels_left, labels_left,
                   labels_right, labels_right, labels_right,
                   labels_right, labels_right, labels_right, labels_right]
    label_name_lists = [DKT31_names for x in range (len(titles))]
    label_name_lists[0] = label_name_list
    label_name_lists[1] = label_name_list

    # ------------------------------------------------------------------------
    # Colors:
    # ------------------------------------------------------------------------
    colors = viridis_colormap()
    #from matplotlib import cm as cmaps
    #import matplotlib.pyplot as plt
    #plt.register_cmap(name='viridis', cmap=cmaps.viridis)
    #plt.set_cmap(cmaps.viridis)

    scale_rect = 40

    # ------------------------------------------------------------------------
    # Subjects with second scans:
    # ------------------------------------------------------------------------
    groups = ['OASIS-TRT-20', 'MMRR-21']
    numbers = [20, 21]
    nsubjects = sum(numbers)
    subjects = []
    subjects2 = []
    for igroup, group in enumerate(groups):
        for n in range(1, numbers[igroup]+1):
            subjects.append(group+'-'+str(n))
            subjects2.append(group+'-rescan-'+str(n))

    # ------------------------------------------------------------------------
    # Loop through tables:
    # ------------------------------------------------------------------------
    data_summaries = np.zeros((len(titles), 8))
    for ititle, title in enumerate(titles):
        table = tables[ititle]
        name = names[ititle]
        index = column_indices[ititle]
        labels = label_lists[ititle]
        label_names = label_name_lists[ititle]

        # --------------------------------------------------------------------
        # Loop through subjects:
        # --------------------------------------------------------------------
        subject_shapes = np.zeros((nsubjects, len(labels)))
        subject2_shapes = np.zeros((nsubjects, len(labels)))
        for isubject, subject in enumerate(subjects):
            subject2 = subjects2[isubject]
            table_file = os.path.join(table_dir, subject, table)
            table_file2 = os.path.join(table_dir, subject2, table)
            columns = pd.read_csv(table_file, sep=",", index_col='name')
            columns2 = pd.read_csv(table_file2, sep=",", index_col='name')

            # ----------------------------------------------------------------
            # Loop through labels:
            # ----------------------------------------------------------------
            for ilabel, label in enumerate(labels):
                for irow in range(columns.shape[0]):
                    if int(columns.iloc[irow][0]) == int(label):
                        value = columns.iloc[irow][index]
                        value2 = columns2.iloc[irow][index]
                        subject_shapes[isubject, ilabel] = value
                        subject2_shapes[isubject, ilabel] = value2

        # --------------------------------------------------------------------
        # Save csv files:
        # --------------------------------------------------------------------
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

        data_summaries[ititle, :] = data_summary.max(axis=1)

        # ignore_columns = []
        # nbins = 100
        # axis_limits = []
        # histograms_of_lists(subject_shape_diffs, title, ignore_columns,
        #                     nbins, axis_limits, [title])

        # --------------------------------------------------------------------
        # Plot heatmap for labels X subjects array for each table:
        # --------------------------------------------------------------------
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

    data_summaries_df = pd.DataFrame(data_summaries,
                                     index=names,
                                     columns=data_summary.index)
    data_summaries_df.to_csv('summary_of_rescan_fractional_abs_shape_differences.csv')

