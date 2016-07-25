#!/usr/bin/env python

"""
Compute surface and volume label overlaps.

Compute the Dice and Jaccard overlap measures for each labeled region
of two labeled surfaces or image volumes, for example one that has been
manually labeled and one that has been automatically labeled.

Results have been saved to https://osf.io/7gmeq/


Authors:
    - Arno Klein, 2012-2015 (arno@mindboggle.info)  http://binarybottle.com

Copyright 2015,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def evaluate_volume_overlaps(labels, file1, file2, 
                             output_file='', save_output=True):
    """
    Compute overlap between individual label regions
    in source and target nifti (nii.gz) images.

    Parameters
    ----------
    labels : list
        label indices
    file1 : string
        source image, consisting of index-labeled pixels/voxels
    file2 : string
        target image, consisting of index-labeled pixels/voxels
    output_file : string
        (optional) output file name
    save_output : bool
        save output file?

    Returns
    -------
    dice_overlaps : numpy array
        Dice overlap values
    jacc_overlaps : numpy array
        Jaccard overlap values
    output_file : string
        output text file name with overlap values

    Examples
    --------
    >>> # Compare FreeSurfer and ants labels for the same brain:
    >>> import os
    >>> from mindboggle.evaluate.evaluate_labels import evaluate_volume_overlaps
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> file1 = fetch_data(urls['freesurfer_labels'], '', '.nii.gz')
    >>> file2 = fetch_data(urls['ants_labels'], '', '.nii.gz')
    >>> dkt = DKTprotocol()
    >>> labels = dkt.cerebrum_cortex_DKT31_numbers
    >>> output_file = ''
    >>> save_output = True
    >>> evaluate_volume_overlaps(labels, file1, file2,
    ...     output_file=output_file, save_output=save_output) # doctest: +SKIP

    """
    import nibabel as nb

    from mindboggle.guts.compute import compute_overlaps

    # Load labeled image volumes:
    list1 = nb.load(file1).get_data().ravel()
    list2 = nb.load(file2).get_data().ravel()

    dice_overlaps, jacc_overlaps, output_file = compute_overlaps(labels,
        list1, list2, output_file=output_file, save_output=save_output)

    return dice_overlaps, jacc_overlaps, output_file


def evaluate_surface_overlaps(labels, index, table1, table2,
                              output_file='', save_output=True):
    """
    Measure surface overlap per label by comparing Mindboggle vertices tables.

    Parameters
    ----------
    labels : list
        cortical label indices to measure surface overlap
    index : integer
        index (starting from zero) to column of table containing label indices
    table1 : string
        table with index labels for scalar values
    table2 : string
        table with index labels for scalar values
    output_file : string
        (optional) output file name

    Returns
    -------
    dice_overlaps : numpy array
        Dice overlap values
    jacc_overlaps : numpy array
        Jaccard overlap values
    output_file : string
        (optional) output file name
    save_output : bool
        save output file?

    Examples
    --------
    >>> # Compare volume label overlaps in trivial case: brain with itself:
    >>> import os
    >>> from mindboggle.evaluate.evaluate_labels import evaluate_surface_overlaps
    >>> from mindboggle.mio.labels import DKTprotocol
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> table1 = fetch_data(urls['left_vertices_table'], '', '.csv')
    >>> table2 = fetch_data(urls['left_vertices_table'], '', '.csv')
    >>> dkt = DKTprotocol()
    >>> labels = dkt.cerebrum_cortex_DKT31_numbers
    >>> index = 1
    >>> output_file = ''
    >>> save_output = True
    >>> evaluate_surface_overlaps(labels, index, table1, table2,
    ...     output_file=output_file, save_output=save_output) # doctest: +SKIP

    """
    import pandas as pd

    from mindboggle.guts.compute import compute_overlaps

    # Load surface label tables:
    df1 = pd.read_csv(table1)
    df2 = pd.read_csv(table2)
    list1 = df1.iloc[:, index]
    list2 = df2.iloc[:, index]
    print(list1)
    dice_overlaps, jacc_overlaps, output_file = compute_overlaps(labels,
        list1, list2, output_file=output_file, save_output=save_output)

    return dice_overlaps, jacc_overlaps, output_file


def evaluate_surface_overlaps_cpp(command, labels_file1, labels_file2,
                                  output_file):
    """
    Measure surface overlap using Joachim Giard's code.

    Note: Fails if two files have different number of vertices.

    Parameters
    ----------
    command : string
        surface overlap C++ executable command
    labels_file1 : string
        ``vtk file`` with index labels for scalar values
    labels_file2 : string
        ``vtk file`` with index labels for scalar values
    output_file : string
        (optional) output file name

    Returns
    -------
    output_file : string
        name of output text file with overlap results

    Examples
    --------
    >>> # Compare surface label overlaps in trivial case: brain with itself:
    >>> import os
    >>> from mindboggle.evaluate.evaluate_labels import evaluate_surface_overlaps_cpp
    >>> from mindboggle.mio.fetch_data import fetch_data
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> file1 = fetch_data(urls['left_freesurfer_labels'], '', '.nii.gz')
    >>> file2 = fetch_data(urls['left_freesurfer_labels'], '', '.nii.gz')
    >>> ccode_path = os.environ['MINDBOGGLE_TOOLS'] # doctest: +SKIP
    >>> command = os.path.join(ccode_path, 'surface_overlap', 'SurfaceOverlapMain') # doctest: +SKIP
    >>> output_file = ''
    >>> evaluate_surface_overlaps_cpp(command, file1, file2, output_file) # doctest: +SKIP

    """
    import os
    from nipype.interfaces.base import CommandLine

    if not output_file:
        output_file = os.path.basename(labels_file1) + '_and_' + \
                           os.path.basename(labels_file2) + '.txt'
    output_file = os.path.join(os.getcwd(), output_file)
    cli = CommandLine(command = command)
    cli.inputs.args = ' '.join([labels_file1, labels_file2, output_file])
    cli.cmdline
    cli.run()

    return output_file


# ----------------------------------------------------------------------------
# Run evaluate_labels.py on Mindboggle-101 data
# to compare manual and automated volume labels and surface labels.
# ----------------------------------------------------------------------------
if __name__ == "__main__":

    import os
    import numpy as np
    import pandas as pd

    from mindboggle.mio.labels import DKTprotocol
    from mindboggle.evaluate.evaluate_labels import evaluate_volume_overlaps
    from mindboggle.evaluate.evaluate_labels import evaluate_surface_overlaps

    from math import pi
    from bokeh.models import HoverTool
    from bokeh.plotting import ColumnDataSource, figure, show, output_file
    from mindboggle.mio.colors import viridis_colormap

    measure_volume_overlaps = False
    measure_surface_overlaps = False
    plot_volume_overlaps = True

    # ------------------------------------------------------------------------
    # Plot setup:
    # ------------------------------------------------------------------------
    if plot_volume_overlaps:
        concat_tables = False  # Run first step (concatenate tables):
        dirnum = 2  # Choose directory of tables
        vol_overlap_path = '/Users/arno/Data/tables_mindboggle_vs_manual_label_overlaps'
        if dirnum == 1:
            vol_overlap_dir = os.path.join(vol_overlap_path,
                'volume_overlap_mindboggle_vs_ants_filled_labels')
            file_append = '_ants_filled_labels_overlap.csv'
            title = "Volume overlap: Manual labels vs. Mindboggle with ANTs-filled labels"
        elif dirnum == 2:
            vol_overlap_dir = os.path.join(vol_overlap_path,
                'volume_overlap_mindboggle_vs_freesurfer_filled_labels')
            file_append = '_freesurfer_wmparc_filled_labels_overlap.csv'
            title = "Volume overlap: Manual labels vs. Mindboggle with FreeSurfer-filled labels"
        elif dirnum == 3:
            vol_overlap_dir = os.path.join(vol_overlap_path,
                'volume_overlap_mindboggle_vs_freesurfer_filled_labels_no_ants')
            file_append = '_freesurfer_wmparc_filled_labels_in_freesurfer_segmentation_overlap.csv'
            title = "Volume overlap: Manual labels vs. Mindboggle with FreeSurfer-filled labels (no ANTs segmentation)"

    # ------------------------------------------------------------------------
    # Subjects:
    # ------------------------------------------------------------------------
    subject_list = '/Users/arno/Data/subject_list_Mindboggle101.txt'
    fid = open(subject_list, 'r')
    subjects = [x.strip() for x in fid.readlines()]

    # ------------------------------------------------------------------------
    # Alternating left, right cortex label numbers:
    # ------------------------------------------------------------------------
    dkt = DKTprotocol()
    labels_left = dkt.left_cerebrum_cortex_DKT31_numbers
    labels_right = dkt.right_cerebrum_cortex_DKT31_numbers
    DKT31_names = dkt.DKT31_names
    labels = []
    label_names = []
    for ilabel, label_left in enumerate(labels_left):
        labels.append(label_left)
        labels.append(labels_right[ilabel])
        label_names.append(DKT31_names[ilabel] + ' (left)')
        label_names.append(DKT31_names[ilabel] + ' (right)')
    labels = [str(x) for x in labels]

    # ------------------------------------------------------------------------
    # Colors:
    # ------------------------------------------------------------------------
    colors = viridis_colormap()
    #from matplotlib import cm as cmaps
    #import matplotlib.pyplot as plt
    #plt.register_cmap(name='viridis', cmap=cmaps.viridis)
    #plt.set_cmap(cmaps.viridis)

    # ------------------------------------------------------------------------
    # Settings:
    # ------------------------------------------------------------------------
    label_method = 'freesurfer'  # 'ants'
    use_ants_segmentation = True  # False

    # ------------------------------------------------------------------------
    # File names, paths:
    # ------------------------------------------------------------------------
    if use_ants_segmentation:
        ants_str = ''
    else:
        ants_str = '_no_ants'
    if label_method == 'ants':
        volstem = 'ants_filled_labels'
    else:
        volstem = 'freesurfer_wmparc_filled_labels'
        if not use_ants_segmentation:
            volstem += '_in_freesurfer_segmentation'
    volfile = volstem + '.nii.gz'

    names = ['OASIS-TRT-20', 'MMRR-21', 'NKI-RS-22', 'NKI-TRT-20',
             'Afterthought', 'Colin27', 'Twins-2', 'MMRR-3T7T-2', 'HLN-12']
    numbers = [20,21,22,20, 1,1,2,2,12]
    mindboggled = '/mnt/nfs-share/Mindboggle101/mindboggled/auto' + ants_str
    labels_dir = '/mnt/nfs-share/Mindboggle101/mindboggled/manual' + ants_str

    # ------------------------------------------------------------------------
    # Evaluate surface labels:
    # ------------------------------------------------------------------------
    if measure_surface_overlaps:

        surfs = ['left_cortical_surface', 'right_cortical_surface']
        index = 1
        for iname, name in enumerate(names):
            number = numbers[iname]
            for n in range(1, number+1):
                subject = name+'-'+str(n)
                for surf in surfs:
                    file1 = os.path.join(mindboggled, subject, 'tables',
                                         surf, 'vertices.csv')
                    file2 = os.path.join(labels_dir, subject, 'tables',
                                         surf, 'vertices.csv')
                    print(file1)
                    print(file2)
                    output_file = "{0}_{1}_surface_label_overlaps.csv".\
                        format(subject, surf)
                    evaluate_surface_overlaps(dkt.cerebrum_cortex_DKT31_numbers,
                                              index, file1, file2, output_file)

    # ------------------------------------------------------------------------
    # Evaluate volume labels:
    # ------------------------------------------------------------------------
    if measure_volume_overlaps:

        for iname, name in enumerate(names):
            number = numbers[iname]
            for n in range(1, number+1):
                subject = name+'-'+str(n)
                file1 = os.path.join(mindboggled, subject, 'labels', volfile)
                file2 = os.path.join(labels_dir, subject, 'labels', volfile)
                print(file1)
                print(file2)
                output_file = "{0}_{1}_volume_label_overlaps.csv".format(
                    subject, volstem)
                evaluate_volume_overlaps(dkt.label_numbers,
                                         file1, file2, output_file)

    # ------------------------------------------------------------------------
    # Plot heatmap of overlap of manual and automated volume labels:
    # ------------------------------------------------------------------------
    if plot_volume_overlaps:

        # --------------------------------------------------------------------
        # Concatenate volume overlap files:
        # --------------------------------------------------------------------
        html_file = vol_overlap_dir + '.html'
        jaccard_file = os.path.join(vol_overlap_path, 'subjects_by_labels_jaccards' + file_append)
        dice_file = os.path.join(vol_overlap_path, 'subjects_by_labels_dices' + file_append)
        vol_overlap_path = os.path.join(vol_overlap_path, vol_overlap_dir)
        files = os.listdir(vol_overlap_path)
        subjects = [x.strip(file_append) for x in files]
        if concat_tables:
            jaccards = np.zeros((len(files), len(labels)))
            dices = jaccards.copy()
            for ifile, file in enumerate(files):
                file = os.path.join(vol_overlap_path, file)
                columns = pd.read_csv(file)
                for ilabel, label in enumerate(labels):
                    for irow in range(columns.shape[0]):
                        if int(columns.iloc[irow][0].split()[0]) == label:
                            row = columns.iloc[irow][0].split()[1::]
                            dices[ifile, ilabel] = np.float(row[0])
                            jaccards[ifile, ilabel] = np.float(row[1])

            df_jaccards = pd.DataFrame(jaccards, index=subjects, columns=labels)
            df_dices = pd.DataFrame(dices, index=subjects, columns=labels)
            df_jaccards.to_csv(jaccard_file)
            df_dices.to_csv(dice_file)

        # --------------------------------------------------------------------
        # Plot heatmap for labels X subjects array:
        # --------------------------------------------------------------------
        data_file = dice_file #jaccard_file
        data = pd.read_csv(data_file)
        data = data.set_index('Unnamed: 0')

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
                valuex.append(data[label][subject])
                rgb = [np.int(255 * x) for x in
                       colors[np.int(255 * data.iloc[isubject, ilabel])]]
                hex = "#%02x%02x%02x" % tuple(rgb)
                colorx.append(hex)

        output_file(html_file, title="Volume overlap")
        source = ColumnDataSource(dict(subject=subjectx, label=labelx,
                                       label_name=label_namex,
                                       color=colorx, value=valuex))
        TOOLS = "hover,save,pan,box_zoom,wheel_zoom"

        plot_width = len(subjects) * 20
        plot_height = (len(labels)) * 20
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

        show(p)      # show the plot
