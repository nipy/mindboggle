#!/usr/bin/env python
"""
Plotting functions.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Contributors:
    - Hal Canary <http://cs.unc.edu/~hal>: vtkviewer.py called by vtkviewer()

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


# # Example Pysurfer plot of FreeSurfer surface + VTK scalars
# ipython
# %gui qt
# import numpy as np
# import surfer
# from mindboggle.utils.io_vtk import read_scalars as rs
# d,n=rs('/drop/MB/data/arno/shapes/travel_depth_rescaled.vtk')
# br = surfer.Brain('Twins-2-1', 'lh', 'inflated')
# br.add_data(np.array(d), min=0, max=1, alpha=0.5)


def vtkviewer(vtk_file_list, colormap_file=None):
    """
    Use vtkviewer to visualize one or more VTK surface files.

    Parameters
    ----------
    vtk_file_list : string or list of strings
        name of VTK surface mesh file or list of file names
    colormap_file : string
        name of Paraview-style XML colormap file

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.plots import vtkviewer
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT31.manual.vtk')
    >>> colormap_file = os.path.join(os.environ['MINDBOGGLE_TOOLS'], 'colormap.xml')
    >>> vtkviewer(vtk_file, colormap_file)

    """
    import os
    import glob
    import mindboggle.utils.vtkviewer as vv

    if isinstance(vtk_file_list, str):
        vtk_file_list = [vtk_file_list]

    if colormap_file:
        vtk_colormap = vv.VTKViewer.LoadColorMap(colormap_file)
    else:
        vtk_colormap = None

    vtkviewer = vv.VTKViewer()
    for vtk_file in vtk_file_list:
        fileNames = glob.glob(vtk_file)
        if len(fileNames) == 0:
            print "what:", vtk_file
        else:
            for fileName in fileNames:
                if os.path.isfile(fileName):
                    vtkviewer.AddFile(fileName,vtk_colormap)
                else:
                    print "what:", fileName
    vtkviewer.Start()


def plot_surfaces(vtk_file, mask_file='', mask_background=-1,
                  masked_output='', program='vtkviewer', colormap_file=None,
                  background_value=-1):
    """
    Use vtkviewer or mayavi2 to visualize VTK surface mesh data.

    Parameters
    ----------
    vtk_file : string
        name of VTK surface mesh file
    mask_file : string
        name of VTK surface mesh file to mask vtk_file vertices
    mask_background : integer
        mask background value
    masked_output : string
        temporary masked output file name
    program : string {'vtkviewer', 'mayavi2'}
        program to visualize VTK file
    colormap_file : string
        name of Paraview-style XML colormap file (for use with vtkviewer)
    background_value : integer
        background value

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.plots import plot_surfaces
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> mask_file = ''#os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> mask_background = -1
    >>> masked_output = ''
    >>> program = 'vtkviewer'
    >>> colormap_file = None
    >>> plot_surfaces(vtk_file, mask_file, mask_background, masked_output, program, colormap_file)

    """
    from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars
    from mindboggle.utils.utils import execute
    from mindboggle.utils.plots import vtkviewer

    if not program:
        program = 'vtkviewer'

    # Filter mesh with non-background values from a second (same-size) mesh:
    if mask_file:
        scalars, name = read_scalars(vtk_file, True, True)
        mask, name = read_scalars(mask_file, True, True)
        scalars[mask == mask_background] = -1
        if not masked_output:
            masked_output = 'temp.vtk'
        rewrite_scalars(vtk_file, masked_output, scalars) #, 'masked', mask)
        file_to_plot = masked_output
    else:
        file_to_plot = vtk_file

    # Display with vtkviewer.py:
    if program == 'vtkviewer':
        vtkviewer(file_to_plot, colormap_file)

    # Display with mayavi2:
    elif program == 'mayavi2':
        cmd = ["mayavi2", "-d", file_to_plot, "-m", "Surface", "&"]
        execute(cmd, 'os')


def plot_volumes(volume_files):
    """
    Use fslview to visualize image volume data.

    Parameters
    ----------
    volume_files : list of strings
        names of image volume files

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.plots import plot_volumes
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> volume_file1 = os.path.join(path, 'arno', 'mri', 't1weighted.nii.gz')
    >>> volume_file2 = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> volume_files = [volume_file1, volume_file2]
    >>> plot_volumes(volume_files)

    """
    from mindboggle.utils.utils import execute

    if isinstance(volume_files, str):
        volume_files = [volume_files]
    elif not isinstance(volume_files, list):
        import sys
        sys.error('plot_volumes() requires volume_files to be a list or string.')

    cmd = ["fslview"]
    cmd.extend(volume_files)
    cmd.extend('&')
    execute(cmd, 'os')


def plot_scalar_histogram(vtk_file, nbins=100):
    """
    Plot histogram of VTK surface mesh scalar values.

    Parameters
    ----------
    vtk_file : string
        name of VTK file with scalar values to plot
    nbins : integer
        number of histogram bins

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.plots import plot_scalar_histogram
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> plot_scalar_histogram(vtk_file, nbins=500)

    """
    import matplotlib.pyplot as plt
    from mindboggle.utils.io_vtk import read_scalars

    # Load values:
    values, name = read_scalars(vtk_file)

    # Histogram:
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.hist(values, nbins, normed=False, facecolor='gray', alpha=0.5)
    plt.show()


def plot_histograms(columns, column_name='', ignore_columns=[],
                    nbins=100, axis_limits=[], titles=[]):
    """
    Construct a histogram for each table column.

    Parameters
    ----------
    columns : list of lists
        list of lists of floats or integers
    column_name :  string
        column name
    ignore_columns : list of integers
        indices to table columns or sublists to exclude
    nbins : integer
        number of histogram bins
    axis_limits : list of four integers
        range of x- and y-axis ranges: [x-low, x-high, y-low, y-high]
    y_lim : list of two integers
        range of y-axis values
    titles : list of strings (length = number of columns - 1)
        histogram titles (if empty, use column headers)

    Examples
    --------
    >>> from mindboggle.utils.plots import plot_histograms
    >>> columns = [[1,1,2,2,2,2,2,2,3,3,3,4,4,8],[2,2,3,3,3,3,5,6,7]]
    >>> column_name = 'label: thickness: median (weighted)'
    >>> ignore_columns = []
    >>> nbins = 100
    >>> axis_limits = []
    >>> titles = ['title1','title2']
    >>> plot_histograms(columns, column_name, ignore_columns, nbins, axis_limits, titles)

    """
    import numpy as np
    import matplotlib.pyplot as plt

    ncolumns = len(columns)
    if ncolumns < 9:
        nplotrows = 1
        nplotcols = ncolumns
    else:
        nplotrows = np.ceil(np.sqrt(ncolumns))
        nplotcols = nplotrows

    #-------------------------------------------------------------------------
    # Construct a histogram from each column and display:
    #-------------------------------------------------------------------------
    fig = plt.figure()
    for icolumn, column in enumerate(columns):
        if icolumn not in ignore_columns:
            ax = fig.add_subplot(nplotrows, nplotcols, icolumn + 1)
            column = [np.float(x) for x in column]
            ax.hist(column, nbins, normed=False, facecolor='gray', alpha=0.5)
            plt.xlabel(column_name, fontsize='small')
            plt.ylabel('frequency')
            if len(titles) == ncolumns:
                plt.title(titles[icolumn], fontsize='small')
            else:
                plt.title(column_name + ' histogram', fontsize='small')
            if axis_limits:
                ax.axis(axis_limits)
    plt.show()


def plot_columns(columns, x_column, ignore_columns=[], plot_line=True,
                 connect_markers=True, title='', x_label='', y_label='',
                 legend=True, legend_labels=[]):
    """
    Scatter plot columns against the values of one of the columns.

    Parameters
    ----------
    columns : list of lists of numbers
        columns of data (all of the same length)
    x_column : list of numbers
        column of numbers against which other columns are plotted
    ignore_columns : list of integers
        indices to columns to exclude
    plot_line : Boolean
        plot identity line?
    connect_markers : Boolean
        connect markers?
    title :  string
        title
    x_label : string
        description of x_column
    y_label : string
        description of other columns
    legend : Boolean
        plot legend?
    legend_labels : list of strings (length = number of columns)
        legend labels

    Examples
    --------
    >>> from mindboggle.utils.plots import plot_columns
    >>> columns = [[1,1,2,2,2,3,3,4,4,8],[2,2,3,3,3,3,5,6,7,7]]
    >>> x_column = [1,1.5,2.1,1.8,2.2,3,3.1,5,7,6]
    >>> ignore_columns = []
    >>> plot_line = True
    >>> connect_markers = True
    >>> title = 'title'
    >>> x_label = 'xlabel'
    >>> y_label = 'ylabel'
    >>> legend = True
    >>> legend_labels = ['mark1','mark2']
    >>> plot_columns(columns, x_column, ignore_columns, plot_line, connect_markers, title, x_label, y_label, legend, legend_labels)

    """
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    import numpy as np

    ncolumns = len(columns)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    #color_indices = np.linspace(0, 1, num=ncolumns, endpoint=False)
    colors = ['b','r','c','m','k','g','y']
    #-------------------------------------------------------------------------
    # Scatter plot:
    #-------------------------------------------------------------------------
    hold = True
    if plot_line:
        min_value = np.inf
        max_value = -np.inf
    for icolumn, column in enumerate(columns):
        column = [np.float(x) for x in column]
        if icolumn not in ignore_columns:
            #color_index = color_indices[icolumn]
            #color_value = plt.cm.hsv(color_index * np.ones(len(column)))
            color_value = colors[icolumn]
            if len(legend_labels) == ncolumns:
                color_text = legend_labels[icolumn]
                if connect_markers and not plot_line:
                    plt.plot(x_column, column, '-', marker='o',
                             color=color_value, hold=hold,
                             label=color_text)
                else:
                    plt.scatter(x_column, column, marker='o', s=50,
                                color=color_value, hold=hold,
                                label=color_text)
            else:
                if connect_markers and not plot_line:
                    plt.plot(x_column, column, '-', marker='o',
                             color=color_value, hold=hold)
                else:
                    plt.scatter(x_column, column, marker='o', s=50,
                                color=color_value, hold=hold)

        if plot_line:
            if min(column) < min_value:
                min_value = min(column)
            if max(column) > max_value:
                max_value = max(column)

    #-------------------------------------------------------------------------
    # Add legend and display:
    #-------------------------------------------------------------------------
    if plot_line:
        plt.plot(range(int(min_value), int(max_value) + 2))
    if legend:
        fontP = FontProperties()
        fontP.set_size('small')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='lower right', prop=fontP)
    ax.grid()
    if x_label:
        plt.xlabel(x_label)
    if y_label:
        plt.ylabel(y_label)
    plt.title(title)
    plt.show()


if __name__== '__main__':

    import os
    from mindboggle.utils.io_table import select_column_from_mindboggle_tables
    from mindboggle.utils.plots import plot_histograms
    from mindboggle.utils.plots import plot_columns

    plot_hist = 0#True
    plot_scat = True

    hemi = 'left'
    tables_dir = os.path.join(os.environ['HOME'], 'mindboggled', 'tables')
    label_name = 'label'
    write_table = False
    output_table = ''
    delimiter = ','
    ignore_columns = []
    legend = True

    if plot_hist:

        subjects = ['UM0029_2R1_full',
                    #'UM0029UMMR2R1_repositioned',
                    'UM0029UMMR2R1_FS11212',
                    #'UM0029_2R1_half',
                    'UM0029UMMR2R1_segmented',
                    'UM0029UMMR1R1_FS11762',
                    'UM0029UMMR2R1_antsCorticalThickness',
                    'UM0029UMMR1R1_antsCorticalThickness']

        table_name = "vertices.csv"
        column_name = "thickness"
        tables, columns, column_name, row_names, row_names_title, \
        output_table = select_column_from_mindboggle_tables(subjects, hemi,
            tables_dir, table_name, column_name, label_name,
            write_table, output_table, delimiter)

        nbins = 100
        axis_limits = [0, 5, 0, 6000]
        plot_histograms(columns, column_name, ignore_columns, nbins, axis_limits, titles=subjects)

    if plot_scat:

        subjects = ['UM0029_2R1_full',
                    'UM0029UMMR2R1_FS11212',
                    'UM0029UMMR2R1_segmented',
                    'UM0029UMMR1R1_FS11762']
        subjects = ['UM0029UMMR2R1_antsCorticalThickness',
                    'UM0029UMMR1R1_antsCorticalThickness']

        table_name = "label_shapes.csv"
        column_name = "label: thickness: mean (weighted)"
        tables, columns, column_name, row_names, row_names_title, \
        output_table = select_column_from_mindboggle_tables(subjects, hemi,
            tables_dir, table_name, column_name, label_name,
            write_table, output_table, delimiter)

        # Columns against one column, with identity line:
        x_column = columns[0]
        columns1 = columns[1::]
        legend_labels = subjects[1::]
        plot_line = True
        connect_markers = False
        title = 'Median cortical region thicknesses against repositioned results'
        x_label = 'median thickness of repositioned'
        y_label = 'median thickness'
        plot_columns(columns1, x_column, ignore_columns, plot_line, connect_markers, title, x_label, y_label, legend, legend_labels)

        # Columns against labels, with no identity line:
        x_column = row_names
        legend_labels = subjects
        plot_line = False
        connect_markers = True
        title = 'Median cortical region thicknesses per label'
        x_label = label_name
        y_label = 'median thickness'
        plot_columns(columns, x_column, ignore_columns, plot_line, connect_markers, title, x_label, y_label, legend, legend_labels)
