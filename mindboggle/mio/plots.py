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
# from mindboggle.mio.vtks import read_scalars as rs
# d,n=rs('/drop/MB/data/arno/shapes/travel_depth_rescaled.vtk')
# br = surfer.Brain('Twins-2-1', 'lh', 'inflated')
# br.add_data(np.array(d), min=0, max=1, alpha=0.5)


def plot_surfaces(vtk_files, use_colormap=False, colormap_file=''):
    """
    Use vtkviewer to visualize one or more VTK surface files.

    Optionally provide colormap file or set $COLORMAP environment variable.

    Parameters
    ----------
    vtk_files : string or list of strings
        name of VTK surface mesh file or list of file names
    use_colormap : Boolean
        use Paraview-style XML colormap file?
    colormap_file : string
        use colormap in given file if use_colormap==True?  if empty and
        use_colormap==True, use file set by $COLORMAP environment variable

    Examples
    --------
    >>> import os
    >>> from mindboggle.mio.plots import plot_surfaces
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> label_file = fetch_data(urls['freesurfer_labels'])
    >>> os.rename(label_file, label_file + '.nii.gz')
    >>> label_file = label_file + '.nii.gz'
    >>> use_colormap = True
    >>> colormap_file = '/software/surface_cpp_tools/colormap.xml' # doctest: +SKIP
    >>> plot_surfaces(vtk_files, use_colormap, colormap_file) # doctest: +SKIP

    """
    import os
    import sys

    import mindboggle.thirdparty.vtkviewer as vtkviewer

    if isinstance(vtk_files, str):
        vtk_files = [vtk_files]

    vv = vtkviewer.VTKViewer()

    colormap = None
    if use_colormap:
        if colormap_file and os.path.isfile(colormap_file):
            colormap = vv.LoadColorMap(colormap_file)
        elif "COLORMAP" in os.environ:
            colormap = vv.LoadColorMap(os.environ["COLORMAP"])

    for vtk_file in vtk_files:
        if os.path.isfile(vtk_file):
            vv.AddFile(vtk_file, colormap)
        else:
            sys.exit("Huh?: {0}".format(vtk_file))

    vv.Start()


def plot_mask_surface(vtk_file, mask_file='', nonmask_value=-1,
                      masked_output='', remove_nonmask=False,
                      program='vtkviewer',
                      use_colormap=False, colormap_file=''):
    """
    Use vtkviewer or mayavi2 to visualize VTK surface mesh data.

    If a mask_file is provided, a temporary masked file is saved,
    and it is this file that is viewed.

    If using vtkviewer, optionally provide colormap file
    or set $COLORMAP environment variable.

    Parameters
    ----------
    vtk_file : string
        name of VTK surface mesh file
    mask_file : string
        name of VTK surface mesh file to mask vtk_file vertices
    nonmask_value : integer
        nonmask (usually background) value
    masked_output : string
        temporary masked output file name
    remove_nonmask : Boolean
        remove vertices that are not in mask? (otherwise assign nonmask_value)
    program : string {'vtkviewer', 'mayavi2'}
        program to visualize VTK file
    use_colormap : Boolean
        use Paraview-style XML colormap file set by $COLORMAP env variable?
    colormap_file : string
        use colormap in given file if use_colormap==True?  if empty and
        use_colormap==True, use file set by $COLORMAP environment variable

    Examples
    --------
    >>> import os
    >>> from mindboggle.mio.plots import plot_mask_surface
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['freesurfer_labels'])
    >>> os.rename(vtk_file, vtk_file + '.nii.gz')
    >>> vtk_file = vtk_file + '.nii.gz'
    >>> mask_file = ''
    >>> nonmask_value = 0 #-1
    >>> masked_output = ''
    >>> remove_nonmask = True
    >>> program = 'vtkviewer'
    >>> use_colormap = True
    >>> colormap_file = ''
    >>> plot_mask_surface(vtk_file, mask_file, nonmask_value, masked_output,
    ...     remove_nonmask, program, use_colormap, colormap_file) # doctest: +SKIP

    """
    import os
    import numpy as np

    from mindboggle.guts.mesh import remove_faces, reindex_faces_points
    from mindboggle.guts.utilities import execute
    from mindboggle.mio.plots import plot_surfaces
    from mindboggle.mio.vtks import read_scalars, rewrite_scalars, \
                                        read_vtk, write_vtk

    #-------------------------------------------------------------------------
    # Filter mesh with non-background values from a second (same-size) mesh:
    #-------------------------------------------------------------------------
    if mask_file:
        mask, name = read_scalars(mask_file, True, True)
        if not masked_output:
            masked_output = os.path.join(os.getcwd(), 'temp.vtk')
        file_to_plot = masked_output

        #---------------------------------------------------------------------
        # Remove nonmask-valued vertices:
        #---------------------------------------------------------------------
        if remove_nonmask:
            #-----------------------------------------------------------------
            # Load VTK files:
            #-----------------------------------------------------------------
            points, indices, lines, faces, scalars, scalar_names, npoints, \
                input_vtk = read_vtk(vtk_file, True, True)
            #-----------------------------------------------------------------
            # Find mask indices, remove nonmask faces, and reindex:
            #-----------------------------------------------------------------
            Imask = [i for i,x in enumerate(mask) if x != nonmask_value]
            mask_faces = remove_faces(faces, Imask)
            mask_faces, points, \
            original_indices = reindex_faces_points(mask_faces, points)
            #-----------------------------------------------------------------
            # Write VTK file with scalar values:
            #-----------------------------------------------------------------
            if np.ndim(scalars) == 1:
                scalar_type = type(scalars[0]).__name__
            elif np.ndim(scalars) == 2:
                scalar_type = type(scalars[0][0]).__name__
            else:
                print("Undefined scalar type!")
            write_vtk(file_to_plot, points, [], [], mask_faces,
                      scalars[original_indices].tolist(), scalar_names,
                      scalar_type=scalar_type)
        else:
            scalars, name = read_scalars(vtk_file, True, True)
            scalars[mask == nonmask_value] = nonmask_value
            rewrite_scalars(vtk_file, file_to_plot, scalars)
    else:
        file_to_plot = vtk_file

    #-------------------------------------------------------------------------
    # Display with vtkviewer.py:
    #-------------------------------------------------------------------------
    if program == 'vtkviewer':
        plot_surfaces(file_to_plot, use_colormap=use_colormap,
                      colormap_file=colormap_file)
    #-------------------------------------------------------------------------
    # Display with mayavi2:
    #-------------------------------------------------------------------------
    elif program == 'mayavi2':
        cmd = ["mayavi2", "-d", file_to_plot, "-m", "Surface", "&"]
        execute(cmd, 'os')


def plot_volumes(volume_files, command='fslview'):
    """
    Use fslview to visualize image volume data.

    Parameters
    ----------
    volume_files : list of strings
        names of image volume files
    command : string
        plotting software command

    Examples
    --------
    >>> import os
    >>> from mindboggle.mio.plots import plot_volumes
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> label_file1 = fetch_data(urls['freesurfer_labels'])
    >>> os.rename(label_file1, label_file1 + '.nii.gz')
    >>> label_file1 = label_file1 + '.nii.gz'
    >>> label_file2 = fetch_data(urls['freesurfer_labels'])
    >>> os.rename(label_file2, label_file2 + '.nii.gz')
    >>> label_file2 = label_file2 + '.nii.gz'
    >>> volume_files = [label_file1, label_file2]
    >>> command = 'fslview'
    >>> command = '/Applications/ITK-SNAP.app/Contents/MacOS/InsightSNAP'
    >>> plot_volumes(volume_files, command=command) # doctest: +SKIP

    """
    from mindboggle.guts.utilities import execute

    if isinstance(volume_files, str):
        volume_files = [volume_files]
    elif not isinstance(volume_files, list):
        raise(IOError('plot_volumes() requires volume_files to be a list or string.'))

    if not isinstance(command, str):
        raise(IOError('plot_volumes() requires command to be a string.'))
    else:
        command = [command]

    command.extend(volume_files)
    command.extend('&')
    execute(command, 'os')


def histogram_of_vtk_scalars(vtk_file, nbins=100):
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
    >>> from mindboggle.mio.plots import histogram_of_vtk_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> vtk_file = fetch_data(urls['left_mean_curvature'])
    >>> os.rename(vtk_file, vtk_file + '.nii.gz')
    >>> vtk_file = vtk_file + '.nii.gz'
    >>> histogram_of_vtk_scalars(vtk_file, nbins=500) # doctest: +SKIP

    """
    import matplotlib.pyplot as plt
    from mindboggle.mio.vtks import read_scalars

    # Load values:
    values, name = read_scalars(vtk_file)

    # Histogram:
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.hist(values, nbins, normed=False, facecolor='gray', alpha=0.5)
    plt.show()


def histograms_of_lists(columns, column_name='', ignore_columns=[],
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
    >>> from mindboggle.mio.plots import histograms_of_lists
    >>> columns = [[1,1,2,2,2,2,2,2,3,3,3,4,4,8],[2,2,3,3,3,3,5,6,7]]
    >>> column_name = 'label: thickness: median (weighted)'
    >>> ignore_columns = []
    >>> nbins = 100
    >>> axis_limits = []
    >>> titles = ['title1','title2']
    >>> histograms_of_lists(columns, column_name, ignore_columns, nbins,
    ...                     axis_limits, titles) # doctest: +SKIP

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
            if len(titles) == ncolumns:
                plt.title(titles[icolumn], fontsize='small')
            else:
                plt.title(column_name + ' histogram', fontsize='small')
            if axis_limits:
                ax.axis(axis_limits)
    plt.show()


def boxplots_of_lists(columns, xlabel='', ylabel='', ylimit=None, title=''):
    """
    Construct a box plot for each table column.

    Parameters
    ----------
    columns : list of lists
        list of lists of floats or integers
    xlabel : str
        x-axis label
    ylabel : str
        y-axis label
    ylimit : float
        maximum y-value
    title : str
        title

    Examples
    --------
    >>> from mindboggle.mio.plots import boxplots_of_lists
    >>> columns = [[1,1,2,2,2,2,2,2,3,3,3,4,4,8],[2,2,3,3,3,3,5,6,7],
    ...            [2,2,2.5,2,2,2,3,3,3,3,5,6,7]]
    >>> xlabel = 'xlabel'
    >>> ylabel = 'ylabel'
    >>> ylimit = None
    >>> title = 'title'
    >>> boxplots_of_lists(columns, xlabel, ylabel, ylimit, title) # doctest: +SKIP

    """
    import matplotlib.pyplot as plt

    #-------------------------------------------------------------------------
    # Construct a box plot from each column and display:
    #-------------------------------------------------------------------------
    plt.figure()
    plt.boxplot(columns, 1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ylim([0,ylimit])
    plt.title(title)
    plt.show()


def scatterplot_lists(y_columns, x_column, ignore_columns=[], plot_line=True,
                      connect_markers=True, mstyle='o', msize=1,
                      title='', x_label='', y_label='',
                      legend=True, legend_labels=[]):
    """
    Scatter plot columns against the values of one of the columns.

    Parameters
    ----------
    y_columns : list of lists of numbers
        columns of data (all of the same length)
    x_column : list of numbers
        column of numbers against which other columns are plotted
    ignore_columns : list of integers
        indices to y_columns to exclude
    plot_line : Boolean
        plot identity line?
    connect_markers : Boolean
        connect markers?
    mstyle : string
        marker style
    msize : integer
        marker size
    title :  string
        title
    x_label : string
        description of x_column
    y_label : string
        description of y_columns
    legend : Boolean
        plot legend?
    legend_labels : list of strings (length = number of y_columns)
        legend labels

    Examples
    --------
    >>> from mindboggle.mio.plots import scatterplot_lists
    >>> y_columns = [[1,1,2,2,2,3,3,4,4,8],[2,2,3,3,3,3,5,6,7,7]]
    >>> x_column = [1,1.5,2.1,1.8,2.2,3,3.1,5,7,6]
    >>> ignore_columns = []
    >>> plot_line = True
    >>> connect_markers = True
    >>> mstyle = 'o'
    >>> msize = 10
    >>> title = 'title'
    >>> x_label = 'xlabel'
    >>> y_label = 'ylabel'
    >>> legend = True
    >>> legend_labels = ['mark1','mark2']
    >>> scatterplot_lists(y_columns, x_column, ignore_columns, plot_line,
    ...                   connect_markers, mstyle, msize, title, x_label,
    ...                   y_label, legend, legend_labels) # doctest: +SKIP

    """
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    import matplotlib.cm as cm
    import numpy as np

    ncolumns = len(y_columns)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    #colors = ['b','r','c','m','k','g','y']
    colors = iter(cm.hsv(np.linspace(0, 1, ncolumns)))
    #-------------------------------------------------------------------------
    # Scatter plot:
    #-------------------------------------------------------------------------
    hold = True
    if plot_line:
        min_value = np.inf
        max_value = -np.inf
    for icolumn, column in enumerate(y_columns):
        column = [np.float(x) for x in column]
        if icolumn not in ignore_columns:
            color = next(colors)
            #color = colors[icolumn]
            if len(legend_labels) == ncolumns:
                color_text = legend_labels[icolumn]
                if connect_markers and not plot_line:
                    plt.plot(x_column, column, '-', marker=mstyle, s=msize,
                             facecolors='none', edgecolors=color, hold=hold,
                             label=color_text)
                else:
                    plt.scatter(x_column, column, marker=mstyle, s=msize,
                                facecolors='none', edgecolors=color, hold=hold,
                                label=color_text)
            else:
                if connect_markers and not plot_line:
                    plt.plot(x_column, column, '-', marker=mstyle, s=msize,
                             facecolors='none', edgecolors=color, hold=hold)
                else:
                    plt.scatter(x_column, column, marker=mstyle, s=msize,
                                facecolors='none', edgecolors=color, hold=hold)

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


def scatterplot_list_pairs(columns, ignore_first_column=False, plot_line=True,
                           connect_markers=True, mstyle='o', msize=1,
                           mcolor='', title='', x_label='', y_label='',
                           limit=None, legend=True, legend_labels=[]):
    """
    Scatter plot pairs of columns.

    Parameters
    ----------
    columns : list of lists of numbers
        alternating columns of data (all of the same length)
    ignore_first_column : Boolean
        exclude first column?
    plot_line : Boolean
        plot identity line?
    connect_markers : Boolean
        connect markers?
    mstyle : string
        marker style
    msize : integer
        marker size
    mcolor : string
        marker color (if empty, generate range of colors)
    title :  string
        title
    x_label : string
        description of x_column
    y_label : string
        description of other columns
    limit : float
        x- and y-axis extent
    legend : Boolean
        plot legend?
    legend_labels : list of strings (length = number of columns)
        legend labels

    Examples
    --------
    >>> from mindboggle.mio.plots import scatterplot_list_pairs
    >>> columns = [['labels'], [1,1,2,2,2,3,3,4,4,8],[2,2,3,3,3,3,5,6,7,7],
    ...            [1,1.5,2.1,1.8,2.2,3,3.1,5,7,6],
    ...            [1.2,0.5,2,1.3,1.2,3,1,5.2,4,4.5]]
    >>> ignore_first_column = True
    >>> plot_line = True
    >>> connect_markers = True
    >>> mstyle = 'o'
    >>> msize = 10
    >>> mcolor = ''
    >>> title = 'title'
    >>> x_label = 'xlabel'
    >>> y_label = 'ylabel'
    >>> limit = None
    >>> legend = True
    >>> legend_labels = ['mark1','mark2']
    >>> scatterplot_list_pairs(columns, ignore_first_column, plot_line,
    ...                        connect_markers, mstyle, msize, mcolor, title,
    ...                        x_label, y_label, limit, legend, legend_labels) # doctest: +SKIP

    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.font_manager import FontProperties
    import numpy as np

    ncolumns = len(columns)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    if not mcolor:
        colors = iter(cm.hsv(np.linspace(0, 1, ncolumns)))
    #-------------------------------------------------------------------------
    # Scatter plot:
    #-------------------------------------------------------------------------
    hold = True
    if plot_line:
        min_value = np.inf
        max_value = -np.inf
    if ignore_first_column:
        columns = columns[1::]
    columns1 = [x for i,x in enumerate(columns) if np.mod(i,2) == 1]
    columns2 = [x for i,x in enumerate(columns) if np.mod(i,2) == 0]
    if not limit:
        limit = np.ceil(np.max([np.max(columns1), np.max(columns1)]))
    for icolumn, column1 in enumerate(columns1):
        column2 = columns2[icolumn]
        column1 = [np.float(x) for x in column1]
        column2 = [np.float(x) for x in column2]
        if mcolor:
            color = mcolor
        else:
            color = next(colors)
        if len(legend_labels) == ncolumns:
            if connect_markers and not plot_line:
                plt.plot(column1, column2, '-', marker=mstyle,
                         color=color, hold=hold,
                         label=legend_labels[icolumn])
            else:
                plt.scatter(column1, column2, marker=mstyle, s=msize,
                            facecolors='none', edgecolors=color, hold=hold,
                            label=legend_labels[icolumn])
        else:
            if connect_markers and not plot_line:
                plt.plot(column1, column2, '-', marker=mstyle,
                         color=color, hold=hold)
            else:
                plt.scatter(column1, column2, marker=mstyle, s=msize,
                            facecolors='none', edgecolors=color, hold=hold)

        if plot_line:
            if min(column1) < min_value:
                min_value = min(column1)
            if max(column1) > max_value:
                max_value = max(column1)

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
    plt.xlim([0,limit])
    plt.ylim([0,limit])
    ax.grid()
    ax.set_aspect(aspect='equal')
    if x_label:
        plt.xlabel(x_label)
    if y_label:
        plt.ylabel(y_label)
    plt.title(title)
    plt.show()


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod()


# #-----------------------------------------------------------------------------
# # Example: Plot scan-rescan thickness values from table (alternating columns)
# #-----------------------------------------------------------------------------
# if __name__ == '__main__':
#
#     import os
#     import numpy as np
#     from mindboggle.mio.plots import scatterplot_list_pairs
#     from mindboggle.mio.plots import boxplots_of_lists
#     from mindboggle.mio.plots import histograms_of_lists
#
#     #-------------------------------------------------------------------------
#     # Load thickness values from table (alternating columns are scan-rescan):
#     #-------------------------------------------------------------------------
#     tablename = '/desk/th_embarc40_antslabels.csv'
#     title = 'Thickinthehead, ANTs labels propagated through FS+ANTs gray (62 labels, 40 EMBARC controls)'
#
#     f1 = open(tablename,'r')
#     f1 = f1.readlines()
#     columns = [[] for x in f1[0].split()]
#     for row in f1:
#         row = row.split()
#         for icolumn, column in enumerate(row):
#             columns[icolumn].append(np.float(column))
#     columns1 = [columns[0]]
#     for i in range(1, len(columns)):
#         if np.mod(i,2) == 1:
#             columns1.append(columns[i])
#
#     #-------------------------------------------------------------------------
#     # Scatter plot:
#     #-------------------------------------------------------------------------
#     scat = True
#     if scat:
#         ignore_first_column = True
#         plot_line = False
#         connect_markers = False
#         mstyle = 'o'
#         msize = 10
#         mcolor = 'black'
#         xlabel = 'scan thickness (mm)'
#         ylabel = 'rescan thickness (mm)'
#         limit = 6.5
#         legend = True
#         legend_labels = ['mark1','mark2']
#         scatterplot_list_pairs(columns, ignore_first_column, plot_line,
#                                connect_markers, mstyle, msize, mcolor, title,
#                                xlabel, ylabel, limit, legend, legend_labels)
#
#     #-------------------------------------------------------------------------
#     # Box plot per label:
#     #-------------------------------------------------------------------------
#     box_per_label = True
#     if box_per_label:
#         rows1 = []
#         for row in f1:
#             rows1.append([np.float(x) for x in row.split()[1::]])
#         xlabel = 'label index'
#         ylabel = 'thickness (mm)'
#         ylimit = 6.5
#         boxplots_of_lists(rows1, xlabel, ylabel, ylimit, title)
#
#     #-------------------------------------------------------------------------
#     # Box plot per scan:
#     #-------------------------------------------------------------------------
#     box_per_scan = True #False
#     if box_per_scan:
#         xlabel = 'scan index'
#         ylabel = 'thickness (mm)'
#         ylimit = 6.5
#         boxplots_of_lists(columns1[1::], xlabel, ylabel, ylimit, title)
#
#     #-------------------------------------------------------------------------
#     # Histogram:
#     #-------------------------------------------------------------------------
#     hist = False
#     if hist:
#         ignore_columns = [0]
#         nbins = 10
#         axis_limits = [0, 5, 0, 10]
#         titles = []
#         histograms_of_lists(columns, 'thickness', ignore_columns, nbins,
#                             axis_limits, titles)
