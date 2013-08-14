#!/usr/bin/env python
"""
Plotting functions.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


#-----------------------------------------------------------------------------
# Example Pysurfer plot of FreeSurfer surface + VTK scalars
#-----------------------------------------------------------------------------
"""
ipython
%gui qt
import numpy as np
import surfer
from mindboggle.utils.io_vtk import read_scalars as rs
d,n=rs('/drop/MB/data/arno/shapes/travel_depth_rescaled.vtk')
br = surfer.Brain('Twins-2-1', 'lh', 'inflated')
br.add_data(np.array(d), min=0, max=1, alpha=0.5)
"""


def plot_vtk(vtk_file, mask_file='', masked_output=''):
    """
    Use mayavi2 to visualize VTK surface mesh data.

    Inputs
    ------
    vtk_file : string
        name of VTK surface mesh file

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.plots import plot_vtk
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> mask_file = os.path.join(path, 'arno', 'features', 'folds.vtk')
    >>> masked_output = ''
    >>> plot_vtk(vtk_file, mask_file, masked_output)

    """
    from mindboggle.utils.utils import execute

    # Filter mesh with the non -1 values from a second (same-size) mesh:
    if mask_file:

        from mindboggle.utils.io_vtk import read_scalars, rewrite_scalars

        scalars, name = read_scalars(vtk_file)
        mask, name = read_scalars(mask_file)
        if not masked_output:
            masked_output = 'temp.vtk'
        rewrite_scalars(vtk_file, masked_output, scalars, 'masked', mask)

        cmd = ["mayavi2", "-d", masked_output, "-m", "Surface"]

    else:

        cmd = ["mayavi2", "-d", vtk_file, "-m", "Surface"]
    cmd.extend('&')
    execute(cmd, 'os')


def plot_volumes(volume_files):
    """
    Use fslview to visualize image volume data.

    Inputs
    ------
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

    Inputs
    ------
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


def plot_histograms(data, ignore_entries=[0], delimiter=',', nbins=100,
                    legend=True):
    """
    Construct a histogram for each table column.

    Inputs
    ------
    data : string or list of lists
        name of comma-separated table file or list of lists of floats
    ignore_entries : list of integers
        indices to table columns or sublists to exclude
    delimiter : string
        delimiter for table file
    nbins : integer
        number of histogram bins
    legend : Boolean
        plot legend?

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.plots import plot_histograms
    >>> data = '/drop/LAB/left_label_thickness_medians_CLEAN.csv'
    >>> ignore_entries = [0]
    >>> delimiter = ','
    >>> nbins = 100
    >>> legend = True #False
    >>> plot_histograms(data, ignore_entries, delimiter, nbins, legend)

    """
    import os
    import sys
    import csv
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    from matplotlib.font_manager import FontProperties

    # Extract columns from the table:
    if isinstance(data, str):
        if os.path.exists(data):
            reader = csv.reader(open(data, 'rb'),
                                delimiter=delimiter, quotechar='"')
            columns = [list(x) for x in zip(*reader)]
        else:
            raise(IOError(data + " not found"))
    elif isinstance(data, list) and isinstance(data[0], list):
        columns = data
    else:
        sys.exit('Data need to be in a .csv table file or list of lists.')

    ncolumns = len(columns)
    if ncolumns < 4:
        nplotrows = 1
        nplotcols = ncolumns
    else:
        nplotrows = np.ceil(np.sqrt(len(columns)))
        nplotrows = nplotrows

    #-------------------------------------------------------------------------
    # Construct a histogram from each column and display:
    #-------------------------------------------------------------------------
    fig = plt.figure()
    for icolumn, column in enumerate(columns):
        if icolumn not in ignore_entries:

            ax = fig.add_subplot(nplotrows, nplotcols, icolumn)
            plot_column = [np.float(x) for x in column[1::]]
            ax.hist(plot_column[1::], nbins, normed=False, facecolor='gray',
                    alpha=0.5)
            plt.title(column[0], fontsize='small')
    plt.show()


def scatter_plot_from_table(data, x_column=0, ignore_columns=[0],
                            delimiter=',', legend=True):
    """
    Scatter plot table columns against the values of one of the columns.

    Inputs
    ------
    data : string or list of lists
        name of comma-separated table file or list of lists of floats
    x_column : integer
        index of column against which other columns are plotted
    ignore_columns : list of integers
        indices to columns to exclude
    delimiter : string
        delimiter for table_file
    legend : Boolean
        plot legend?

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.plots import scatter_plot_from_table
    >>> #path = os.path.join(os.environ['HOME'], 'mindboggled', 'tables')
    >>> #table_file = os.path.join('/drop/subjects_label_shapes.csv')
    >>> hemi = 'right'
    >>> s = '_CLEAN_bad'
    >>> data = os.path.join('/drop/'+hemi+'_label_thickness_medians'+s+'.csv')
    >>> x_column = 2
    >>> ignore_columns = [0,1]
    >>> delimiter = ','
    >>> legend = True #False
    >>> data = '/drop/test.csv'
    >>> scatter_plot_from_table(data, x_column, ignore_columns, delimiter, legend)

    """
    import os
    import sys
    import csv
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    from matplotlib.font_manager import FontProperties
    import numpy as np

    # Extract columns from the table:
    if isinstance(data, str):
        if os.path.exists(data):
            reader = csv.reader(open(data, 'rb'),
                                delimiter=delimiter, quotechar='"')
            columns = [list(x) for x in zip(*reader)]
        else:
            raise(IOError(data + " not found"))
    elif isinstance(data, list) and isinstance(data[0], list):
        columns = data
    else:
        sys.exit('Data need to be in a .csv table file or list of lists.')

    ncolumns = len(columns)

    keep_columns = []
    for icolumn, column in enumerate(columns):
        if icolumn not in ignore_columns:
            keep_columns.extend(column[1::])

    column_xaxis = columns[x_column][1::]

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    jet = plt.get_cmap('jet', ncolumns)
    color_values = range(ncolumns)
    cNorm = colors.Normalize(vmin=0, vmax=color_values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    #-------------------------------------------------------------------------
    # Scatter plot:
    #-------------------------------------------------------------------------
    hold = True
    for icolumn, column in enumerate(columns):
        if icolumn not in ignore_columns:

            color_text = column[0]
            color_value = scalarMap.to_rgba(color_values[icolumn])
            plt.plot(column_xaxis, column[1::], 'o', hold=hold,
                     color=color_value, label=color_text)
            if icolumn == x_column:
                plt.xlabel(color_text)

    #-------------------------------------------------------------------------
    # Add legend and display:
    #-------------------------------------------------------------------------
    if x_column not in ignore_columns:
        min_value = int(min([np.float(x) for x in keep_columns]))
        max_value = int(max([np.float(x) for x in keep_columns])) + 2
        plt.plot(range(min_value, max_value))
    if legend:
        fontP = FontProperties()
        fontP.set_size('small')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='lower right', prop=fontP)
    ax.grid()
    plt.show()
