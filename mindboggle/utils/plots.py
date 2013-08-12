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
    ax = fig.add_subplot(111)
    ax.hist(values, nbins, normed=False, facecolor='gray', alpha=0.5)
    plt.show()


def scatter_plot_table(table_file, x_column=1, ignore_columns=[0]):
    """
    Scatter plot table columns against the values of one of the columns.

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.plots import scatter_plot_table
    >>> #path = os.path.join(os.environ['HOME'], 'mindboggled', 'tables')
    >>> table_file = os.path.join('/drop/test.csv')
    >>> x_column = 1
    >>> ignore_columns = [0]
    >>> scatter_plot_table(table_file, x_column, ignore_columns)

    """
    import os
    import csv
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import numpy as np

    #--------------------------------------------------------------------------
    # Extract columns from the table
    #--------------------------------------------------------------------------
    if not os.path.exists(table_file):
        raise(IOError(table_file + " not found"))

    reader = csv.reader(open(table_file, 'rb'),
                        delimiter='\t', quotechar='"')
    columns = [list(x) for x in zip(*reader)]
    ncolumns = len(columns)

    keep_columns = []
    for icolumn, column in enumerate(columns):
        if icolumn not in ignore_columns:
            keep_columns.extend(column[1::])
    max_value = int(max([np.float(x) for x in keep_columns])) + 1

    column_xaxis = columns[x_column][1::]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    jet = plt.get_cmap('jet', ncolumns)
    color_values = range(ncolumns)
    cNorm = colors.Normalize(vmin=0, vmax=color_values[-1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    hold = True
    for icolumn, column in enumerate(columns):
        if icolumn not in ignore_columns:

            color_text = column[0]
            color_value = scalarMap.to_rgba(color_values[icolumn])
            plt.plot(column_xaxis, column[1::], 'o', hold=hold,
                     color=color_value, label=color_text)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper left')
    ax.grid()
    plt.plot(range(max_value))
    plt.show()
