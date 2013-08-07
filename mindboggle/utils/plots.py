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


#-----------------------------------------------------------------------------
# Plot VTK surface mesh
#-----------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------
# Plot image volume
#-----------------------------------------------------------------------------
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


#-----------------------------------------------------------------------------
# Plot histogram of VTK surface mesh scalar values
#-----------------------------------------------------------------------------
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

