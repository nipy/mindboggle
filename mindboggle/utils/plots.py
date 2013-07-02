#!/usr/bin/env python
"""
Plotting functions.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

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
    import os
#    import subprocess

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

# Note: subprocess won't allow me to put the command in the background:
#    p = subprocess.Popen(cmd)
#    p.communicate()
    cmd = ' '.join(cmd) + ' &'
    print(cmd)
    os.system(cmd)


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
    >>> from mindboggle.utils.plots import plot_volume
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> volume_file1 = os.path.join(path, 'arno', 'mri', 't1weighted.nii.gz')
    >>> volume_file2 = os.path.join(path, 'arno', 'mri', 't1weighted_brain.nii.gz')
    >>> volume_files = [volume_file1, volume_file2]
    >>> plot_volume(volume_files)

    """
    import os
#    import subprocess

    cmd = ["fslview", " ".join(volume_files)]

# Note: subprocess won't allow me to put the command in the background:
#    p = subprocess.Popen(cmd)
#    p.communicate()
    cmd = ' '.join(cmd) + ' &'
    print(cmd)
    os.system(cmd)


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
    ax.hist(values, nbins, normed=1, facecolor='gray', alpha=0.1)
    plt.show()

