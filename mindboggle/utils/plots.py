#!/usr/bin/env python
"""
Plotting functions.

Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#------------------------------------------------------------------------------
# Plot VTK surface mesh
#------------------------------------------------------------------------------
def plot_vtk(vtk_file):
    """
    Use mayavi2 to visualize VTK surface mesh data.

    Inputs
    ------
    vtk_file : string
        name of VTK surface mesh file
    """
    import os
#    import subprocess

    cmd = ["mayavi2", "-d", vtk_file, "-m", "Surface"]
    print(' '.join(cmd))

# Note: subprocess won't allow me to put the command in the background:
#    p = subprocess.Popen(cmd)
#    p.communicate()
    c = ' '.join(cmd) + ' &'
    print(c); os.system(c)


#------------------------------------------------------------------------------
# Plot histogram of VTK surface mesh scalar values
#------------------------------------------------------------------------------
def plot_scalar_histogram(vtk_file, nbins=100):
    """
    Plot histogram of VTK surface mesh scalar values.

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.mesh import plot_scalar_histogram
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> vtk_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
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

