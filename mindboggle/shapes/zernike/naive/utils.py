'''VALID'''

from mindboggle.utils.io_vtk import read_vtk as old_read_vtk
from numpy import array

def read_vtk(filename) :

    faces, lines, indices, points, n_points, scalars, scalar_names, input_file = old_read_vtk(filename)
    return array(points), array(faces)
