"""
Example
-------

from pipeline import create_mindboggle_flow
wf = create_mindboggle_flow()
wf.inputs.feature_extractor.curvature_file = '/Users/arno/Documents/Projects/mindboggle/data/lh.curv'
wf.inputs.feature_extractor.surface_file = '/Users/arno/Documents/Projects/mindboggle/data/lh.pial'
#wf.inputs.convex_hull.surface_file = '/Users/arno/Documents/Projects/mindboggle/data/lh.basin.vtk'
wf.run() # doctest: +SKIP

"""

import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
import numpy as np

def extract_structures(surface_file, curvature_file):
    """Extract structures: basins, pits, and fundi

    extract_structures('../data/lh.curv','../data/lh.pial')
    """
    from glob import glob
    import subprocess as sp
    cmd = '/projects/mindboggle/mindboggle/feature_extraction/basins_fundi_pits/extract.py'
    cmd = ['python', cmd, '-T', '%s'%curvature_file, '%s'%surface_file]
    proc = sp.Popen(cmd)
    o, e = proc.communicate()
    if proc.returncode > 0 :
        raise Exception('\n'.join(['extract.py failed', o, e]))
    #output_file = glob('file1.vtk').pop()
    structure_files = glob('*.vtk')
    return structure_files

def construct_convex_hull(structure_file):
    """Generate convex hull from surface mesh

    convex_hull('../data/lh.basin.vtk')
    """
    import subprocess as sp
    cmd = '/projects/mindboggle/mindboggle/convex_hull.py'
    cmd = ['python', cmd, '%s'%structure_file]
    proc = sp.Popen(cmd)
    o, e = proc.communicate()
    if proc.returncode > 0 :
        raise Exception('\n'.join(['convex_hull.py failed', o, e]))
    convex_hull = 'convex_hull'
    return convex_hull

def measure_position(object):
    """Measure

    measure_()
    """
    from measure_shape.py import measure_position
    if type(object) is np.ndarray:
        measurement = measure_position(object)
    return measurement

def measure_extent(object):
    """Measure

    measure_()
    """
    from measure_shape.py import measure_extent
    if type(object) is np.ndarray:
        measurement = measure_extent(object)
    return measurement

def measure_depth(object):
    """Measure

    measure_()
    """
    from measure_shape.py import measure_depth
    if type(object) is np.ndarray:
        measurement = measure_(object)
    return measurement

def measure_curvature(object):
    """Measure

    measure_()
    """
    from measure_shape.py import measure_curvature
    if type(object) is np.ndarray:
        measurement = measure_curvature(object)
    return measurement

def measure_LaplaceBeltrami(object):
    """Measure

    measure_()
    """
    from measure_shape.py import measure_LaplaceBeltrami
    if type(object) is np.ndarray:
        measurement = measure_LaplaceBeltrami(object)
    return measurement

def measure_texture(object):
    """Measure

    measure_()
    """
    from measure_shape.py import measure_texture
    if type(object) is np.ndarray:
        measurement = measure_texture(object)
    return measurement

"""
def measure_(object):
    ""Measure

    measure_()
    ""
    from measure_shape.py import measure_
    if type(object) is np.ndarray:
        measurement = measure_(object)
    return measurement
"""

def write_structures_to_database(object):
    """Write to database

    write_to_database()
    """
    from write_to_database.py import structures_to_database
    if type(object) is np.ndarray:
        structures_to_database(object)
    success = 'True'
    return success

def write_measures_to_database(measurement):
    """Write to database

    write_to_database()
    """
    from write_to_database.py import measures_to_database
    if type(measurement) is np.ndarray:
        measures_to_database(measurement)
    success = 'True'
    return success

def create_mindboggle_flow():
    """Create a Mindboggle workflow"""

    # Create an instance of a workflow and indicate that it should operate in the current directory:
    flow = pe.Workflow(name='pipeline')
    flow.base_dir = '.'

    # Define nodes
    extract_structures_node = pe.Node(util.Function(input_names = ['surface_file', 'curvature_file'],
                                              output_names=['structure_files'],
                                              function = extract_structures),
                                name='extract_structures_node')
    measure_position_node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=['position'],
                                              function = measure_position),
                                       iterfield=['object'],
                                name='measure_position_node')
    measure_extent_node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=['extent'],
                                              function = measure_extent),
                                       iterfield=['object'],
                                name='measure_extent_node')
    measure_curvature_node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=['curvature'],
                                              function = measure_curvature),
                                       iterfield=['object'],
                                name='measure_curvature_node')
    measure_LaplaceBeltrami_node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=['LaplaceBeltrami'],
                                              function = measure_LaplaceBeltrami),
                                       iterfield=['object'],
                                name='measure_LaplaceBeltrami_node')
    measure_texture_node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=['texture'],
                                              function = measure_texture),
                                       iterfield=['object'],
                                name='measure_texture_node')
    """
    measure__node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=[''],
                                              function = measure_),
                                       iterfield=['object'],
                                name='position_')
    """
    construct_convex_hull_node = pe.MapNode(util.Function(input_names = ['structure_file'],
                                              output_names=['convex_hull'],
                                              function = construct_convex_hull),
                                       iterfield=['structure_file'],
                                name='construct_convex_hull_node')

    measure_hull_extent_node = pe.MapNode(util.Function(input_names = ['convex_hull'],
                                              output_names=['extent'],
                                              function = measure_extent),
                                       iterfield=['convex_hull'],
                                name='measure_hull_extent_node')
    measure_hull_depth_node = pe.MapNode(util.Function(input_names = ['convex_hull'],
                                              output_names=['depth'],
                                              function = measure_depth),
                                       iterfield=['convex_hull'],
                                name='measure_hull_depth_node')

    write_structures_to_database_node = pe.MapNode(util.Function(input_names = ['object'],
                                              output_names=['success'],
                                              function = write_structures_to_database),
                                       iterfield=['object'],
                                name='write_structures_to_database_node')
    write_measures_to_database_node = pe.MapNode(util.Function(input_names = ['measurement'],
                                              output_names=['success'],
                                              function = write_measures_to_database),
                                       iterfield=['measurement'],
                                name='write_measures_to_database_node')

# Define connections
    #flow.add_nodes([extract_structures_node])

    flow.connect(extract_structures_node, 'structure_files', write_structures_to_database_node, 'object')

    flow.connect(extract_structures_node, 'structure_files', measure_position_node, 'object')
    flow.connect(extract_structures_node, 'structure_files', measure_extent_node, 'object')
    flow.connect(extract_structures_node, 'structure_files', measure_curvature_node, 'object')
    flow.connect(extract_structures_node, 'structure_files', measure_LaplaceBeltrami_node, 'object')
    flow.connect(extract_structures_node, 'structure_files', measure_texture_node, 'object')
    #flow.connect(extract_structures_node, 'structure_files', measure__node, 'object')

    flow.connect(extract_structures_node, 'structure_files', construct_convex_hull_node, 'structure_file')

    flow.connect(construct_convex_hull_node, 'convex_hull', measure_hull_extent_node, 'convex_hull')
    flow.connect(construct_convex_hull_node, 'convex_hull', measure_hull_depth_node, 'convex_hull')

    """
    flow.connect(measure_position_node, 'position', write_measures_to_database_node, 'measurement')
    flow.connect(measure_extent_node, 'extent', write_measures_to_database_node, 'measurement')
    flow.connect(measure_depth_node, 'depth', write_measures_to_database_node, 'measurement')
    flow.connect(measure_curvature_node, 'curvature', write_measures_to_database_node, 'measurement')
    flow.connect(measure_LaplaceBeltrami_node, 'LaplaceBeltrami', write_measures_to_database_node, 'measurement')
    flow.connect(measure_texture_node, 'texture', write_measures_to_database_node, 'measurement')
    #flow.connect(measure__node, '', write_measures_to_database_node, 'measurement')
    """

    return flow

if __name__ == '__main__':
    flow = create_mindboggle_flow()
    flow.inputs.extract_structures_node.curvature_file = '/Users/arno/Documents/Projects/mindboggle/data/lh.curv'
    flow.inputs.extract_structures_node.surface_file = '/Users/arno/Documents/Projects/mindboggle/data/lh.pial'
    flow.write_graph()
    #flow.run()

