#!/usr/bin/python

"""
This is Mindboggle's NiPype pipeline!

Example usage:

wf.run() # doctest: +SKIP


Authors:  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

import os
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.utility as util     # utility
import nipype.interfaces.io as nio
import numpy as np

from atlases import register_template, register_atlas, multilabel
from features import *

use_freesurfer_surfaces = 1
hemis = ['lh','rh']
surface_types = ['pial'] #,'inflated']
template_id = 'KKI'
template_name = template_id + '_2.tif'
template_reg_name = 'sphere_to_' + template_id + '_template.reg'
atlas_annot_name = 'aparcNMMjt.annot'

# Subjects
subjects_list = ['KKI2009-11', 'KKI2009-14']

use_linux_paths = 1
if use_linux_paths:
    subjects_path = '/usr/local/freesurfer/subjects'
else:
    subjects_path = '/Applications/freesurfer/subjects'

# Paths
templates_path = '/projects/mindboggle/data/templates_freesurfer'
atlases_path = subjects_path

# Output directory
results_path = '/projects/mindboggle/results/'
working_path = results_path + 'workingdir'
if not os.path.isdir(results_path):
    os.makedirs(results_path)
if not os.path.isdir(working_path):
    os.makedirs(working_path)

# Commands
depth_command = './measure/surface_measures/bin/travel_depth/TravelDepthMain'
curvature_command = './measure/surface_measures/bin/curvature/CurvatureMain'
extract_fundi_command = './extract/fundi/vtk_extract.py'

# List of atlas subjects
print("\nTEST ATLAS LIST\n")
atlas_list_file = os.path.join(atlases_path, 'MMRR.txt')
f = open(atlas_list_file)
atlas_list_lines = f.readlines()
atlas_names = [a.strip("\n") for a in atlas_list_lines if a.strip("\n")]

##############################################################################
#
#   Mindboggle workflow combining:
#   * Multi-atlas registration-based labeling workflow
#   * Feature-based labeling and shape analysis workflow
#   * Analytics
#
##############################################################################
mbflow = pe.Workflow(name='Mindboggle_workflow')
mbflow.base_dir = working_path

# Iterate inputs over subjects, hemispheres, surface types
infosource = pe.Node(interface=util.IdentityInterface(fields=['subject_id',
                                                              'hemi',
                                                              'surface_type']),
                     name = 'Inputs')
infosource.iterables = ([('subject_id', subjects_list),
                         ('hemi', hemis)])
datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id',
                                                         'hemi'],
                                               outfields=['surface_files',
                                               'inf_surface_files',
                                               'sph_surface_files']),
                     name = 'Surfaces')

# Specify the location and structure of the inputs and outputs
datasource.inputs.base_directory = subjects_path
datasource.inputs.template = '%s/surf/%s.%s'
datasource.inputs.template_args['surface_files'] = [['subject_id', 
                                                     'hemi', 
                                                     'pial']]
datasource.inputs.template_args['inf_surface_files'] = [['subject_id', 
                                                         'hemi', 
                                                         'inflated']]
datasource.inputs.template_args['sph_surface_files'] = [['subject_id', 
                                                         'hemi', 
                                                         'sphere']]
datasink = pe.Node(nio.DataSink(), name = 'Results')
datasink.inputs.base_directory = results_path
datasink.inputs.container = 'output'

# Connect input nodes
mbflow.connect([(infosource, datasource, 
                 [('subject_id','subject_id'),
                  ('hemi','hemi')])])

##############################################################################
#   Surface input and conversion
##############################################################################

# Convert FreeSurfer surfaces to VTK format
if use_freesurfer_surfaces:

    import nipype.interfaces.freesurfer as fs

    surface_conversion = pe.MapNode(fs.MRIsConvert(out_datatype='vtk'),
                                                   iterfield=['in_file'],
                                    name = 'Convert_surface')
    mbflow.connect([(datasource, surface_conversion, 
                     [('surface_files','in_file')])])

##############################################################################
#
#   Multi-atlas registration-based labeling workflow
#
##############################################################################
flo1 = pe.Workflow(name='Atlas_workflow')

##############################################################################
#   Multi-atlas registration
##############################################################################

# Template registration
template_reg = pe.Node(util.Function(input_names=['hemi',
                                                  'sph_surface_file',
                                                  'template_name',
                                                  'templates_path',
                                                  'template_reg_name'],
                                     output_names=['template_reg_name'],
                                     function = register_template),
                       name = 'Register_template')
template_reg.inputs.template_name = template_name
template_reg.inputs.templates_path = templates_path
template_reg.inputs.template_reg_name = template_reg_name

# Atlas registration
atlas_reg = pe.MapNode(util.Function(input_names=['hemi',
                                                  'subject_id',
                                                  'template_reg_name',
                                                  'atlas_name',
                                                  'atlases_path',
                                                  'atlas_annot_name'],
                                     output_names=['output_file'],
                                     function = register_atlas),
                       iterfield = ['atlas_name'],
                       name='Register_atlases')
atlas_reg.inputs.atlas_name = atlas_names
atlas_reg.inputs.atlases_path = atlases_path
atlas_reg.inputs.atlas_annot_name = atlas_annot_name

# Majority vote labeling
majority_vote = pe.Node(util.Function(input_names=['surface_file',
                                                   'annot_files'],
                                      output_names=['output_files'],
                                      function = multilabel),
                        name='Vote_majority')

# Add and connect the above nodes
flo1.add_nodes([template_reg, atlas_reg, majority_vote])

mbflow.connect([(infosource, flo1,
                 [('hemi', 'Register_template.hemi')]),
                (datasource, flo1, 
                 [('sph_surface_files',
                   'Register_template.sph_surface_file')])])

mbflow.connect([(infosource, flo1, 
                 [('hemi', 'Register_atlases.hemi'),
                  ('subject_id', 'Register_atlases.subject_id')])])
flo1.connect([(template_reg, atlas_reg, 
               [('template_reg_name', 'template_reg_name')])])

if use_freesurfer_surfaces:
    mbflow.connect([(surface_conversion, flo1, 
                     [('converted', 'Vote_majority.surface_file')])])
else:
    mbflow.connect([(datasource, flo1, 
                     [('surface_files','Vote_majority.surface_file')])])
flo1.connect([(atlas_reg, majority_vote,
               [('output_file', 'annot_files')])])
#mbflow.connect([(flo1, datasink,
#                 [('Register_atlases.output_file', 'atlas_registration')])])
mbflow.connect([(flo1, datasink,
                 [('Register_template.template_reg_name', 'test')])])
mbflow.connect([(flo1, datasink,
                 [('Vote_majority.output_files', 'maxlabels')])])

##############################################################################
#
#   Feature-based labeling and shape analysis workflow
#
##############################################################################

flo2 = pe.Workflow(name='Feature_workflow')

##############################################################################
#   Surface calculations
##############################################################################
"""
# Measure surface depth and curvature nodes
surface_depth = pe.MapNode(util.Function(input_names = ['command',
                                                        'surface_file'],
                                         output_names = ['depth_file'],
                                         function = measure_surface_depth),
                           iterfield = ['surface_file'],
                           name='Measure_depth')
surface_depth.inputs.command = depth_command

surface_curvature = pe.MapNode(util.Function(input_names = ['command',
                                                            'surface_file'],
                                    output_names = ['mean_curvature_file',
                                                    'gauss_curvature_file',
                                                    'max_curvature_file',
                                                    'min_curvature_file'],
                                    function = measure_surface_curvature),
                               iterfield = ['surface_file'],
                               name='Measure_curvature')
surface_curvature.inputs.command = curvature_command
"""
"""
# Add and connect nodes
flo2.add_nodes([surface_depth, surface_curvature])

if use_freesurfer_surfaces:
    flo2.connect([(surface_conversion, surface_depth, 
                   [('converted', 'surface_file')])])
    flo2.connect([(surface_conversion, surface_curvature, 
                   [('converted', 'surface_file')])])
else:
    # Connect input to surface depth and curvature nodes
    mbflow.connect([(flo1, flo2, 
                     [('Surfaces.surface_files',
                       'Measure_depth.surface_file')])])
    mbflow.connect([(flo1, flo2, 
                     [('Surfaces.surface_files',
                       'Measure_curvature.surface_file')])])

# Save
mbflow.connect([(flo1, datasink,
                 [('Vote_majority.output_files', 'surfaces.@depth')])])
mbflow.connect([(flo2, datasink,
                 [('Measure_curvature.mean_curvature_file', 
                   'surfaces.@mean_curvature'),
                  ('Measure_curvature.gauss_curvature_file', 
                   'surfaces.@gauss_curvature'),
                  ('Measure_curvature.max_curvature_file', 
                   'surfaces.@max_curvature'),
                  ('Measure_curvature.min_curvature_file', 
                   'surfaces.@min_curvature')])])

##############################################################################
#   Feature extraction
##############################################################################

# Extract features
fundus_extraction = pe.Node(util.Function(input_names = ['command',
                                                         'depth_file'],
                                          output_names = ['fundi'],
                                          function = extract_fundi),
                            name='Extract_fundi')
fundus_extraction.inputs.command = extract_fundi_command
"""
"""
sulcus_extraction = pe.Node(util.Function(input_names = ['depth_file',
                                                         'mean_curv_file',
                                                         'gauss_curv_file'],
                                          output_names = ['sulci'],
                                          function = extract_sulci),
                            name='Extract_sulci')

midaxis_extraction = pe.Node(util.Function(input_names = ['depth_file',
                                                          'mean_curv_file',
                                                          'gauss_curv_file'],
                                           output_names = ['midaxis'],
                                           function = extract_midaxis),
                             name='Extract_midaxis')

"""
"""
# Connect surface depth to feature extraction nodes
flo2.connect([(surface_depth, fundus_extraction, 
               [('depth_file', 'depth_file')])])
flo2.connect([(surface_depth, datasink, 
               [('depth_file', 'surface_depth')])])
"""
"""
flo2.connect([(surfaces, sulcus_extraction, 
               [('depth_file', 'depth_file'),
                ('mean_curv_file', 'mean_curv_file'),
                ('gauss_curv_file', 'gauss_curv_file')])])
flo2.connect([(surfaces, midaxis_extraction, 
               [('depth_file', 'depth_file'),
                ('mean_curv_file', 'mean_curv_file'),
                ('gauss_curv_file', 'gauss_curv_file')])])

##############################################################################
#   Label propagation
##############################################################################

# Label propagation node
label_propagation = pe.Node(util.Function(input_names=['labels', 'fundi'],
                                          output_names=['labels'],
                                          function = propagate_labels),
                            name='Propagate_labels')

# Volume label propagation node
volume_propagation = pe.Node(util.Function(input_names=['labels'],
                                           output_names=['labels'],
                                           function = propagate_volume_labels),
                             name='Propagate_volume_labels')

# Labeled surface patch and volume extraction nodes
patch_extraction = pe.Node(util.Function(input_names=['labels'],
                                         output_names=['patches'],
                                         function = extract_patches),
                           name='Extract_patches')

region_extraction = pe.Node(util.Function(input_names=['labels'],
                                          output_names=['regions'],
                                          function = extract_regions),
                            name='Extract_regions')

# Connect multiatlas registration(-based labeling) to label propagation nodes
flo2.connect([(atlas_reg, label_propagation, [('labels','labels')]),
              (fundus_extraction, label_propagation, [('fundi','fundi')])])

# Connect label propagation to labeled surface patch and volume extraction nodes
flo2.connect([(label_propagation, volume_propagation, [('labels', 'labels')])])
flo2.connect([(volume_propagation, region_extraction, [('labels', 'labels')])])
flo2.connect([(label_propagation, patch_extraction, [('labels', 'labels')])])

##############################################################################
#   Feature segmentation / identification
##############################################################################

# Feature segmentation nodes
sulcus_segmentation = pe.Node(util.Function(input_names=['sulci','labels'],
                                            output_names=['segmented_sulci'],
                                            function = segment_sulci),
                              name='Segment_sulci')

fundus_segmentation = pe.Node(util.Function(input_names=['fundi','labels'],
                                            output_names=['segmented_fundi'],
                                            function = segment_fundi),
                              name='Segment_fundi')

midaxis_segmentation = pe.Node(util.Function(input_names=['midaxis','labels'],
                                             output_names=['segmented_midaxis'],
                                             function = segment_midaxis),
                               name='Segment_midaxis')

# Connect feature and feature segmentation nodes
flo2.connect([(sulcus_extraction, sulcus_segmentation, [('sulci','sulci')]),
              (fundus_extraction, fundus_segmentation, [('fundi','fundi')]),
              (midaxis_extraction, midaxis_segmentation, [('midaxis','midaxis')])])

# Connect multiatlas registration(-based labeling) and feature segmentation nodes
flo2.connect([(label_propagation, sulcus_segmentation, [('labels','labels')]),
              (label_propagation, fundus_segmentation, [('labels','labels')]),
              (sulcus_segmentation, midaxis_segmentation, [('segmented_sulci','labels')])])
              
##############################################################################
#   Shape measurement
##############################################################################

# Shape measurement nodes
positions = pe.Node(util.Function(input_names = ['segmented_sulci', 
                                                 'segmented_fundi',
                                                 'segmented_midaxis', 
                                                 'pits', 
                                                 'patches', 
                                                 'regions'],
                                  output_names=['positions_sulci', 
                                                'positions_fundi',
                                                'positions_midaxis', 
                                                'positions_pits', 
                                                'positions_patches', 
                                                'positions_regions'],
                                  function = measure_positions),
                    name='Measure_positions')

extents = pe.Node(util.Function(input_names = ['segmented_sulci', 
                                               'segmented_fundi',
                                               'segmented_midaxis', 
                                               'pits', 
                                               'patches', 
                                               'regions'],
                                output_names=['extents_sulci', 
                                              'extents_fundi',
                                              'extents_midaxis', 
                                              'extents_pits',
                                              'extents_patches', 
                                              'extents_regions'],
                                function = measure_extents),
                    name='Measure_extents')

curvatures = pe.Node(util.Function(input_names = ['segmented_sulci', 
                                                  'segmented_fundi',
                                                  'segmented_midaxis', 
                                                  'pits', 
                                                  'patches', 
                                                  'regions'],
                                   output_names=['curvatures_sulci', 
                                                 'curvatures_fundi',
                                                 'curvatures_midaxis', 
                                                 'curvatures_pits',
                                                 'curvatures_patches', 
                                                 'curvatures_regions'],
                                   function = measure_curvatures),
                    name='Measure_curvatures')

depths = pe.Node(util.Function(input_names = ['segmented_sulci',
                                              'segmented_fundi',
                                              'segmented_midaxis',
                                              'pits', 
                                              'patches', 
                                              'regions'],
                               output_names=['depths_sulci', 
                                             'depths_fundi',
                                             'depths_midaxis',
                                             'depths_pits',
                                             'depths_patches', 
                                             'depths_regions'],
                               function = measure_depths),
                 name='Measure_depths')

spectra = pe.Node(util.Function(input_names = ['segmented_sulci',
                                               'segmented_fundi',
                                               'segmented_midaxis',
                                               'pits', 
                                               'patches', 
                                               'regions'],
                                output_names=['spectra_sulci',
                                              'spectra_fundi',
                                              'spectra_midaxis',
                                              'spectra_pits',
                                              'spectra_patches', 
                                              'spectra_regions'],
                                function = measure_spectra),
                  name='Measure_spectra')

# Connect labeled surface patches and volumes to shape measurement nodes
flo2.connect([(patch_extraction,  positions, [('patches', 'patches')])])
flo2.connect([(region_extraction, positions, [('regions', 'regions')])])
flo2.connect([(patch_extraction,  extents, [('patches', 'patches')])])
flo2.connect([(region_extraction, extents, [('regions', 'regions')])])
flo2.connect([(patch_extraction,  depths, [('patches', 'patches')])])
flo2.connect([(region_extraction, depths, [('regions', 'regions')])])
flo2.connect([(patch_extraction,  curvatures, [('patches', 'patches')])])
flo2.connect([(region_extraction, curvatures, [('regions', 'regions')])])
flo2.connect([(patch_extraction,  spectra, [('patches', 'patches')])])
flo2.connect([(region_extraction, spectra, [('regions', 'regions')])])

# Connect feature to shape measurement nodes
flo2.connect([(sulcus_segmentation, positions, [('segmented_sulci', 'segmented_sulci')])])
flo2.connect([(fundus_segmentation, positions, [('segmented_fundi', 'segmented_fundi')])])
flo2.connect([(pit_extraction, positions, [('pits', 'pits')])])
flo2.connect([(midaxis_segmentation, positions, [('segmented_midaxis', 'segmented_midaxis')])])

flo2.connect([(sulcus_segmentation, extents, [('segmented_sulci', 'segmented_sulci')])])
flo2.connect([(fundus_segmentation, extents, [('segmented_fundi', 'segmented_fundi')])])
flo2.connect([(midaxis_segmentation, extents, [('segmented_midaxis', 'segmented_midaxis')])])

flo2.connect([(sulcus_segmentation, curvatures, [('segmented_sulci', 'segmented_sulci')])])
flo2.connect([(fundus_segmentation, curvatures, [('segmented_fundi', 'segmented_fundi')])])
flo2.connect([(pit_extraction, curvatures, [('pits', 'pits')])])
flo2.connect([(midaxis_segmentation, curvatures, [('segmented_midaxis', 'segmented_midaxis')])])

flo2.connect([(sulcus_segmentation, depths, [('segmented_sulci', 'segmented_sulci')])])
flo2.connect([(fundus_segmentation, depths, [('segmented_fundi', 'segmented_fundi')])])
flo2.connect([(pit_extraction, depths, [('pits', 'pits')])])
flo2.connect([(midaxis_segmentation, depths, [('segmented_midaxis', 'segmented_midaxis')])])

flo2.connect([(sulcus_segmentation, spectra, [('segmented_sulci', 'segmented_sulci')])])
flo2.connect([(fundus_segmentation, spectra, [('segmented_fundi', 'segmented_fundi')])])
flo2.connect([(midaxis_segmentation, spectra, [('segmented_midaxis', 'segmented_midaxis')])])

##############################################################################
#    Store surface maps, features, and measures in database
##############################################################################

# Database nodes
maps_database = pe.Node(util.Function(input_names = ['depth_file',
                                                     'mean_curv_file',
                                                     'gauss_curv_file'],
                                      output_names=['success'],
                                      function = write_surfaces_to_database),
                        name='Write_surfaces_to_database')

features_database = pe.Node(util.Function(input_names = ['segmented_sulci',
                                                         'segmented_fundi',
                                                         'pits',
                                                         'segmented_midaxis'],
                                          output_names=['success'],
                                          function = write_features_to_database),
                            name='Write_features_to_database')

measures_database = pe.Node(util.Function(input_names = ['positions_sulci',
                                                         'positions_fundi',
                                                         'positions_pits',
                                                         'positions_midaxis',
                                                         'positions_patches',
                                                         'positions_regions',
                                                         'extents_sulci',
                                                         'extents_fundi',
                                                         'extents_midaxis',
                                                         'extents_patches',
                                                         'extents_regions',
                                                         'curvatures_sulci',
                                                         'curvatures_fundi',
                                                         'curvatures_pits',
                                                         'curvatures_midaxis',
                                                         'curvatures_patches',
                                                         'curvatures_regions',
                                                         'depths_sulci',
                                                         'depths_fundi',
                                                         'depths_pits',
                                                         'depths_midaxis',
                                                         'depths_patches',
                                                         'depths_regions',
                                                         'spectra_sulci',
                                                         'spectra_fundi',
                                                         'spectra_midaxis',
                                                         'spectra_patches',
                                                         'spectra_regions'],
                                          output_names=['measures'],
                                          function = write_measures_to_database),
                            name='Write_measures_to_database')

measures_table = pe.Node(util.Function(input_names = ['measures'],
                                       output_names=['success'],
                                       function = write_measures_to_table),
                         name='Write_measures_to_table')

# Connect surface maps to database nodes
flo2.connect([(surfaces, maps_database, [('depth_file','depth_file'),
                                ('mean_curv_file','mean_curv_file'),
                                ('gauss_curv_file','gauss_curv_file')])])

# Connect feature to database nodes
flo2.connect([(sulcus_segmentation, features_database, [('segmented_sulci', 'segmented_sulci')]),
              (fundus_segmentation, features_database, [('segmented_fundi', 'segmented_fundi')]),
              (pit_extraction, features_database, [('pits', 'pits')]),
              (midaxis_segmentation, features_database, 
                          [('segmented_midaxis', 'segmented_midaxis')])])

# Connect feature measures to database nodes
flo2.connect([(positions, measures_database, [('positions_sulci', 'positions_sulci'),
                                              ('positions_fundi', 'positions_fundi'),
                                              ('positions_pits', 'positions_pits'),
                                              ('positions_midaxis', 'positions_midaxis')]),
              (extents, measures_database, [('extents_sulci', 'extents_sulci'),
                                            ('extents_fundi', 'extents_fundi'),
                                            ('extents_midaxis', 'extents_midaxis')]),
              (curvatures, measures_database, [('curvatures_sulci', 'curvatures_sulci'),
                                               ('curvatures_fundi', 'curvatures_fundi'),
                                               ('curvatures_pits', 'curvatures_pits'),
                                               ('curvatures_midaxis', 'curvatures_midaxis')]),
              (depths, measures_database, [('depths_sulci', 'depths_sulci'),
                                           ('depths_fundi', 'depths_fundi'),
                                           ('depths_pits', 'depths_pits'),
                                           ('depths_midaxis', 'depths_midaxis')]),
              (spectra, measures_database, [('spectra_sulci', 'spectra_sulci'),
                                            ('spectra_fundi', 'spectra_fundi'),
                                            ('spectra_midaxis', 'spectra_midaxis')])])

# Connect label measures to database nodes
flo2.connect([(positions, measures_database, [('positions_patches', 'positions_patches'),
                                              ('positions_regions', 'positions_regions')]),
              (extents, measures_database, [('extents_patches', 'extents_patches'),
                                            ('extents_regions', 'extents_regions')]),
              (curvatures, measures_database, [('curvatures_patches', 'curvatures_patches'),
                                               ('curvatures_regions', 'curvatures_regions')]),
              (depths, measures_database, [('depths_patches', 'depths_patches'),
                                           ('depths_regions', 'depths_regions')]),
              (spectra, measures_database, [('spectra_patches', 'spectra_patches'),
                                            ('spectra_regions', 'spectra_regions')])])

# Connect measure to table nodes
flo2.connect([(measures_database, measures_table, [('measures', 'measures')])])

"""

##############################################################################
#    Run workflow
##############################################################################
if __name__== '__main__':

    mbflow.write_graph(graph2use='flat')
    mbflow.write_graph(graph2use='hierarchical')
    #flo1.write_graph(graph2use='flat')
    #flo1.write_graph(graph2use='hierarchical')
    mbflow.run(updatehash=False)  #mbflow.run(updatehash=True)

