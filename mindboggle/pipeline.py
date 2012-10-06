#!/usr/bin/python
"""
This is Mindboggle's nipype software pipeline!

Examples
--------
python mindboggle.py <output path> <1 or more subject names>
python mindboggle.py output HLN-12-1 HLN-12-2

.. note::
  Mindboggle assumes a file tree like FreeSurfer's,
  and for label initialization, assumes that subjects have been processed
  by FreeSurfer (autorecon -all), so subject names correspond to directory
  names in FreeSurfer's subjects directory.

For more information about Mindboggle,
see the website: http://www.mindboggle.info and
read the documentation: http://mindboggle.info/software/documentation.html

For information on Nipype (http://www.nipy.org/nipype/):
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3159964/


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

################################################################################
#
#   Mindboggle workflow:
#   * Multi-atlas labeling
#   * Feature extraction
#   * Feature-based labeling
#   * Shape measurement
#
#   Followed by:
#
#   Label Volume workflow:
#   * Volume-filling labels
#   * Label evaluation
#
################################################################################

#===============================================================================
#  Command line arguments
#===============================================================================
import sys, os

args = sys.argv[:]
if len(args) < 3:
    print("\n\t Please provide the names of an output directory\n" +
          " \t and one or more subjects corresponding to the names\n" +
          " \t of directories within FreeSurfer's subjects directory.\n")
    print("\t Example: python " + args[0] + " output HLN-12-1 HLN-12-2\n")
    sys.exit()
else:
    output_path = str(args[1])
    subjects = list(args[2::])

#===============================================================================
#  User settings
#===============================================================================
input_vtk = False  # Load my VTK surfaces directly (not FreeSurfer surfaces)
fill_volume = True  # Fill (gray matter) volumes with surface labels
include_thickness = True  # Include FreeSurfer's thickness measure
include_convexity = True  # Include FreeSurfer's convexity measure (sulc.pial)
#-------------------------------------------------------------------------------
# Labeling protocol used by Mindboggle:
# 'DKT31': 'Desikan-Killiany-Tourville (DKT) protocol with 31 labeled regions
# 'DKT25': 'fundus-friendly' version of the DKT protocol following fundi
#-------------------------------------------------------------------------------
protocol = 'DKT25'
#-------------------------------------------------------------------------------
# Initialize labels with:
# 'DKatlas': the standard FreeSurfer classifier atlas trained on the DK protocol
# 'DKTatlas': a FreeSurfer-style classifier atlas trained on the DKT protocol
# 'max': maximum probability (majority vote) labels from multiple atlases
# 'manual': process manual labels (atlas)
#-------------------------------------------------------------------------------
init_labels = 'manual'
#-------------------------------------------------------------------------------
# Labeling source:
# 'manual': manual edits
# FUTURE:
# <'adjusted': manual edits after automated alignment to fundi>
#-------------------------------------------------------------------------------
label_method = 'manual'
hemis = ['lh','rh']  # Prepend ('lh.'/'rh.') indicating left/right surfaces
#-------------------------------------------------------------------------------
# Evaluation options
#-------------------------------------------------------------------------------
evaluate_surface_labels = 0 #False  # Surface overlap: auto vs. manual labels
evaluate_volume_labels = 0 #False  # Volume overlap: auto vs. manual labels
run_atlasflow = 0#True
run_measureflow = True
run_featureflow = 0#True
run_shapeflow = 0#True

#===============================================================================
#  Setup: import libraries, set file paths, and initialize main workflow
#===============================================================================
#-------------------------------------------------------------------------------
# Import system and nipype Python libraries
#-------------------------------------------------------------------------------
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import Function as Fn
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import DataGrabber, DataSink
#-------------------------------------------------------------------------------
# Import Mindboggle Python libraries
#-------------------------------------------------------------------------------
from mindboggle.utils.io_vtk import rewrite_scalars, write_mean_shapes_table, \
     load_scalar, freesurface_to_vtk, freecurvature_to_vtk, freeannot_to_vtk, \
     vtk_to_freelabels
from mindboggle.utils.io_file import read_columns
from mindboggle.utils.io_free import labels_to_annot, labels_to_volume
from mindboggle.utils.mesh_operations import find_neighbors
from mindboggle.label.multiatlas_labeling import register_template,\
     transform_atlas_labels, majority_vote_label
from mindboggle.label.relabel import relabel_volume
from mindboggle.label.label_functions import label_with_classifier
from mindboggle.measure.measure_functions import compute_area, compute_depth, \
     compute_curvature
from mindboggle.extract.extract_folds import extract_folds
from mindboggle.extract.extract_fundi import extract_fundi
from mindboggle.evaluate.evaluate_labels import measure_surface_overlap, \
     measure_volume_overlap
#from mindboggle import get_info
#-------------------------------------------------------------------------------
# Paths
#-------------------------------------------------------------------------------
subjects_path = os.environ['SUBJECTS_DIR']  # FreeSurfer subjects directory
data_path = os.environ['MINDBOGGLE_DATA']  # Mindboggle data directory
temp_path = os.path.join(output_path, 'workspace')  # Where to save temp files
ccode_path = os.environ['MINDBOGGLE_TOOLS']
#info_path = os.path.join(get_info()['pkg_path'], 'info')
info_path = os.path.join(os.environ['MINDBOGGLE'], 'info')
atlases_path = subjects_path
# Label with classifier atlas
templates_path = os.path.join(subjects_path, 'MindboggleTemplates')
# Label with classifier atlas
classifier_path = os.path.join(subjects_path, 'MindboggleClassifierAtlases')
classifier_atlas = 'DKTatlas40.gcs'
#-------------------------------------------------------------------------------
# Initialize main workflow
#-------------------------------------------------------------------------------
mbflow = Workflow(name='Mindboggle')
mbflow.base_dir = temp_path
if not os.path.isdir(temp_path):  os.makedirs(temp_path)

#===============================================================================
#   Inputs and outputs
#===============================================================================
#-------------------------------------------------------------------------------
# Iterate inputs over subjects, hemispheres
# (surfaces are assumed to take the form: lh.pial or lh.pial.vtk)
#-------------------------------------------------------------------------------
info = Node(name = 'Inputs',
            interface = IdentityInterface(fields=['subject', 'hemi']))
info.iterables = ([('subject', subjects), ('hemi', hemis)])
#-------------------------------------------------------------------------------
# Location and structure of the surface inputs
#-------------------------------------------------------------------------------
surf = Node(name = 'Surfaces',
            interface = DataGrabber(infields=['subject', 'hemi'],
                                    outfields=['surface_files', 'sphere_files']))
surf.inputs.base_directory = subjects_path
surf.inputs.template = '%s/surf/%s.%s'
surf.inputs.template_args['surface_files'] = [['subject', 'hemi', 'pial']]
surf.inputs.template_args['sphere_files'] = [['subject', 'hemi', 'sphere']]
if include_thickness:
    surf.inputs.template_args['thickness_files'] = [['subject', 'hemi', 'thickness']]
if include_convexity:
    surf.inputs.template_args['convexity_files'] = [['subject', 'hemi', 'sulc']]
mbflow.connect([(info, surf, [('subject','subject'), ('hemi','hemi')])])
#-------------------------------------------------------------------------------
# Outputs
#-------------------------------------------------------------------------------
sink = Node(DataSink(), name = 'Results')
sink.inputs.base_directory = output_path
sink.inputs.container = 'results'
if not os.path.isdir(output_path):  os.makedirs(output_path)
#-------------------------------------------------------------------------------
# Convert surfaces to VTK
#-------------------------------------------------------------------------------
if not input_vtk:
    convertsurf = Node(name = 'Surf_to_VTK',
                       interface = Fn(function = freesurface_to_vtk,
                                      input_names = ['surface_file'],
                                      output_names = ['vtk_file']))
    mbflow.connect([(surf, convertsurf, [('surface_files','surface_file')])])
#-------------------------------------------------------------------------------
# Evaluation inputs: location and structure of atlas surfaces
#-------------------------------------------------------------------------------
if evaluate_surface_labels or init_labels == 'manual':
    atlas = Node(name = 'Atlases',
                 interface = DataGrabber(infields=['subject','hemi'],
                                         outfields=['atlas_file']))
    atlas.inputs.base_directory = atlases_path

    atlas.inputs.template = '%s/label/%s.labels.' +\
                            protocol + '.' + label_method + '.vtk'
    atlas.inputs.template_args['atlas_file'] = [['subject','hemi']]

    mbflow.connect([(info, atlas, [('subject','subject'),('hemi','hemi')])])
#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------
ctx_labels_file = os.path.join(info_path, 'labels.surface.' + protocol + '.txt')
ctx_label_numbers, ctx_label_names, RGBs = read_columns(ctx_labels_file,
                                                n_columns=3, trail=True)

################################################################################
#
#   Multi-atlas labeling workflow
#
################################################################################
if run_atlasflow:

    atlasflow = Workflow(name='Label_initialization')

    #===========================================================================
    #   Initialize labels with FreeSurfer's standard DK classifier atlas
    #===========================================================================
    if init_labels == 'DKatlas':
        freelabels = Node(name = 'DK_annot_to_VTK',
                        interface = Fn(function = freeannot_to_vtk,
                                       input_names = ['surface_file',
                                                      'hemi',
                                                      'subject',
                                                      'subjects_path',
                                                      'annot_name'],
                                       output_names = ['labels',
                                                       'vtk_file']))
        atlasflow.add_nodes([freelabels])
        if input_vtk:
            mbflow.connect([(surf, atlasflow,
                             [('surface_files',
                               'DK_annot_to_VTK.surface_file')])])
        else:
            mbflow.connect([(convertsurf, atlasflow,
                             [('vtk_file',
                               'DK_annot_to_VTK.surface_file')])])
        mbflow.connect([(info, atlasflow,
                         [('hemi', 'DK_annot_to_VTK.hemi'),
                          ('subject', 'DK_annot_to_VTK.subject')])])
        freelabels.inputs.subjects_path = subjects_path
        freelabels.inputs.annot_name = 'aparc.annot'
    #===========================================================================
    #   Initialize labels with the DKT classifier atlas
    #===========================================================================
    elif init_labels == 'DKTatlas':
        """
        Label a brain with the DKT atlas using FreeSurfer's mris_ca_label.
        """
        classifier = Node(name = 'Label_with_DKTatlas',
                          interface = Fn(function = label_with_classifier,
                                         input_names = ['hemi',
                                                        'subject',
                                                        'subjects_path',
                                                        'sphere_file',
                                                        'classifier_path',
                                                        'classifier_atlas'],
                                         output_names = ['annot_name',
                                                         'annot_file']))
        atlasflow.add_nodes([classifier])
        mbflow.connect([(info, atlasflow,
                         [('hemi', 'Label_with_DKTatlas.hemi'),
                          ('subject', 'Label_with_DKTatlas.subject')])])
        classifier.inputs.subjects_path = subjects_path
        mbflow.connect([(surf, atlasflow,
                         [('sphere_files',
                           'Label_with_DKTatlas.sphere_file')])])
        classifier.inputs.classifier_path = classifier_path
        classifier.inputs.classifier_atlas = classifier_atlas

        # Convert .annot file to .vtk format
        classifier2vtk = Node(name = 'DKT_annot_to_VTK',
                              interface = Fn(function = freeannot_to_vtk,
                                             input_names = ['surface_file',
                                                            'hemi',
                                                            'subject',
                                                            'subjects_path',
                                                            'annot_name'],
                                             output_names = ['labels',
                                                             'vtk_file']))
        atlasflow.add_nodes([classifier2vtk])
        if input_vtk:
            mbflow.connect([(surf, atlasflow,
                             [('surface_files',
                               'DKT_annot_to_VTK.surface_file')])])
        else:
            mbflow.connect([(convertsurf, atlasflow,
                             [('vtk_file',
                               'DKT_annot_to_VTK.surface_file')])])
        mbflow.connect([(info, atlasflow,
                         [('hemi', 'DKT_annot_to_VTK.hemi'),
                          ('subject', 'DKT_annot_to_VTK.subject')])])
        classifier2vtk.inputs.subjects_path = subjects_path
        atlasflow.connect([(classifier, classifier2vtk,
                            [('annot_name', 'annot_name')])])
    #===========================================================================
    #   Initialize labels using multi-atlas registration
    #===========================================================================
    elif init_labels == 'max':
        #-----------------------------------------------------------------------
        # Register surfaces to average template
        #-----------------------------------------------------------------------
        free_template = 'OASIS-TRT-20'  # FreeSurfer template

        register = Node(name = 'Register_template',
                        interface = Fn(function = register_template,
                                       input_names = ['hemi',
                                                      'sphere_file',
                                                      'transform',
                                                      'templates_path',
                                                      'template'],
                                       output_names = ['transform']))
        atlasflow.add_nodes([register])
        mbflow.connect([(info, atlasflow, [('hemi', 'Register_template.hemi')]),
                        (surf, atlasflow, [('sphere_files',
                                            'Register_template.sphere_file')])])
        register.inputs.transform = 'sphere_to_' + template + '_template.reg'
        register.inputs.templates_path = os.path.join(templates_path, 'freesurfer')
        register.inputs.template = free_template + '.tif'
        #-----------------------------------------------------------------------
        # Register atlases to subject via template
        #-----------------------------------------------------------------------
        # Load atlas list
        atlas_list_file = os.path.join(info_path, 'atlases.txt')
        atlas_list = read_columns(atlas_list_file, 1)[0]

        transform = MapNode(name = 'Transform_labels',
                            iterfield = ['atlas'],
                            interface = Fn(function = transform_atlas_labels,
                                           input_names = ['hemi',
                                                          'subject',
                                                          'transform',
                                                          'subjects_path',
                                                          'atlas',
                                                          'atlas_string'],
                                           output_names = ['output_file']))
        atlasflow.add_nodes([transform])
        mbflow.connect([(info, atlasflow,
                         [('hemi', 'Transform_labels.hemi'),
                          ('subject', 'Transform_labels.subject')])])
        atlasflow.connect([(register, transform, [('transform', 'transform')])])
        #transform.inputs.transform = 'sphere_to_' + template + '_template.reg'
        transform.inputs.subjects_path = subjects_path
        transform.inputs.atlas = atlas_list
        transform.inputs.atlas_string = 'labels.' + protocol + '.' + label_method
        #-----------------------------------------------------------------------
        # Majority vote label
        #-----------------------------------------------------------------------
        vote = Node(name='Label_vote',
                    interface = Fn(function = majority_vote_label,
                                   input_names = ['surface_file',
                                                  'annot_files'],
                                   output_names = ['labels_max',
                                                   'label_counts',
                                                   'label_votes',
                                                   'consensus_vertices',
                                                   'maxlabel_file',
                                                   'labelcounts_file',
                                                   'labelvotes_file']))
        atlasflow.add_nodes([vote])
        if input_vtk:
            mbflow.connect([(surf, atlasflow,
                             [('surface_files', 'Label_vote.surface_file')])])
        else:
            mbflow.connect([(convertsurf, atlasflow,
                             [('vtk_file', 'Label_vote.surface_file')])])
        atlasflow.connect([(transform, vote, [('output_file', 'annot_files')])])
        mbflow.connect([(atlasflow, sink,
                         [('Label_vote.maxlabel_file', 'labels.@max'),
                          ('Label_vote.labelcounts_file', 'labels.@counts'),
                          ('Label_vote.labelvotes_file', 'labels.@votes')])])
    #===========================================================================
    #   Skip label initialization and process manual (atlas) labels
    #===========================================================================
    elif init_labels == 'manual':
        atlaslabels = Node(name = 'Atlas_labels',
                           interface = Fn(function = load_scalar,
                                          input_names = ['filename',
                                                         'return_arrays'],
                                          output_names = ['points',
                                                          'faces',
                                                          'scalars',
                                                          'n_vertices']))
        atlasflow.add_nodes([atlaslabels])
        mbflow.connect([(atlas, atlasflow,
                         [('atlas_file', 'Atlas_labels.filename')])])
        atlaslabels.inputs.return_arrays = 0  # 0: return lists instead of arrays


################################################################################
#
#   Surface measurement workflow
#
################################################################################
if run_measureflow:

    measureflow = Workflow(name='Surface_measurement')

    #===========================================================================
    #   Surface measurements
    #===========================================================================
    #---------------------------------------------------------------------------
    # Measure surface area
    #---------------------------------------------------------------------------
    area = Node(name='Area',
                interface = Fn(function = compute_area,
                               input_names = ['command',
                                              'surface_file'],
                               output_names = ['area_file']))
    area_command = os.path.join(ccode_path, 'area', 'PointAreaMain')
    area.inputs.command = area_command
    #---------------------------------------------------------------------------
    # Measure surface depth
    #---------------------------------------------------------------------------
    depth = Node(name='Depth',
                 interface = Fn(function = compute_depth,
                                input_names = ['command',
                                               'surface_file'],
                                output_names = ['depth_file']))
    depth_command = os.path.join(ccode_path, 'travel_depth', 'TravelDepthMain')
    depth.inputs.command = depth_command
    #---------------------------------------------------------------------------
    # Measure surface curvature
    #---------------------------------------------------------------------------
    curvature = Node(name='Curvature',
                     interface = Fn(function = compute_curvature,
                                    input_names = ['command',
                                                   'surface_file'],
                                    output_names = ['mean_curvature_file',
                                                    'gauss_curvature_file',
                                                    'max_curvature_file',
                                                    'min_curvature_file',
                                                    'min_curvature_vector_file']))
    curvature_command = os.path.join(ccode_path, 'curvature', 'CurvatureMain')
    curvature.inputs.command = curvature_command
    #---------------------------------------------------------------------------
    # Convert FreeSurfer surface measures to VTK
    #---------------------------------------------------------------------------
    if include_thickness:
        convertthickness = Node(name = 'Thickness_to_VTK',
                           interface = Fn(function = freecurvature_to_vtk,
                                          input_names = ['file_string',
                                                         'surface_file',
                                                         'hemi',
                                                         'subject',
                                                         'subjects_path'],
                                          output_names = ['vtk_file']))
        measureflow.add_nodes([convertthickness])
        convertthickness.inputs.file_string = 'thickness'
        if not input_vtk:
            mbflow.connect([(convertsurf, measureflow,
                             [('vtk_file','Thickness_to_VTK.surface_file')])])
        else:
            mbflow.connect([(surf, measureflow,
                             [('surface_files','Thickness_to_VTK.surface_file')])])
        mbflow.connect([(info, measureflow,
                         [('hemi','Thickness_to_VTK.hemi'),
                          ('subject','Thickness_to_VTK.subject')])])
        convertthickness.inputs.subjects_path = subjects_path
        mbflow.connect([(measureflow, sink,
                         [('Thickness_to_VTK.vtk_file', 'measures.@thickness')])])
    if include_convexity:
        convertconvexity = Node(name = 'Convexity_to_VTK',
                           interface = Fn(function = freecurvature_to_vtk,
                                          input_names = ['file_string',
                                                         'surface_file',
                                                         'hemi',
                                                         'subject',
                                                         'subjects_path'],
                                          output_names = ['vtk_file']))
        measureflow.add_nodes([convertconvexity])
        convertconvexity.inputs.file_string = 'sulc'
        if not input_vtk:
            mbflow.connect([(convertsurf, measureflow,
                             [('vtk_file','Convexity_to_VTK.surface_file')])])
        else:
            mbflow.connect([(surf, measureflow,
                             [('surface_files','Convexity_to_VTK.surface_file')])])
        mbflow.connect([(info, measureflow, [('hemi','Convexity_to_VTK.hemi'),
                        ('subject','Convexity_to_VTK.subject')])])
        convertconvexity.inputs.subjects_path = subjects_path
        mbflow.connect([(measureflow, sink,
                         [('Convexity_to_VTK.vtk_file', 'measures.@convexity')])])
    #---------------------------------------------------------------------------
    # Add and connect nodes, save output files
    #---------------------------------------------------------------------------
    measureflow.add_nodes([area, depth, curvature])
    if input_vtk:
        mbflow.connect([(surf, measureflow,
                         [('surface_files','Area.surface_file')])])
        mbflow.connect([(surf, measureflow,
                         [('surface_files','Depth.surface_file')])])
        mbflow.connect([(surf, measureflow,
                         [('surface_files','Curvature.surface_file')])])
    else:
        mbflow.connect([(convertsurf, measureflow,
                         [('vtk_file', 'Area.surface_file')])])
        mbflow.connect([(convertsurf, measureflow,
                         [('vtk_file', 'Depth.surface_file')])])
        mbflow.connect([(convertsurf, measureflow,
                         [('vtk_file', 'Curvature.surface_file')])])
    mbflow.connect([(measureflow, sink,
                     [('Area.area_file', 'measures.@area')])])
    mbflow.connect([(measureflow, sink,
                     [('Depth.depth_file', 'measures.@depth')])])
    mbflow.connect([(measureflow, sink,
                     [('Curvature.mean_curvature_file',
                       'measures.@mean_curvature'),
                      ('Curvature.gauss_curvature_file',
                       'measures.@gauss_curvature'),
                      ('Curvature.max_curvature_file',
                       'measures.@max_curvature'),
                      ('Curvature.min_curvature_file',
                       'measures.@min_curvature'),
                      ('Curvature.min_curvature_vector_file',
                       'measures.@min_curvature_vectors')])])

################################################################################
#
#   Feature extraction workflow
#
################################################################################
if run_featureflow:

    featureflow = Workflow(name='Feature_extraction')

    #===========================================================================
    #   Feature extraction
    #===========================================================================
    #---------------------------------------------------------------------------
    # Find all neighbors
    #---------------------------------------------------------------------------
    load_surface = Node(name = 'Load_surface',
                        interface = Fn(function = load_scalar,
                                       input_names = ['filename',
                                                      'return_arrays'],
                                       output_names = ['points',
                                                       'faces',
                                                       'scalars',
                                                       'n_vertices']))
    featureflow.add_nodes([load_surface])
    if input_vtk:
        mbflow.connect([(surf, featureflow,
                         [('surface_files','Load_surface.filename')])])
    else:
        mbflow.connect([(convertsurf, featureflow,
                         [('vtk_file', 'Load_surface.filename')])])
    load_surface.inputs.return_arrays = 0  # 0: return lists instead of arrays

    neighbors = Node(name='Neighbors',
                     interface = Fn(function = find_neighbors,
                                    input_names = ['faces', 'n_vertices'],
                                    output_names = ['neighbor_lists']))
    featureflow.add_nodes([neighbors])
    featureflow.connect([(load_surface, neighbors,
                          [('faces','faces'), ('n_vertices','n_vertices')])])
    #---------------------------------------------------------------------------
    # Extract folds
    #---------------------------------------------------------------------------
    fraction_folds = 0.5
    min_fold_size = 50

    folds = Node(name='Folds',
                 interface = Fn(function = extract_folds,
                                input_names = ['depth_file',
                                               'neighbor_lists',
                                               'fraction_folds',
                                               'min_fold_size'],
                                output_names = ['folds',
                                                'n_folds']))
    featureflow.add_nodes([folds])
    mbflow.connect([(measureflow, featureflow,
                     [('Depth.depth_file','Folds.depth_file')])])
    featureflow.connect([(neighbors, folds,
                          [('neighbor_lists','neighbor_lists')])])
    folds.inputs.fraction_folds = fraction_folds
    folds.inputs.min_fold_size = min_fold_size
    #---------------------------------------------------------------------------
    # Extract sulci from folds
    #---------------------------------------------------------------------------
    """
    sulci = Node(name='Sulci',
                 interface = Fn(function = rewrite_scalars,
                                input_names = ['folds',
                                               'sulci_list'],
                                output_names = ['sulci']))
    featureflow.add_nodes([sulci])
    mbflow.connect([(folds, sulci,
                     [('folds','folds')])])
    sulci.inputs.output_vtk = 'folds.vtk'
    featureflow.connect([(folds, save_folds, [('folds','new_scalars')])])
    featureflow.connect([(folds, save_folds, [('folds','filter_scalars')])])
    mbflow.connect([(featureflow, sink,
                     [('Save_folds.output_vtk','features.@folds')])])
    """
    #---------------------------------------------------------------------------
    # Extract fundi (curves at the bottoms of sulci)
    #---------------------------------------------------------------------------
    thr = 0.5
    min_distance = 5.0

    fundi = Node(name='Fundi',
                 interface = Fn(function = extract_fundi,
                                input_names = ['folds',
                                               'n_folds',
                                               'neighbor_lists',
                                               'depth_file',
                                               'mean_curvature_file',
                                               'min_curvature_vector_file',
                                               'min_fold_size',
                                               'min_distance',
                                               'thr'],
                                output_names = ['fundi']))
    featureflow.connect([(folds, fundi, [('folds','folds'),
                                         ('n_folds','n_folds')]),
                         (neighbors, fundi,
                          [('neighbor_lists','neighbor_lists')])])
    mbflow.connect([(measureflow, featureflow,
                     [('Depth.depth_file','Fundi.depth_file'),
                      ('Curvature.mean_curvature_file',
                       'Fundi.mean_curvature_file'),
                      ('Curvature.min_curvature_vector_file',
                       'Fundi.min_curvature_vector_file')])])
    fundi.inputs.min_fold_size = min_fold_size
    fundi.inputs.min_distance = min_distance
    fundi.inputs.thr = thr
    #---------------------------------------------------------------------------
    # Write sulci and fundi to VTK files
    #---------------------------------------------------------------------------
    save_sulci = Node(name='Save_sulci',
                      interface = Fn(function = rewrite_scalars,
                                     input_names = ['input_vtk',
                                                    'output_vtk',
                                                    'new_scalars',
                                                    'filter_scalars'],
                                     output_names = ['output_vtk']))
    featureflow.add_nodes([save_sulci])
    mbflow.connect([(measureflow, featureflow,
                     [('Depth.depth_file','Save_sulci.input_vtk')])])
    save_sulci.inputs.output_vtk = 'sulci.vtk'
    featureflow.connect([(folds, save_sulci, [('folds','new_scalars')])])
    featureflow.connect([(folds, save_sulci, [('folds','filter_scalars')])])
    mbflow.connect([(featureflow, sink,
                     [('Save_sulci.output_vtk','features.@sulci')])])

    save_fundi = save_sulci.clone('Save_fundi')
    featureflow.add_nodes([save_fundi])
    mbflow.connect([(measureflow, featureflow,
                     [('Depth.depth_file','Save_fundi.input_vtk')])])
    save_fundi.inputs.output_vtk = 'fundi.vtk'
    featureflow.connect([(fundi, save_fundi, [('fundi','new_scalars')])])
    featureflow.connect([(folds, save_fundi, [('folds','filter_scalars')])])
    mbflow.connect([(featureflow, sink,
                     [('Save_fundi.output_vtk','features.@fundi')])])

################################################################################
#
#   Shape analysis workflow
#
################################################################################
if run_shapeflow:

    shapeflow = Workflow(name='Shape_analysis')
    column_names = ['labels', 'area', 'depth', 'mean_curvature',
                    'gauss_curvature', 'max_curvature', 'min_curvature']
    if include_thickness:
        column_names.append('thickness')
    if include_convexity:
        column_names.append('convexity')
    shape_files = [x + '_file' for x in column_names[2::]]
    input_names = ['filename', 'column_names', 'labels', 'area_file']
    input_names.extend(shape_files)
    #===========================================================================
    # Labeled surface patch shapes
    #===========================================================================
    labeltable = Node(name='Label_table',
                      interface = Fn(function = write_mean_shapes_table,
                                     input_names = input_names,
                                     output_names = ['means_file',
                                                     'norm_means_file']))
    shapeflow.add_nodes([labeltable])
    labeltable.inputs.filename = 'label_shapes.txt'
    labeltable.inputs.column_names = column_names
    mbflow.connect([(measureflow, shapeflow,
                     [('Area.area_file','Label_table.area_file')])])
    mbflow.connect([(measureflow, shapeflow,
                     [('Depth.depth_file','Label_table.depth_file')])])
    mbflow.connect([(measureflow, shapeflow,
                     [('Curvature.mean_curvature_file',
                       'Label_table.mean_curvature_file')])])
    mbflow.connect([(measureflow, shapeflow,
                     [('Curvature.gauss_curvature_file',
                       'Label_table.gauss_curvature_file')])])
    mbflow.connect([(measureflow, shapeflow,
                     [('Curvature.max_curvature_file',
                       'Label_table.max_curvature_file')])])
    mbflow.connect([(measureflow, shapeflow,
                     [('Curvature.min_curvature_file',
                       'Label_table.min_curvature_file')])])
    if include_thickness:
        mbflow.connect([(measureflow, shapeflow,
                         [('Thickness_to_VTK.vtk_file', 'Label_table.thickness_file')])])
    if include_convexity:
        mbflow.connect([(measureflow, shapeflow,
                         [('Convexity_to_VTK.vtk_file', 'Label_table.convexity_file')])])
    #---------------------------------------------------------------------------
    # Use initial labels assigned by classifier atlas (ex: FreeSurfer labels)
    #---------------------------------------------------------------------------
    if init_labels == 'DKatlas':
        mbflow.connect([(atlasflow, shapeflow,
                         [('DK_annot_to_VTK.labels','Label_table.labels')])])
    #---------------------------------------------------------------------------
    # Use initial labels assigned by multi-atlas registration
    #---------------------------------------------------------------------------
    elif init_labels == 'DKTatlas':
        mbflow.connect([(atlasflow, shapeflow,
                         [('DKT_annot_to_VTK.labels','Label_table.labels')])])
    #---------------------------------------------------------------------------
    # Use initial labels assigned by multi-atlas registration
    #---------------------------------------------------------------------------
    elif init_labels == 'max':
        mbflow.connect([(atlasflow, shapeflow,
                         [('Label_vote.labels_max','Label_table.labels')])])
    #---------------------------------------------------------------------------
    # Use manual (atlas) labels
    #---------------------------------------------------------------------------
    elif init_labels == 'manual':
        mbflow.connect([(atlasflow, shapeflow,
                         [('Atlas_labels.scalars','Label_table.labels')])])
    # Save results
    mbflow.connect([(shapeflow, sink,
                     [('Label_table.means_file', 'shapes.@labels'),
                      ('Label_table.norm_means_file', 'shapes.@labels_norm')])])
    #===========================================================================
    # Sulcus fold shapes
    #===========================================================================
    if run_featureflow:
        foldtable = labeltable.clone('Fold_table')
        shapeflow.add_nodes([foldtable])
        foldtable.inputs.filename = 'fold_shapes.txt'
        foldtable.inputs.column_names = column_names
        mbflow.connect([(measureflow, shapeflow,
                         [('Area.area_file','Fold_table.area_file')])])
        mbflow.connect([(featureflow, shapeflow,
                         [('Folds.folds','Fold_table.labels')])])
        mbflow.connect([(measureflow, shapeflow,
                         [('Depth.depth_file','Fold_table.depth_file')])])
        mbflow.connect([(measureflow, shapeflow,
                         [('Curvature.mean_curvature_file',
                           'Fold_table.mean_curvature_file')])])
        mbflow.connect([(measureflow, shapeflow,
                         [('Curvature.gauss_curvature_file',
                           'Fold_table.gauss_curvature_file')])])
        mbflow.connect([(measureflow, shapeflow,
                         [('Curvature.max_curvature_file',
                           'Fold_table.max_curvature_file')])])
        mbflow.connect([(measureflow, shapeflow,
                         [('Curvature.min_curvature_file',
                           'Fold_table.min_curvature_file')])])
        if include_thickness:
            mbflow.connect([(measureflow, shapeflow,
                             [('Thickness_to_VTK.vtk_file',
                               'Fold_table.thickness_file')])])
        if include_convexity:
            mbflow.connect([(measureflow, shapeflow,
                             [('Convexity_to_VTK.vtk_file',
                               'Fold_table.convexity_file')])])
        # Save results
        mbflow.connect([(shapeflow, sink,
                         [('Fold_table.means_file', 'shapes.@folds'),
                          ('Fold_table.norm_means_file', 'shapes.@folds_norm')])])
    #===========================================================================
    # Fundus shapes
    #===========================================================================
    if run_featureflow:
        fundustable = labeltable.clone('Fundus_table')
        shapeflow.add_nodes([fundustable])
        fundustable.inputs.filename = 'fundus_shapes.txt'
        fundustable.inputs.column_names = column_names
        mbflow.connect([(measureflow, shapeflow,
                         [('Area.area_file','Fundus_table.area_file')])])
        mbflow.connect([(featureflow, shapeflow,
                         [('Fundi.fundi','Fundus_table.labels')])])
        mbflow.connect([(measureflow, shapeflow,
                         [('Depth.depth_file','Fundus_table.depth_file')])])
        mbflow.connect([(measureflow, shapeflow,
                         [('Curvature.mean_curvature_file',
                           'Fundus_table.mean_curvature_file')])])
        mbflow.connect([(measureflow, shapeflow,
                         [('Curvature.gauss_curvature_file',
                           'Fundus_table.gauss_curvature_file')])])
        mbflow.connect([(measureflow, shapeflow,
                         [('Curvature.max_curvature_file',
                           'Fundus_table.max_curvature_file')])])
        mbflow.connect([(measureflow, shapeflow,
                         [('Curvature.min_curvature_file',
                           'Fundus_table.min_curvature_file')])])
        if include_thickness:
            mbflow.connect([(measureflow, shapeflow,
                             [('Thickness_to_VTK.vtk_file',
                               'Fundus_table.thickness_file')])])
        if include_convexity:
            mbflow.connect([(measureflow, shapeflow,
                             [('Convexity_to_VTK.vtk_file',
                               'Fundus_table.convexity_file')])])
        # Save results
        mbflow.connect([(shapeflow, sink,
                         [('Fundus_table.means_file', 'shapes.@fundi'),
                          ('Fundus_table.norm_means_file', 'shapes.@fundi_norm')])])

################################################################################
#
#   Surface label evaluation
#
################################################################################
if evaluate_surface_labels:

    #---------------------------------------------------------------------------
    # Evaluate surface labels
    #---------------------------------------------------------------------------
    eval_surf_labels = Node(name='Evaluate_surface_labels',
                            interface = Fn(function = measure_surface_overlap,
                                           input_names = ['command',
                                                          'labels_file1',
                                                          'labels_file2'],
                                           output_names = ['overlaps']))
    mbflow.add_nodes([eval_surf_labels])
    surface_overlap_command = os.path.join(ccode_path,
        'surface_overlap', 'SurfaceOverlapMain')
    eval_surf_labels.inputs.command = surface_overlap_command
    mbflow.connect([(atlas, eval_surf_labels, [('atlas_file','labels_file1')])])
    if init_labels == 'DKatlas':
        mbflow.connect([(atlasflow, eval_surf_labels,
                         [('DK_annot_to_VTK.vtk_file','labels_file2')])])
    elif init_labels == 'DKTatlas':
        mbflow.connect([(atlasflow, eval_surf_labels,
                         [('DKT_annot_to_VTK.vtk_file','labels_file2')])])
    elif init_labels == 'max':
        mbflow.connect([(atlasflow, eval_surf_labels,
                         [('Label_vote.maxlabel_file','labels_file2')])])
    elif init_labels == 'manual':
        mbflow.connect([(atlas, eval_surf_labels,
                         [('atlas_file','labels_file2')])])
    #mbflow.connect([(eval_surf_labels, sink,
    #                 [('overlaps', 'evaluate_labels')])])

################################################################################
#
#   Fill volume prep workflow:
#   Convert labels from VTK to .annot format
#
################################################################################
if fill_volume:

    annotflow = Workflow(name='Fill_volume_prep')

    #===========================================================================
    #   Convert VTK labels to .annot format
    #===========================================================================
    #---------------------------------------------------------------------------
    # Write .label files for surface vertices
    #---------------------------------------------------------------------------
    writelabels = Node(name='Write_label_files',
                       interface = Fn(function = vtk_to_freelabels,
                                      input_names = ['hemi',
                                                     'surface_file',
                                                     'label_numbers',
                                                     'label_names',
                                                     'RGBs',
                                                     'scalar_name'],
                                      output_names = ['label_files',
                                                      'colortable']))
    annotflow.add_nodes([writelabels])
    mbflow.connect([(info, annotflow, [('hemi', 'Write_label_files.hemi')])])
    writelabels.inputs.label_numbers = ctx_label_numbers
    writelabels.inputs.label_names = ctx_label_names
    writelabels.inputs.RGBs = RGBs
    if init_labels == 'DKatlas':
        writelabels.inputs.scalar_name = 'Labels'
        mbflow.connect([(atlasflow, annotflow,
                         [('DK_annot_to_VTK.vtk_file',
                           'Write_label_files.surface_file')])])
    if init_labels == 'DKTatlas':
        writelabels.inputs.scalar_name = 'Labels'
        mbflow.connect([(atlasflow, annotflow,
                         [('DKT_annot_to_VTK.vtk_file',
                           'Write_label_files.surface_file')])])
    elif init_labels == 'max':
        writelabels.inputs.scalar_name = 'Max_(majority_labels)'
        mbflow.connect([(atlasflow, annotflow,
                         [('Label_vote.maxlabel_file',
                           'Write_label_files.surface_file')])])
    elif init_labels == 'manual':
        writelabels.inputs.scalar_name = 'Labels'
        mbflow.connect([(atlas, annotflow,
                         [('atlas_file',
                           'Write_label_files.surface_file')])])
    #---------------------------------------------------------------------------
    # Write .annot file from .label files
    # NOTE:  incorrect labels to be corrected below!
    #---------------------------------------------------------------------------
    writeannot = Node(name='Write_annot_file',
                      interface = Fn(function = labels_to_annot,
                                     input_names = ['hemi',
                                                    'subjects_path',
                                                    'subject',
                                                    'label_files',
                                                    'colortable',
                                                    'annot_name'],
                                     output_names = ['annot_name',
                                                     'annot_file']))
    writeannot.inputs.annot_name = 'labels.' + init_labels
    writeannot.inputs.subjects_path = subjects_path
    annotflow.add_nodes([writeannot])
    mbflow.connect([(info, annotflow,
                     [('hemi', 'Write_annot_file.hemi')])])
    mbflow.connect([(info, annotflow,
                     [('subject', 'Write_annot_file.subject')])])
    annotflow.connect([(writelabels, writeannot,
                      [('label_files','label_files')])])
    annotflow.connect([(writelabels, writeannot,
                      [('colortable','colortable')])])

################################################################################
#
#   Label volumes workflow:
#   * Fill volume
#   * Label evaluation
#
################################################################################
mbflow2 = Workflow(name='Label_volumes')
mbflow2.base_dir = temp_path

#-------------------------------------------------------------------------------
# Iterate inputs over subjects
#-------------------------------------------------------------------------------
info2 = info.clone('Inputs2')
info2.iterables = ([('subject', subjects)])
sink2 = sink.clone('Results2')

#-------------------------------------------------------------------------------
# Fill volume mask with surface vertex labels from .annot file
#-------------------------------------------------------------------------------
if fill_volume:

    fillvolume = Node(name='Fill_volume',
                      interface = Fn(function = labels_to_volume,
                                     input_names = ['subject',
                                                    'annot_name'],
                                     output_names = ['output_file']))
    mbflow2.add_nodes([fillvolume])
    fillvolume.inputs.annot_name = 'labels.' + init_labels
    mbflow2.connect([(info2, fillvolume, [('subject', 'subject')])])
    #---------------------------------------------------------------------------
    # Relabel file, replacing colortable labels with real labels
    #---------------------------------------------------------------------------
    relabel = Node(name='Correct_labels',
                   interface = Fn(function = relabel_volume,
                                  input_names = ['input_file',
                                                 'old_labels',
                                                 'new_labels'],
                                  output_names = ['output_file']))
    mbflow2.add_nodes([relabel])
    mbflow2.connect([(fillvolume, relabel,
                          [('output_file', 'input_file')])])
    relabel_file = os.path.join(info_path,
                                'label_volume_errors.' + protocol + '.txt')
    old_labels, new_labels = read_columns(relabel_file, 2)
    relabel.inputs.old_labels = old_labels
    relabel.inputs.new_labels = new_labels
    mbflow2.connect([(relabel, sink2,
                      [('output_file', 'labels_volume')])])

################################################################################
#
#   Volume label evaluation workflow
#
################################################################################
if evaluate_volume_labels:

    #---------------------------------------------------------------------------
    # Evaluation inputs: location and structure of atlas volumes
    #---------------------------------------------------------------------------
    atlas_vol = Node(name = 'Atlas',
                     interface = DataGrabber(infields=['subject'],
                     outfields=['atlas_vol_file']))
    atlas_vol.inputs.base_directory = atlases_path
    atlas_vol.inputs.template = '%s/mri/labels.' + protocol + '.manual.nii.gz'
    atlas_vol.inputs.template_args['atlas_vol_file'] = [['subject']]
    mbflow2.connect([(info2, atlas_vol, [('subject','subject')])])
    #---------------------------------------------------------------------------
    # Evaluate volume labels
    #---------------------------------------------------------------------------
    eval_vol_labels = Node(name='Evaluate_volume_labels',
                           interface = Fn(function = measure_volume_overlap,
                                          input_names = ['labels',
                                                         'atlas_file',
                                                         'input_file'],
                                          output_names = ['overlaps',
                                                          'out_file']))
    labels_file = os.path.join(info_path, 'labels.volume.' + protocol + '.txt')
    labels = read_columns(labels_file, 1)[0]
    eval_vol_labels.inputs.labels = labels
    mbflow2.add_nodes([eval_vol_labels])
    mbflow2.connect([(atlas_vol, eval_vol_labels,
                      [('atlas_vol_file','atlas_file')])])
    mbflow2.connect([(relabel, eval_vol_labels,
                      [('output_file', 'input_file')])])
    mbflow2.connect([(eval_vol_labels, sink2,
                      [('out_file', 'evaluate_labels_volume')])])

################################################################################
#
#    Run workflow
#
################################################################################
if __name__== '__main__':

    #from nipype import config, logging
    #config.set('logging', 'interface_level', 'DEBUG')
    #config.set('logging', 'workflow_level', 'DEBUG')
    #logging.update_logging(config)

    run_flow1 = True
    run_flow2 = True
    generate_graphs = True
    if generate_graphs:
        if run_flow1:
            mbflow.write_graph(graph2use='flat')
            mbflow.write_graph(graph2use='hierarchical')
        if run_flow2:
            mbflow2.write_graph(graph2use='flat')
            mbflow2.write_graph(graph2use='hierarchical')
    if run_flow1:
        mbflow.run()
    if run_flow2:
        mbflow2.run()
