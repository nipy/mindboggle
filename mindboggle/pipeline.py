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
fill_volume = 0#True  # Fill (gray matter) volumes with surface labels
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
hemis = ['lh'] #,'rh']  # Prepend ('lh.'/'rh.') indicating left/right surfaces
#-------------------------------------------------------------------------------
# Evaluation options
#-------------------------------------------------------------------------------
evaluate_surface_labels = 0 #False  # Surface overlap: auto vs. manual labels
evaluate_volume_labels = 0 #False  # Volume overlap: auto vs. manual labels
run_atlasFlow = True
run_measureFlow = True
run_featureFlow = True
run_shapeFlow = 0#True

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
from mindboggle.utils.io_vtk import rewrite_scalar_lists, write_mean_shapes_table, \
     load_scalars, freesurface_to_vtk, freecurvature_to_vtk, freeannot_to_vtk, \
     vtk_to_freelabels
from mindboggle.utils.io_file import read_columns
from mindboggle.utils.io_free import labels_to_annot, labels_to_volume
from mindboggle.utils.mesh_operations import find_neighbors
from mindboggle.label.multiatlas_labeling import register_template,\
     transform_atlas_labels, majority_vote_label
from mindboggle.info.sulcus_boundaries import sulcus_boundaries
from mindboggle.label.relabel import relabel_volume
from mindboggle.label.label_functions import label_with_classifier
from mindboggle.measure.measure_functions import compute_area, compute_depth, \
     compute_curvature
from mindboggle.extract.extract_folds import extract_folds, extract_sulci
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
mbFlow = Workflow(name='Mindboggle')
mbFlow.base_dir = temp_path
if not os.path.isdir(temp_path):  os.makedirs(temp_path)

#===============================================================================
#   Inputs and outputs
#===============================================================================
#-------------------------------------------------------------------------------
# Iterate inputs over subjects, hemispheres
# (surfaces are assumed to take the form: lh.pial or lh.pial.vtk)
#-------------------------------------------------------------------------------
Info = Node(name = 'Inputs',
            interface = IdentityInterface(fields=['subject', 'hemi']))
Info.iterables = ([('subject', subjects), ('hemi', hemis)])
#-------------------------------------------------------------------------------
# Location and structure of the surface inputs
#-------------------------------------------------------------------------------
Surf = Node(name = 'Surfaces',
            interface = DataGrabber(infields=['subject', 'hemi'],
                                    outfields=['surface_files', 'sphere_files']))
Surf.inputs.base_directory = subjects_path
Surf.inputs.template = '%s/surf/%s.%s'
Surf.inputs.template_args['surface_files'] = [['subject', 'hemi', 'pial']]
Surf.inputs.template_args['sphere_files'] = [['subject', 'hemi', 'sphere']]
if include_thickness:
    Surf.inputs.template_args['thickness_files'] = [['subject', 'hemi', 'thickness']]
if include_convexity:
    Surf.inputs.template_args['convexity_files'] = [['subject', 'hemi', 'sulc']]
mbFlow.connect([(Info, Surf, [('subject','subject'), ('hemi','hemi')])])
#-------------------------------------------------------------------------------
# Outputs
#-------------------------------------------------------------------------------
Sink = Node(DataSink(), name = 'Results')
Sink.inputs.base_directory = output_path
Sink.inputs.container = 'results'
if not os.path.isdir(output_path):  os.makedirs(output_path)
#-------------------------------------------------------------------------------
# Convert surfaces to VTK
#-------------------------------------------------------------------------------
if not input_vtk:
    ConvertSurf = Node(name = 'Surf_to_VTK',
                       interface = Fn(function = freesurface_to_vtk,
                                      input_names = ['surface_file'],
                                      output_names = ['vtk_file']))
    mbFlow.connect([(Surf, ConvertSurf, [('surface_files','surface_file')])])
#-------------------------------------------------------------------------------
# Evaluation inputs: location and structure of atlas surfaces
#-------------------------------------------------------------------------------
if evaluate_surface_labels or init_labels == 'manual':
    Atlas = Node(name = 'Atlases',
                 interface = DataGrabber(infields=['subject','hemi'],
                                         outfields=['atlas_file']))
    Atlas.inputs.base_directory = atlases_path

    Atlas.inputs.template = '%s/label/%s.labels.' +\
                            protocol + '.' + label_method + '.vtk'
    Atlas.inputs.template_args['atlas_file'] = [['subject','hemi']]

    mbFlow.connect([(Info, Atlas, [('subject','subject'),('hemi','hemi')])])
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
if run_atlasFlow:

    atlasFlow = Workflow(name='Label_initialization')

    #===========================================================================
    #   Initialize labels with FreeSurfer's standard DK classifier atlas
    #===========================================================================
    if init_labels == 'DKatlas':
        FreeLabels = Node(name = 'DK_annot_to_VTK',
                        interface = Fn(function = freeannot_to_vtk,
                                       input_names = ['surface_file',
                                                      'hemi',
                                                      'subject',
                                                      'subjects_path',
                                                      'annot_name'],
                                       output_names = ['labels',
                                                       'vtk_file']))
        atlasFlow.add_nodes([FreeLabels])
        if input_vtk:
            mbFlow.connect([(Surf, atlasFlow,
                             [('surface_files',
                               'DK_annot_to_VTK.surface_file')])])
        else:
            mbFlow.connect([(ConvertSurf, atlasFlow,
                             [('vtk_file',
                               'DK_annot_to_VTK.surface_file')])])
        mbFlow.connect([(Info, atlasFlow,
                         [('hemi', 'DK_annot_to_VTK.hemi'),
                          ('subject', 'DK_annot_to_VTK.subject')])])
        FreeLabels.inputs.subjects_path = subjects_path
        FreeLabels.inputs.annot_name = 'aparc.annot'
    #===========================================================================
    #   Initialize labels with the DKT classifier atlas
    #===========================================================================
    elif init_labels == 'DKTatlas':
        """
        Label a brain with the DKT atlas using FreeSurfer's mris_ca_label.
        """
        Classifier = Node(name = 'Label_with_DKTatlas',
                          interface = Fn(function = label_with_classifier,
                                         input_names = ['hemi',
                                                        'subject',
                                                        'subjects_path',
                                                        'sphere_file',
                                                        'classifier_path',
                                                        'classifier_atlas'],
                                         output_names = ['annot_name',
                                                         'annot_file']))
        atlasFlow.add_nodes([Classifier])
        mbFlow.connect([(Info, atlasFlow,
                         [('hemi', 'Label_with_DKTatlas.hemi'),
                          ('subject', 'Label_with_DKTatlas.subject')])])
        Classifier.inputs.subjects_path = subjects_path
        mbFlow.connect([(Surf, atlasFlow,
                         [('sphere_files',
                           'Label_with_DKTatlas.sphere_file')])])
        Classifier.inputs.classifier_path = classifier_path
        Classifier.inputs.classifier_atlas = classifier_atlas

        # Convert .annot file to .vtk format
        Classifier2vtk = Node(name = 'DKT_annot_to_VTK',
                              interface = Fn(function = freeannot_to_vtk,
                                             input_names = ['surface_file',
                                                            'hemi',
                                                            'subject',
                                                            'subjects_path',
                                                            'annot_name'],
                                             output_names = ['labels',
                                                             'vtk_file']))
        atlasFlow.add_nodes([Classifier2vtk])
        if input_vtk:
            mbFlow.connect([(Surf, atlasFlow,
                             [('surface_files',
                               'DKT_annot_to_VTK.surface_file')])])
        else:
            mbFlow.connect([(ConvertSurf, atlasFlow,
                             [('vtk_file',
                               'DKT_annot_to_VTK.surface_file')])])
        mbFlow.connect([(Info, atlasFlow,
                         [('hemi', 'DKT_annot_to_VTK.hemi'),
                          ('subject', 'DKT_annot_to_VTK.subject')])])
        Classifier2vtk.inputs.subjects_path = subjects_path
        atlasFlow.connect([(Classifier, Classifier2vtk,
                            [('annot_name', 'annot_name')])])
    #===========================================================================
    #   Initialize labels using multi-atlas registration
    #===========================================================================
    elif init_labels == 'max':
        #-----------------------------------------------------------------------
        # Register surfaces to average template
        #-----------------------------------------------------------------------
        free_template = 'OASIS-TRT-20'  # FreeSurfer template

        Register = Node(name = 'Register_template',
                        interface = Fn(function = register_template,
                                       input_names = ['hemi',
                                                      'sphere_file',
                                                      'transform',
                                                      'templates_path',
                                                      'template'],
                                       output_names = ['transform']))
        atlasFlow.add_nodes([Register])
        mbFlow.connect([(Info, atlasFlow, [('hemi', 'Register_template.hemi')]),
                        (Surf, atlasFlow, [('sphere_files',
                                            'Register_template.sphere_file')])])
        Register.inputs.transform = 'sphere_to_' + template + '_template.reg'
        Register.inputs.templates_path = os.path.join(templates_path, 'freesurfer')
        Register.inputs.template = free_template + '.tif'
        #-----------------------------------------------------------------------
        # Register atlases to subject via template
        #-----------------------------------------------------------------------
        # Load atlas list
        atlas_list_file = os.path.join(info_path, 'atlases.txt')
        atlas_list = read_columns(atlas_list_file, 1)[0]

        Transform = MapNode(name = 'Transform_labels',
                            iterfield = ['atlas'],
                            interface = Fn(function = transform_atlas_labels,
                                           input_names = ['hemi',
                                                          'subject',
                                                          'transform',
                                                          'subjects_path',
                                                          'atlas',
                                                          'atlas_string'],
                                           output_names = ['output_file']))
        atlasFlow.add_nodes([Transform])
        mbFlow.connect([(info, atlasFlow,
                         [('hemi', 'Transform_labels.hemi'),
                          ('subject', 'Transform_labels.subject')])])
        atlasFlow.connect([(Register, Transform, [('transform', 'transform')])])
        #Transform.inputs.transform = 'sphere_to_' + template + '_template.reg'
        Transform.inputs.subjects_path = subjects_path
        Transform.inputs.atlas = atlas_list
        Transform.inputs.atlas_string = 'labels.' + protocol + '.' + label_method
        #-----------------------------------------------------------------------
        # Majority vote label
        #-----------------------------------------------------------------------
        Vote = Node(name='Label_vote',
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
        atlasFlow.add_nodes([Vote])
        if input_vtk:
            mbFlow.connect([(Surf, atlasFlow,
                             [('surface_files', 'Label_vote.surface_file')])])
        else:
            mbFlow.connect([(ConvertSurf, atlasFlow,
                             [('vtk_file', 'Label_vote.surface_file')])])
        atlasFlow.connect([(Transform, Vote, [('output_file', 'annot_files')])])
        mbFlow.connect([(atlasFlow, Sink,
                         [('Label_vote.maxlabel_file', 'labels.@max'),
                          ('Label_vote.labelcounts_file', 'labels.@counts'),
                          ('Label_vote.labelvotes_file', 'labels.@votes')])])
    #===========================================================================
    #   Skip label initialization and process manual (atlas) labels
    #===========================================================================
    elif init_labels == 'manual':
        AtlasLabels = Node(name = 'Atlas_labels',
                           interface = Fn(function = load_scalars,
                                          input_names = ['filename',
                                                         'return_arrays'],
                                          output_names = ['points',
                                                          'faces',
                                                          'scalars',
                                                          'n_vertices']))
        atlasFlow.add_nodes([AtlasLabels])
        mbFlow.connect([(Atlas, atlasFlow,
                         [('atlas_file', 'Atlas_labels.filename')])])
        AtlasLabels.inputs.return_arrays = True  # 0: return lists instead of arrays


################################################################################
#
#   Surface measurement workflow
#
################################################################################
if run_measureFlow:

    measureFlow = Workflow(name='Surface_measurement')

    #===========================================================================
    #   Surface measurements
    #===========================================================================
    #---------------------------------------------------------------------------
    # Measure surface area
    #---------------------------------------------------------------------------
    AreaNode = Node(name='Area',
                interface = Fn(function = compute_area,
                               input_names = ['command',
                                              'surface_file'],
                               output_names = ['area_file']))
    area_command = os.path.join(ccode_path, 'area', 'PointAreaMain')
    AreaNode.inputs.command = area_command
    #---------------------------------------------------------------------------
    # Measure surface depth
    #---------------------------------------------------------------------------
    DepthNode = Node(name='Depth',
                 interface = Fn(function = compute_depth,
                                input_names = ['command',
                                               'surface_file'],
                                output_names = ['depth_file']))
    depth_command = os.path.join(ccode_path, 'travel_depth', 'TravelDepthMain')
    DepthNode.inputs.command = depth_command
    #---------------------------------------------------------------------------
    # Measure surface curvature
    #---------------------------------------------------------------------------
    CurvNode = Node(name='Curvature',
                     interface = Fn(function = compute_curvature,
                                    input_names = ['command',
                                                   'surface_file'],
                                    output_names = ['mean_curvature_file',
                                                    'gauss_curvature_file',
                                                    'max_curvature_file',
                                                    'min_curvature_file',
                                                    'min_curvature_vector_file']))
    curvature_command = os.path.join(ccode_path, 'curvature', 'CurvatureMain')
    CurvNode.inputs.command = curvature_command
    #---------------------------------------------------------------------------
    # Convert FreeSurfer surface measures to VTK
    #---------------------------------------------------------------------------
    if include_thickness:
        ThickNode = Node(name = 'Thickness_to_VTK',
                           interface = Fn(function = freecurvature_to_vtk,
                                          input_names = ['file_string',
                                                         'surface_file',
                                                         'hemi',
                                                         'subject',
                                                         'subjects_path'],
                                          output_names = ['vtk_file']))
        measureFlow.add_nodes([ThickNode])
        ThickNode.inputs.file_string = 'thickness'
        if not input_vtk:
            mbFlow.connect([(ConvertSurf, measureFlow,
                             [('vtk_file','Thickness_to_VTK.surface_file')])])
        else:
            mbFlow.connect([(Surf, measureFlow,
                             [('surface_files','Thickness_to_VTK.surface_file')])])
        mbFlow.connect([(Info, measureFlow,
                         [('hemi','Thickness_to_VTK.hemi'),
                          ('subject','Thickness_to_VTK.subject')])])
        ThickNode.inputs.subjects_path = subjects_path
        mbFlow.connect([(measureFlow, Sink,
                         [('Thickness_to_VTK.vtk_file', 'measures.@thickness')])])
    if include_convexity:
        ConvexNode = Node(name = 'Convexity_to_VTK',
                           interface = Fn(function = freecurvature_to_vtk,
                                          input_names = ['file_string',
                                                         'surface_file',
                                                         'hemi',
                                                         'subject',
                                                         'subjects_path'],
                                          output_names = ['vtk_file']))
        measureFlow.add_nodes([ConvexNode])
        ConvexNode.inputs.file_string = 'sulc'
        if not input_vtk:
            mbFlow.connect([(ConvertSurf, measureFlow,
                             [('vtk_file','Convexity_to_VTK.surface_file')])])
        else:
            mbFlow.connect([(Surf, measureFlow,
                             [('surface_files','Convexity_to_VTK.surface_file')])])
        mbFlow.connect([(Info, measureFlow, [('hemi','Convexity_to_VTK.hemi'),
                        ('subject','Convexity_to_VTK.subject')])])
        ConvexNode.inputs.subjects_path = subjects_path
        mbFlow.connect([(measureFlow, Sink,
                         [('Convexity_to_VTK.vtk_file', 'measures.@convexity')])])
    #---------------------------------------------------------------------------
    # Add and connect nodes, save output files
    #---------------------------------------------------------------------------
    measureFlow.add_nodes([AreaNode, DepthNode, CurvNode])
    if input_vtk:
        mbFlow.connect([(Surf, measureFlow,
                         [('surface_files','Area.surface_file')])])
        mbFlow.connect([(Surf, measureFlow,
                         [('surface_files','Depth.surface_file')])])
        mbFlow.connect([(Surf, measureFlow,
                         [('surface_files','Curvature.surface_file')])])
    else:
        mbFlow.connect([(ConvertSurf, measureFlow,
                         [('vtk_file', 'Area.surface_file')])])
        mbFlow.connect([(ConvertSurf, measureFlow,
                         [('vtk_file', 'Depth.surface_file')])])
        mbFlow.connect([(ConvertSurf, measureFlow,
                         [('vtk_file', 'Curvature.surface_file')])])
    mbFlow.connect([(measureFlow, Sink,
                     [('Area.area_file', 'measures.@area')])])
    mbFlow.connect([(measureFlow, Sink,
                     [('Depth.depth_file', 'measures.@depth')])])
    mbFlow.connect([(measureFlow, Sink,
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
if run_featureFlow:

    featureFlow = Workflow(name='Feature_extraction')

    #===========================================================================
    # Load surface and find all vertex neighbors
    #===========================================================================
    LoadSurf = Node(name = 'Load_surface',
                        interface = Fn(function = load_scalars,
                                       input_names = ['filename',
                                                      'return_arrays'],
                                       output_names = ['points',
                                                       'faces',
                                                       'scalars',
                                                       'n_vertices']))
    featureFlow.add_nodes([LoadSurf])
    if input_vtk:
        mbFlow.connect([(Surf, featureFlow,
                         [('surface_files','Load_surface.filename')])])
    else:
        mbFlow.connect([(ConvertSurf, featureFlow,
                         [('vtk_file', 'Load_surface.filename')])])
    LoadSurf.inputs.return_arrays = True  # 0: return lists instead of arrays

    NbrNode = Node(name='Neighbors',
                     interface = Fn(function = find_neighbors,
                                    input_names = ['faces', 'n_vertices'],
                                    output_names = ['neighbor_lists']))
    featureFlow.add_nodes([NbrNode])
    featureFlow.connect([(LoadSurf, NbrNode,
                          [('faces','faces'), ('n_vertices','n_vertices')])])

    #===========================================================================
    # Extract folds
    #===========================================================================
    fraction_folds = 0.5
    min_fold_size = 50
    FoldsNode = Node(name='Folds',
                 interface = Fn(function = extract_folds,
                                input_names = ['depth_file',
                                               'area_file',
                                               'neighbor_lists',
                                               'fraction_folds',
                                               'min_fold_size',
                                               'do_fill_holes'],
                                output_names = ['folds',
                                                'n_folds']))
    featureFlow.add_nodes([FoldsNode])
    mbFlow.connect([(measureFlow, featureFlow,
                     [('Depth.depth_file','Folds.depth_file')])])
    mbFlow.connect([(measureFlow, featureFlow,
                     [('Area.area_file','Folds.area_file')])])
    featureFlow.connect([(NbrNode, FoldsNode,
                          [('neighbor_lists','neighbor_lists')])])
    FoldsNode.inputs.fraction_folds = fraction_folds
    FoldsNode.inputs.min_fold_size = min_fold_size
    FoldsNode.inputs.do_fill_holes = True

    #===========================================================================
    # Extract sulci from folds
    #===========================================================================
    min_sulcus_size = 50
    LabelPairs = Node(name='Label_pairs',
                       interface = Fn(function = sulcus_boundaries,
                                      input_names = [],
                                      output_names = ['label_pair_lists']))
    featureFlow.add_nodes([LabelPairs])

    SulciNode = Node(name='Sulci',
                 interface = Fn(function = extract_sulci,
                                input_names = ['surface_vtk',
                                               'folds',
                                               'labels',
                                               'neighbor_lists',
                                               'label_pair_lists',
                                               'min_boundary',
                                               'sulcus_names=[]'],
                                output_names = ['sulci',
                                                'n_sulci']))
    featureFlow.add_nodes([SulciNode])
    mbFlow.connect([(measureFlow, featureFlow,
                     [('Depth.depth_file','SulciNode.surface_vtk')])])
    featureFlow.connect([(FoldsNode, SulciNode, [('folds','folds')])])
    #---------------------------------------------------------------------------
    # Use initial labels assigned by FreeSurfer classifier atlas
    if init_labels == 'DKatlas':
        mbFlow.connect([(atlasFlow, featureFlow,
                         [('DK_annot_to_VTK.labels','SulciNode.labels')])])
    # Use initial labels assigned by Mindboggle classifier atlas
    elif init_labels == 'DKTatlas':
        mbFlow.connect([(atlasFlow, featureFlow,
                         [('DKT_annot_to_VTK.labels','SulciNode.labels')])])
    # Use initial labels assigned by multi-atlas registration
    elif init_labels == 'max':
        mbFlow.connect([(atlasFlow, featureFlow,
                         [('Label_vote.labels_max','SulciNode.labels')])])
    # Use manual (atlas) labels
    elif init_labels == 'manual':
        mbFlow.connect([(atlasFlow, featureFlow,
                         [('Atlas_labels.scalars','SulciNode.labels')])])
    #---------------------------------------------------------------------------
    featureFlow.connect([(NbrNode, SulciNode,
                          [('neighbor_lists','neighbor_lists')])])
    featureFlow.connect([(LabelPairs, SulciNode,
                          [('label_pair_lists','label_pair_lists')])])
    SulciNode.inputs.min_boundary = min_boundary
    SulciNode.inputs.sulcus_names = []

    #===========================================================================
    # Extract fundi (curves at the bottoms of sulci)
    #===========================================================================
    thr = 0.5
    min_distance = 5.0
    fundi_from_sulci = True
    FundiNode = Node(name='Fundi',
                     interface = Fn(function = extract_fundi,
                                    input_names = ['folds',
                                                   'neighbor_lists',
                                                   'depth_file',
                                                   'mean_curvature_file',
                                                   'min_curvature_vector_file',
                                                   'min_distance',
                                                   'thr',
                                                   'use_only_endpoints'],
                                    output_names = ['fundi',
                                                    'n_fundi',
                                                    'likelihoods']))
    if fundi_from_sulci:
        featureFlow.connect([(SulciNode, FundiNode, [('sulci','folds')])])
    else:
        featureFlow.connect([(FoldsNode, FundiNode, [('folds','folds')])])
    featureFlow.connect([(NbrNode, FundiNode,
                          [('neighbor_lists','neighbor_lists')])])
    mbFlow.connect([(measureFlow, featureFlow,
                     [('Depth.depth_file','Fundi.depth_file'),
                      ('Curvature.mean_curvature_file',
                       'Fundi.mean_curvature_file'),
                      ('Curvature.min_curvature_vector_file',
                       'Fundi.min_curvature_vector_file')])])
    FundiNode.inputs.min_distance = min_distance
    FundiNode.inputs.thr = thr
    FundiNode.inputs.use_only_endpoints = True

    #---------------------------------------------------------------------------
    # Write folds, sulci, fundi, and likelihoods to VTK files
    #---------------------------------------------------------------------------
    SulciVTK = Node(name='Sulci_to_VTK',
                    interface = Fn(function = rewrite_scalar_lists,
                                   input_names = ['input_vtk',
                                                  'output_vtk',
                                                  'new_scalar_lists',
                                                  'new_scalar_names',
                                                  'filter_scalars'],
                                   output_names = ['output_vtk']))
    # Save sulci
    featureFlow.add_nodes([SulciVTK])
    mbFlow.connect([(measureFlow, featureFlow,
                     [('Depth.depth_file','Sulci_to_VTK.input_vtk')])])
    SulciVTK.inputs.output_vtk = 'sulci.vtk'
    SulciVTK.inputs.new_scalar_names = ['sulci']
    featureFlow.connect([(SulciNode, SulciVTK, [('sulci','new_scalar_lists')])])
    featureFlow.connect([(SulciNode, SulciVTK, [('sulci','filter_scalars')])])
    mbFlow.connect([(featureFlow, Sink,
                     [('Sulci_to_VTK.output_vtk','features.@sulci')])])

    # Save fundi
    FundiVTK = SulciVTK.clone('Fundi_to_VTK')
    featureFlow.add_nodes([FundiVTK])
    mbFlow.connect([(measureFlow, featureFlow,
                     [('Depth.depth_file','Fundi_to_VTK.input_vtk')])])
    FundiVTK.inputs.output_vtk = 'fundi.vtk'
    FundiVTK.inputs.new_scalar_names = ['fundi']
    featureFlow.connect([(FundiNode, FundiVTK, [('fundi','new_scalar_lists')])])
    featureFlow.connect([(FundiNode, FundiVTK, [('fundi','filter_scalars')])])
    mbFlow.connect([(featureFlow, Sink,
                     [('Fundi_to_VTK.output_vtk','features.@fundi')])])

    # Save likelihoods values in sulci (or in folds if fundi_from_sulci==False)
    LikelihoodsVTK = SulciVTK.clone('Likelihoods_to_VTK')
    featureFlow.add_nodes([LikelihoodsVTK])
    mbFlow.connect([(measureFlow, featureFlow,
                     [('Depth.depth_file','Likelihoods_to_VTK.input_vtk')])])
    LikelihoodsVTK.inputs.output_vtk = 'likelihoods.vtk'
    LikelihoodsVTK.inputs.new_scalar_names = ['likelihoods']
    featureFlow.connect([(FundiNode, LikelihoodsVTK, [('likelihoods','new_scalar_lists')])])
    featureFlow.connect([(FundiNode, LikelihoodsVTK, [('likelihoods','filter_scalars')])])
    mbFlow.connect([(featureFlow, Sink,
                     [('Likelihoods_to_VTK.output_vtk','features.@likelihoods')])])

################################################################################
#
#   Shape analysis workflow
#
################################################################################
if run_shapeFlow:

    shapeFlow = Workflow(name='Shape_analysis')
    column_names = ['labels', 'area', 'depth', 'mean_curvature',
                    'gauss_curvature', 'max_curvature', 'min_curvature']
    if include_thickness:
        column_names.append('thickness')
    if include_convexity:
        column_names.append('convexity')
    shape_files = [x + '_file' for x in column_names[2::]]
    input_names = ['filename', 'column_names', 'labels', 'exclude_values', 'area_file']
    input_names.extend(shape_files)

    #===========================================================================
    # Labeled surface patch shapes
    #===========================================================================
    LabelTable = Node(name='Label_table',
                      interface = Fn(function = write_mean_shapes_table,
                                     input_names = input_names,
                                     output_names = ['means_file',
                                                     'norm_means_file']))
    shapeFlow.add_nodes([LabelTable])
    LabelTable.inputs.filename = 'label_shapes.txt'
    LabelTable.inputs.column_names = column_names
    #---------------------------------------------------------------------------
    # Use initial labels assigned by FreeSurfer classifier atlas
    if init_labels == 'DKatlas':
        mbFlow.connect([(atlasFlow, shapeFlow,
                         [('DK_annot_to_VTK.labels','Label_table.labels')])])
    # Use initial labels assigned by Mindboggle classifier atlas
    elif init_labels == 'DKTatlas':
        mbFlow.connect([(atlasFlow, shapeFlow,
                         [('DKT_annot_to_VTK.labels','Label_table.labels')])])
    # Use initial labels assigned by multi-atlas registration
    elif init_labels == 'max':
        mbFlow.connect([(atlasFlow, shapeFlow,
                         [('Label_vote.labels_max','Label_table.labels')])])
    # Use manual (atlas) labels
    elif init_labels == 'manual':
        mbFlow.connect([(atlasFlow, shapeFlow,
                         [('Atlas_labels.scalars','Label_table.labels')])])
    #---------------------------------------------------------------------------
    LabelTable.inputs.exclude_values = [-1,0]
    mbFlow.connect([(measureFlow, shapeFlow,
                     [('Area.area_file','Label_table.area_file')])])
    mbFlow.connect([(measureFlow, shapeFlow,
                     [('Depth.depth_file','Label_table.depth_file')])])
    mbFlow.connect([(measureFlow, shapeFlow,
                     [('Curvature.mean_curvature_file',
                       'Label_table.mean_curvature_file')])])
    mbFlow.connect([(measureFlow, shapeFlow,
                     [('Curvature.gauss_curvature_file',
                       'Label_table.gauss_curvature_file')])])
    mbFlow.connect([(measureFlow, shapeFlow,
                     [('Curvature.max_curvature_file',
                       'Label_table.max_curvature_file')])])
    mbFlow.connect([(measureFlow, shapeFlow,
                     [('Curvature.min_curvature_file',
                       'Label_table.min_curvature_file')])])
    if include_thickness:
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Thickness_to_VTK.vtk_file',
                           'Label_table.thickness_file')])])
    if include_convexity:
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Convexity_to_VTK.vtk_file',
                           'Label_table.convexity_file')])])
    # Save results
    mbFlow.connect([(shapeFlow, Sink,
                     [('Label_table.means_file', 'shapes.@labels'),
                      ('Label_table.norm_means_file', 'shapes.@labels_norm')])])

    #===========================================================================
    # Sulcus fold shapes
    #===========================================================================
    if run_featureFlow:
        SulcusTable = LabelTable.clone('Sulcus_table')
        shapeFlow.add_nodes([SulcusTable])
        SulcusTable.inputs.filename = 'sulcus_shapes.txt'
        SulcusTable.inputs.column_names = column_names
        mbFlow.connect([(featureFlow, shapeFlow,
                         [('Sulci.sulci','Sulcus_table.labels')])])
        SulcusTable.inputs.exclude_values = [-1]
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Area.area_file','Sulcus_table.area_file')])])
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Depth.depth_file','Sulcus_table.depth_file')])])
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Curvature.mean_curvature_file',
                           'Sulcus_table.mean_curvature_file')])])
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Curvature.gauss_curvature_file',
                           'Sulcus_table.gauss_curvature_file')])])
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Curvature.max_curvature_file',
                           'Sulcus_table.max_curvature_file')])])
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Curvature.min_curvature_file',
                           'Sulcus_table.min_curvature_file')])])
        if include_thickness:
            mbFlow.connect([(measureFlow, shapeFlow,
                             [('Thickness_to_VTK.vtk_file',
                               'Sulcus_table.thickness_file')])])
        if include_convexity:
            mbFlow.connect([(measureFlow, shapeFlow,
                             [('Convexity_to_VTK.vtk_file',
                               'Sulcus_table.convexity_file')])])
        # Save results
        mbFlow.connect([(shapeFlow, Sink,
                         [('Sulcus_table.means_file', 'shapes.@sulci'),
                          ('Sulcus_table.norm_means_file', 'shapes.@sulci_norm')])])

    #===========================================================================
    # Fundus shapes
    #===========================================================================
    if run_featureFlow:
        FundusTable = LabelTable.clone('Fundus_table')
        shapeFlow.add_nodes([FundusTable])
        FundusTable.inputs.filename = 'fundus_shapes.txt'
        FundusTable.inputs.column_names = column_names
        mbFlow.connect([(featureFlow, shapeFlow,
                         [('Fundi.fundi','Fundus_table.labels')])])
        FundusTable.inputs.exclude_values = [-1]
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Area.area_file','Fundus_table.area_file')])])
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Depth.depth_file','Fundus_table.depth_file')])])
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Curvature.mean_curvature_file',
                           'Fundus_table.mean_curvature_file')])])
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Curvature.gauss_curvature_file',
                           'Fundus_table.gauss_curvature_file')])])
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Curvature.max_curvature_file',
                           'Fundus_table.max_curvature_file')])])
        mbFlow.connect([(measureFlow, shapeFlow,
                         [('Curvature.min_curvature_file',
                           'Fundus_table.min_curvature_file')])])
        if include_thickness:
            mbFlow.connect([(measureFlow, shapeFlow,
                             [('Thickness_to_VTK.vtk_file',
                               'Fundus_table.thickness_file')])])
        if include_convexity:
            mbFlow.connect([(measureFlow, shapeFlow,
                             [('Convexity_to_VTK.vtk_file',
                               'Fundus_table.convexity_file')])])
        # Save results
        mbFlow.connect([(shapeFlow, Sink,
                         [('Fundus_table.means_file', 'shapes.@fundi'),
                          ('Fundus_table.norm_means_file', 'shapes.@fundi_norm')])])

################################################################################
#
#   Surface label evaluation
#
################################################################################
if evaluate_surface_labels:

    #===========================================================================
    # Evaluate surface labels
    #===========================================================================
    EvalSurfLabels = Node(name='Evaluate_surface_labels',
                            interface = Fn(function = measure_surface_overlap,
                                           input_names = ['command',
                                                          'labels_file1',
                                                          'labels_file2'],
                                           output_names = ['overlap_file']))
    mbFlow.add_nodes([EvalSurfLabels])
    surface_overlap_command = os.path.join(ccode_path,
        'surface_overlap', 'SurfaceOverlapMain')
    EvalSurfLabels.inputs.command = surface_overlap_command
    mbFlow.connect([(Atlas, EvalSurfLabels, [('atlas_file','labels_file1')])])
    if init_labels == 'DKatlas':
        mbFlow.connect([(atlasFlow, EvalSurfLabels,
                         [('DK_annot_to_VTK.vtk_file','labels_file2')])])
    elif init_labels == 'DKTatlas':
        mbFlow.connect([(atlasFlow, EvalSurfLabels,
                         [('DKT_annot_to_VTK.vtk_file','labels_file2')])])
    elif init_labels == 'max':
        mbFlow.connect([(atlasFlow, EvalSurfLabels,
                         [('Label_vote.maxlabel_file','labels_file2')])])
    elif init_labels == 'manual':
        mbFlow.connect([(Atlas, EvalSurfLabels,
                         [('atlas_file','labels_file2')])])
    mbFlow.connect([(EvalSurfLabels, Sink,
                     [('overlap_file', 'evaluate_labels')])])

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
    WriteLabels = Node(name='Write_label_files',
                       interface = Fn(function = vtk_to_freelabels,
                                      input_names = ['hemi',
                                                     'surface_file',
                                                     'label_numbers',
                                                     'label_names',
                                                     'RGBs',
                                                     'scalar_name'],
                                      output_names = ['label_files',
                                                      'colortable']))
    annotflow.add_nodes([WriteLabels])
    mbFlow.connect([(Info, annotflow, [('hemi', 'Write_label_files.hemi')])])
    WriteLabels.inputs.label_numbers = ctx_label_numbers
    WriteLabels.inputs.label_names = ctx_label_names
    WriteLabels.inputs.RGBs = RGBs
    if init_labels == 'DKatlas':
        WriteLabels.inputs.scalar_name = 'Labels'
        mbFlow.connect([(atlasFlow, annotflow,
                         [('DK_annot_to_VTK.vtk_file',
                           'Write_label_files.surface_file')])])
    if init_labels == 'DKTatlas':
        WriteLabels.inputs.scalar_name = 'Labels'
        mbFlow.connect([(atlasFlow, annotflow,
                         [('DKT_annot_to_VTK.vtk_file',
                           'Write_label_files.surface_file')])])
    elif init_labels == 'max':
        WriteLabels.inputs.scalar_name = 'Max_(majority_labels)'
        mbFlow.connect([(atlasFlow, annotflow,
                         [('Label_vote.maxlabel_file',
                           'Write_label_files.surface_file')])])
    elif init_labels == 'manual':
        WriteLabels.inputs.scalar_name = 'Labels'
        mbFlow.connect([(Atlas, annotflow,
                         [('atlas_file',
                           'Write_label_files.surface_file')])])
    #---------------------------------------------------------------------------
    # Write .annot file from .label files
    # NOTE:  incorrect labels to be corrected below!
    #---------------------------------------------------------------------------
    WriteAnnot = Node(name='Write_annot_file',
                      interface = Fn(function = labels_to_annot,
                                     input_names = ['hemi',
                                                    'subjects_path',
                                                    'subject',
                                                    'label_files',
                                                    'colortable',
                                                    'annot_name'],
                                     output_names = ['annot_name',
                                                     'annot_file']))
    WriteAnnot.inputs.annot_name = 'labels.' + init_labels
    WriteAnnot.inputs.subjects_path = subjects_path
    annotflow.add_nodes([WriteAnnot])
    mbFlow.connect([(Info, annotflow,
                     [('hemi', 'Write_annot_file.hemi')])])
    mbFlow.connect([(Info, annotflow,
                     [('subject', 'Write_annot_file.subject')])])
    annotflow.connect([(WriteLabels, WriteAnnot,
                      [('label_files','label_files')])])
    annotflow.connect([(WriteLabels, WriteAnnot,
                      [('colortable','colortable')])])

################################################################################
#
#   Label volumes workflow:
#   * Fill volume
#   * Label evaluation
#
################################################################################
mbFlow2 = Workflow(name='Label_volumes')
mbFlow2.base_dir = temp_path

#-------------------------------------------------------------------------------
# Iterate inputs over subjects
#-------------------------------------------------------------------------------
Info2 = Info.clone('Inputs2')
Info2.iterables = ([('subject', subjects)])
Sink2 = Sink.clone('Results2')

#-------------------------------------------------------------------------------
# Fill volume mask with surface vertex labels from .annot file
#-------------------------------------------------------------------------------
if fill_volume:

    FillVolume = Node(name='Fill_volume',
                      interface = Fn(function = labels_to_volume,
                                     input_names = ['subject',
                                                    'annot_name'],
                                     output_names = ['output_file']))
    mbFlow2.add_nodes([FillVolume])
    FillVolume.inputs.annot_name = 'labels.' + init_labels
    mbFlow2.connect([(Info2, FillVolume, [('subject', 'subject')])])
    #---------------------------------------------------------------------------
    # Relabel file, replacing colortable labels with real labels
    #---------------------------------------------------------------------------
    Relabel = Node(name='Correct_labels',
                   interface = Fn(function = relabel_volume,
                                  input_names = ['input_file',
                                                 'old_labels',
                                                 'new_labels'],
                                  output_names = ['output_file']))
    mbFlow2.add_nodes([Relabel])
    mbFlow2.connect([(FillVolume, Relabel,
                          [('output_file', 'input_file')])])
    relabel_file = os.path.join(info_path,
                                'label_volume_errors.' + protocol + '.txt')
    old_labels, new_labels = read_columns(relabel_file, 2)
    Relabel.inputs.old_labels = old_labels
    Relabel.inputs.new_labels = new_labels
    mbFlow2.connect([(Relabel, Sink2,
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
    AtlasVol = Node(name = 'Atlas_volume',
                     interface = DataGrabber(infields=['subject'],
                     outfields=['atlas_vol_file']))
    AtlasVol.inputs.base_directory = atlases_path
    AtlasVol.inputs.template = '%s/mri/labels.' + protocol + '.manual.nii.gz'
    AtlasVol.inputs.template_args['atlas_vol_file'] = [['subject']]
    mbFlow2.connect([(Info2, AtlasVol, [('subject','subject')])])
    #---------------------------------------------------------------------------
    # Evaluate volume labels
    #---------------------------------------------------------------------------
    EvalVolLabels = Node(name='Evaluate_volume_labels',
                           interface = Fn(function = measure_volume_overlap,
                                          input_names = ['labels',
                                                         'atlas_file',
                                                         'input_file'],
                                          output_names = ['overlaps',
                                                          'out_file']))
    labels_file = os.path.join(info_path, 'labels.volume.' + protocol + '.txt')
    labels = read_columns(labels_file, 1)[0]
    EvalVolLabels.inputs.labels = labels
    mbFlow2.add_nodes([EvalVolLabels])
    mbFlow2.connect([(AtlasVol, EvalVolLabels,
                      [('atlas_vol_file','atlas_file')])])
    mbFlow2.connect([(Relabel, EvalVolLabels,
                      [('output_file', 'input_file')])])
    mbFlow2.connect([(EvalVolLabels, Sink2,
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
            mbFlow.write_graph(graph2use='flat')
            mbFlow.write_graph(graph2use='hierarchical')
        if run_flow2:
            mbFlow2.write_graph(graph2use='flat')
            mbFlow2.write_graph(graph2use='hierarchical')
    if run_flow1:
        mbFlow.run()
    if run_flow2:
        mbFlow2.run()
