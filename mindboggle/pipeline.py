#!/usr/bin/env python
"""
This is Mindboggle's nipype software pipeline!

Examples
--------
$ python pipeline.py <output path> <1 or more subject names>
$ python pipeline.py output HLN-12-1 HLN-12-2

..Mindboggle surface workflows ::

    * Label surfaces:
        - call manual labels OR
        - call labels from FreeSurfer OR
        - multi-atlas labeling OR
        - DKT40 atlas labeling (default)
        -> evaluate surface labels

    * Extract features:
        - folds
        - fundi
        - sulci

    * Measure shapes:
        - travel depth
        - curvatures (mean, Gaussian, min, max, min directions)
        - area
        - thickness (from FreeSurfer)
        - convexity (from FreeSurfer)

    * Construct tables:
        - label shapes
        - sulcus shapes
        - fundus shapes
        - per-vertex shape measures

..Mindboggle volume workflows ::

    * Label volume workflow:
        - fill gray matter with labels
        - measure label volumes
        - construct tables
        -> evaluate label volumes


.. Note::
      Mindboggle assumes a file tree like FreeSurfer's,
      and for label initialization, assumes that subjects have been processed
      by FreeSurfer (autorecon -all), so subject names correspond to directory
      names in FreeSurfer's subjects directory.

For more information about Mindboggle (http://mindboggle.info)
and read the documentation: http://mindboggle.info/software/documentation.html

For information on Nipype (http://www.nipy.org/nipype/):
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3159964/


Authors:
    - Arno Klein, 2011-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#=============================================================================
# Command line arguments
#=============================================================================
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

#=============================================================================
# User settings
#=============================================================================
do_input_vtk = False  # Load VTK surfaces directly (not FreeSurfer surfaces)
do_fundi = False  # Extract fundi
do_sulci = True  # Extract sulci
do_thickness = True  # Include FreeSurfer's thickness measure
do_convexity = True  # Include FreeSurfer's convexity measure (sulc.pial)
do_measure_spectra = False  # Measure Laplace-Beltrami spectra for features
do_vertex_tables = True  # Create per-vertex shape tables
do_fill = True  # Fill (gray matter) volumes with surface labels (FreeSurfer)
do_measure_volume = True  # Measure volumes of labeled regions
do_evaluate_surface = False  # Surface overlap: auto vs. manual labels
do_evaluate_volume = False  # Volume overlap: auto vs. manual labels
#-----------------------------------------------------------------------------
# Mindboggle workflows
#-----------------------------------------------------------------------------
run_labelFlow = True
run_shapeFlow = True
run_featureFlow = True
run_tableFlow = True
run_volumeFlow = True
#-----------------------------------------------------------------------------
# Labeling protocol used by Mindboggle:
# 'DKT31': 'Desikan-Killiany-Tourville (DKT) protocol with 31 labeled regions
# 'DKT25': 'fundus-friendly' version of the DKT protocol following fundi
#-----------------------------------------------------------------------------
protocol = 'DKT25'
#-----------------------------------------------------------------------------
# Initialize labels with:
# 'DKatlas': the standard FreeSurfer classifier atlas trained on the DK protocol
# 'DKTatlas': a FreeSurfer-style classifier atlas trained on the DKT protocol
# 'max': maximum probability (majority vote) labels from multiple atlases
# 'manual': process manual labels (atlas)
#-----------------------------------------------------------------------------
init_labels = 'DKTatlas'
#-----------------------------------------------------------------------------
# Labeling source:
# 'manual': manual edits
# FUTURE:
# <'adjusted': manual edits after automated alignment to fundi>
#-----------------------------------------------------------------------------
label_method = 'manual'
hemis = ['lh','rh']  # Prepend ('lh.'/'rh.') indicating left/right surfaces
#-----------------------------------------------------------------------------
# Evaluation options
#-----------------------------------------------------------------------------

#=============================================================================
# Setup: import libraries, set file paths, and initialize main workflow
#=============================================================================
#-----------------------------------------------------------------------------
# Import system and nipype Python libraries
#-----------------------------------------------------------------------------
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import Function as Fn
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import DataGrabber, DataSink
#-----------------------------------------------------------------------------
# Import Mindboggle Python libraries
#-----------------------------------------------------------------------------
from mindboggle.utils.io_vtk import rewrite_scalars, read_vtk
from mindboggle.utils.io_file import read_columns, write_columns
from mindboggle.utils.io_free import labels_to_annot, labels_to_volume, \
    surface_to_vtk, curvature_to_vtk, annot_to_vtk, vtk_to_labels
from mindboggle.utils.mesh import find_neighbors_from_file
from mindboggle.labels.label_free import register_template,\
     transform_atlas_labels, label_with_classifier
from mindboggle.labels.labels import majority_vote_label
from mindboggle.labels.protocol.sulci_labelpairs_DKT import sulcus_boundaries
from mindboggle.labels.relabel import relabel_volume
from mindboggle.shapes.measure import area, depth, curvature,\
    volume_per_label, rescale_by_neighborhood
from mindboggle.shapes.tabulate import write_mean_shapes_tables, \
    write_vertex_shapes_table
from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_labels
from mindboggle.features.folds import extract_folds
from mindboggle.shapes.likelihood import compute_likelihood
from mindboggle.features.fundi import extract_fundi
from mindboggle.features.sulci import extract_sulci
from mindboggle.evaluate.evaluate_labels import measure_surface_overlap, \
     measure_volume_overlap

#from mindboggle import get_info
#-----------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------
subjects_path = os.environ['SUBJECTS_DIR']  # FreeSurfer subjects directory
atlases_path = subjects_path  # Mindboggle-101 atlases directory
data_path = os.environ['MINDBOGGLE_DATA']  # Mindboggle data directory
ccode_path = os.environ['MINDBOGGLE_TOOLS']  # Mindboggle C++ code directory
protocol_path = os.path.join(os.environ['MINDBOGGLE'], 'labels', 'protocol')
x_path = os.path.join(os.environ['MINDBOGGLE'], 'x')
temp_path = os.path.join(output_path, 'workspace')  # Where to save temp files
# Label with classifier atlas
templates_path = os.path.join(data_path, 'atlases')
classifier_path = os.path.join(data_path, 'atlases')
classifier_atlas = 'DKTatlas40.gcs'  # 'DKTatlas100.gcs'
#-----------------------------------------------------------------------------
# Initialize main workflow
#-----------------------------------------------------------------------------
mbFlow = Workflow(name='Mindboggle')
mbFlow.base_dir = temp_path
if not os.path.isdir(temp_path):  os.makedirs(temp_path)

#=============================================================================
#  Inputs and outputs
#=============================================================================
#-----------------------------------------------------------------------------
# Iterate inputs over subjects, hemispheres
# (surfaces are assumed to take the form: lh.pial or lh.pial.vtk)
#-----------------------------------------------------------------------------
Info = Node(name = 'Inputs',
            interface = IdentityInterface(fields=['subject', 'hemi']))
Info.iterables = ([('subject', subjects), ('hemi', hemis)])
#-----------------------------------------------------------------------------
# Location and structure of the surface inputs
#-----------------------------------------------------------------------------
Surf = Node(name = 'Surfaces',
            interface = DataGrabber(infields=['subject', 'hemi'],
                                    outfields=['surface_files', 'sphere_files']))
Surf.inputs.base_directory = subjects_path
Surf.inputs.template = '%s/surf/%s.%s'
Surf.inputs.template_args['surface_files'] = [['subject', 'hemi', 'pial']]
Surf.inputs.template_args['sphere_files'] = [['subject', 'hemi', 'sphere']]
if do_thickness:
    Surf.inputs.template_args['thickness_files'] = [['subject', 'hemi', 'thickness']]
if do_convexity:
    Surf.inputs.template_args['convexity_files'] = [['subject', 'hemi', 'sulc']]
mbFlow.connect([(Info, Surf, [('subject','subject'), ('hemi','hemi')])])
#-----------------------------------------------------------------------------
# Location and structure of the FreeSurfer label inputs
#-----------------------------------------------------------------------------
Annot = Node(name = 'Annots',
             interface = DataGrabber(infields=['subject', 'hemi'],
                                     outfields=['annot_files']))
Annot.inputs.base_directory = subjects_path
Annot.inputs.template = '%s/label/%s.aparc.annot'
Annot.inputs.template_args['annot_files'] = [['subject','hemi']]
#-----------------------------------------------------------------------------
# Location and structure of the volume inputs
#-----------------------------------------------------------------------------
if do_fill:
    Vol = Node(name = 'Volumes',
        interface = DataGrabber(infields=['subject'],
                                outfields=['original_volume']))
    Vol.inputs.base_directory = subjects_path
    Vol.inputs.template = '%s/mri/orig/001.mgz'
    Vol.inputs.template_args['original_volume'] = [['subject']]
#-----------------------------------------------------------------------------
# Outputs
#-----------------------------------------------------------------------------
Sink = Node(DataSink(), name = 'Results')
Sink.inputs.base_directory = output_path
Sink.inputs.container = 'results'
if not os.path.isdir(output_path):  os.makedirs(output_path)
#-----------------------------------------------------------------------------
# Convert surfaces to VTK
#-----------------------------------------------------------------------------
if not do_input_vtk:
    ConvertSurf = Node(name = 'Surface_to_VTK',
                       interface = Fn(function = surface_to_vtk,
                                      input_names = ['surface_file'],
                                      output_names = ['vtk_file']))
    mbFlow.connect([(Surf, ConvertSurf, [('surface_files','surface_file')])])
#-----------------------------------------------------------------------------
# Evaluation inputs: location and structure of atlas surfaces
#-----------------------------------------------------------------------------
if do_evaluate_surface or init_labels == 'manual':
    Atlas = Node(name = 'Atlases',
                 interface = DataGrabber(infields=['subject','hemi'],
                                         outfields=['atlas_file']))
    Atlas.inputs.base_directory = atlases_path

    Atlas.inputs.template = '%s/label/%s.labels.' +\
                            protocol + '.' + label_method + '.vtk'
    Atlas.inputs.template_args['atlas_file'] = [['subject','hemi']]

    mbFlow.connect([(Info, Atlas, [('subject','subject'),('hemi','hemi')])])
#-----------------------------------------------------------------------------
# Load data
#-----------------------------------------------------------------------------
ctx_labels_file = os.path.join(protocol_path, 'labels.surface.' + protocol + '.txt')
ctx_label_numbers, ctx_label_names, RGBs = read_columns(ctx_labels_file,
                                                n_columns=3, trail=True)

##############################################################################
#
#   Label workflow
#
##############################################################################
if run_labelFlow:

    labelFlow = Workflow(name='Labels')

    #=========================================================================
    # Initialize labels with FreeSurfer's standard DK classifier atlas
    #=========================================================================
    if init_labels == 'DKatlas':
        FreeLabels = Node(name = 'DK_annot_to_VTK',
                          interface = Fn(function = annot_to_vtk,
                                         input_names = ['annot_file',
                                                        'vtk_file'],
                                         output_names = ['labels',
                                                         'output_vtk']))
        labelFlow.add_nodes([FreeLabels])
        mbFlow.connect([(Annot, labelFlow,
                         [('annot_files', 'DK_annot_to_VTK.annot_file')])])
        if do_input_vtk:
            mbFlow.connect([(Surf, labelFlow,
                             [('surface_files',
                               'DK_annot_to_VTK.vtk_file')])])
        else:
            mbFlow.connect([(ConvertSurf, labelFlow,
                             [('vtk_file',
                               'DK_annot_to_VTK.vtk_file')])])
        mbFlow.connect([(labelFlow, Sink,
                         [('DK_annot_to_VTK.output_vtk', 'labels.@DKsurface')])])
        init_labels_plug = 'DK_annot_to_VTK.output_vtk'

    #=========================================================================
    # Initialize labels with the DKT classifier atlas
    #=========================================================================
    elif init_labels == 'DKTatlas':
        #---------------------------------------------------------------------
        # Label a brain with the DKT atlas using FreeSurfer's mris_ca_label
        #---------------------------------------------------------------------
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
        labelFlow.add_nodes([Classifier])
        mbFlow.connect([(Info, labelFlow,
                         [('hemi', 'Label_with_DKTatlas.hemi'),
                          ('subject', 'Label_with_DKTatlas.subject')])])
        Classifier.inputs.subjects_path = subjects_path
        mbFlow.connect([(Surf, labelFlow,
                         [('sphere_files',
                           'Label_with_DKTatlas.sphere_file')])])
        Classifier.inputs.classifier_path = classifier_path
        Classifier.inputs.classifier_atlas = classifier_atlas

        #---------------------------------------------------------------------
        # Convert .annot file to .vtk format
        #---------------------------------------------------------------------
        Classifier2vtk = Node(name = 'DKT_annot_to_VTK',
                              interface = Fn(function = annot_to_vtk,
                                             input_names = ['annot_file',
                                                            'vtk_file'],
                                             output_names = ['labels',
                                                             'output_vtk']))
        labelFlow.add_nodes([Classifier2vtk])
        labelFlow.connect([(Classifier, Classifier2vtk,
                            [('annot_file', 'annot_file')])])
        if do_input_vtk:
            mbFlow.connect([(Surf, labelFlow,
                             [('surface_files',
                               'DKT_annot_to_VTK.vtk_file')])])
        else:
            mbFlow.connect([(ConvertSurf, labelFlow,
                             [('vtk_file',
                               'DKT_annot_to_VTK.vtk_file')])])
        #mbFlow.connect([(labelFlow, Sink,
        #                 [('Classifier2vtk.output_vtk', 'labels.@DKTsurface')])])
        init_labels_plug = 'DKT_annot_to_VTK.output_vtk'

    #=========================================================================
    # Initialize labels using multi-atlas registration
    #=========================================================================
    elif init_labels == 'max':
        #---------------------------------------------------------------------
        # Register surfaces to average template
        #---------------------------------------------------------------------
        free_template = 'OASIS-TRT-20'  # FreeSurfer template

        Register = Node(name = 'Register_template',
                        interface = Fn(function = register_template,
                                       input_names = ['hemi',
                                                      'sphere_file',
                                                      'transform',
                                                      'templates_path',
                                                      'template'],
                                       output_names = ['transform']))
        labelFlow.add_nodes([Register])
        mbFlow.connect([(Info, labelFlow, [('hemi', 'Register_template.hemi')]),
                        (Surf, labelFlow, [('sphere_files',
                                            'Register_template.sphere_file')])])
        Register.inputs.transform = 'sphere_to_' + template + '_template.reg'
        Register.inputs.templates_path = os.path.join(templates_path, 'freesurfer')
        Register.inputs.template = free_template + '.tif'
        #---------------------------------------------------------------------
        # Register atlases to subject via template
        #---------------------------------------------------------------------
        # Load atlas list
        atlas_list_file = os.path.join(x_path, 'mindboggle101_atlases.txt')
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
        labelFlow.add_nodes([Transform])
        mbFlow.connect([(Info, labelFlow,
                         [('hemi', 'Transform_labels.hemi'),
                          ('subject', 'Transform_labels.subject')])])
        labelFlow.connect([(Register, Transform, [('transform', 'transform')])])
        #Transform.inputs.transform = 'sphere_to_' + template + '_template.reg'
        Transform.inputs.subjects_path = subjects_path
        Transform.inputs.atlas = atlas_list
        Transform.inputs.atlas_string = 'labels.' + protocol + '.' + label_method
        #---------------------------------------------------------------------
        # Majority vote label
        #---------------------------------------------------------------------
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
        labelFlow.add_nodes([Vote])
        if do_input_vtk:
            mbFlow.connect([(Surf, labelFlow,
                             [('surface_files', 'Label_vote.surface_file')])])
        else:
            mbFlow.connect([(ConvertSurf, labelFlow,
                             [('vtk_file', 'Label_vote.surface_file')])])
        labelFlow.connect([(Transform, Vote, [('output_file', 'annot_files')])])
        mbFlow.connect([(labelFlow, Sink,
                         [('Label_vote.maxlabel_file', 'labels.@max'),
                          ('Label_vote.labelcounts_file', 'labels.@counts'),
                          ('Label_vote.labelvotes_file', 'labels.@votes')])])
        init_labels_plug = 'Label_vote.maxlabel_file'

    #=========================================================================
    # Skip label initialization and process manual (atlas) labels
    #=========================================================================
    elif init_labels == 'manual':
        AtlasLabels = Node(name = 'Atlas_labels',
                           interface = Fn(function = read_vtk,
                                          input_names = ['input_vtk',
                                                         'return_first',
                                                         'return_array'],
                                          output_names = ['faces',
                                                          'lines',
                                                          'indices',
                                                          'points',
                                                          'npoints',
                                                          'scalars',
                                                          'scalar_names',
                                                          'input_vtk']))
        labelFlow.add_nodes([AtlasLabels])
        mbFlow.connect([(Atlas, labelFlow,
                         [('atlas_file', 'Atlas_labels.input_vtk')])])
        AtlasLabels.inputs.return_first = 'True'
        AtlasLabels.inputs.return_array = 'False'
        init_labels_plug = 'Atlas_labels.input_vtk'

##############################################################################
#
#   Surface shape measurement workflow
#
##############################################################################
if run_shapeFlow:

    shapeFlow = Workflow(name='Shapes')

    #=========================================================================
    # Measure surface area
    #=========================================================================
    AreaNode = Node(name='Area',
                interface = Fn(function = area,
                               input_names = ['command',
                                              'surface_file'],
                               output_names = ['area_file']))
    area_command = os.path.join(ccode_path, 'area', 'PointAreaMain')
    AreaNode.inputs.command = area_command

    #=========================================================================
    # Measure surface depth
    #=========================================================================
    DepthNode = Node(name='Depth',
                     interface = Fn(function = depth,
                                    input_names = ['command',
                                                   'surface_file'],
                                    output_names = ['depth_file']))
    depth_command = os.path.join(ccode_path, 'travel_depth', 'TravelDepthMain')
    DepthNode.inputs.command = depth_command

    #=========================================================================
    # Measure surface curvature
    #=========================================================================
    CurvNode = Node(name='Curvature',
                     interface = Fn(function = curvature,
                                    input_names = ['command',
                                                   'surface_file'],
                                    output_names = ['mean_curvature_file',
                                                    'gauss_curvature_file',
                                                    'max_curvature_file',
                                                    'min_curvature_file',
                                                    'min_curvature_vector_file']))
    curvature_command = os.path.join(ccode_path, 'curvature', 'CurvatureMain') +\
                         ' -m 0 -n 2 '
    CurvNode.inputs.command = curvature_command

    #=========================================================================
    # Convert FreeSurfer surface measures to VTK
    #=========================================================================
    if do_thickness:
        ThickNode = Node(name = 'Thickness_to_VTK',
                         interface = Fn(function = curvature_to_vtk,
                                        input_names = ['surface_file',
                                                       'vtk_file'],
                                        output_names = ['output_vtk']))
        shapeFlow.add_nodes([ThickNode])
        mbFlow.connect([(Surf, shapeFlow,
                         [('thickness_files','Thickness_to_VTK.surface_file')])])
        mbFlow.connect([(ConvertSurf, shapeFlow,
                         [('vtk_file', 'Thickness_to_VTK.vtk_file')])])
        mbFlow.connect([(shapeFlow, Sink,
                         [('Thickness_to_VTK.output_vtk', 'shapes.@thickness')])])
    if do_convexity:
        ConvexNode = Node(name = 'Convexity_to_VTK',
                          interface = Fn(function = curvature_to_vtk,
                                         input_names = ['surface_file',
                                                        'vtk_file'],
                                         output_names = ['output_vtk']))
        shapeFlow.add_nodes([ConvexNode])
        mbFlow.connect([(Surf, shapeFlow,
                         [('convexity_files','Convexity_to_VTK.surface_file')])])
        mbFlow.connect([(ConvertSurf, shapeFlow,
                         [('vtk_file', 'Convexity_to_VTK.vtk_file')])])
        mbFlow.connect([(shapeFlow, Sink,
                         [('Convexity_to_VTK.output_vtk', 'shapes.@convexity')])])
    #-------------------------------------------------------------------------
    # Add and connect nodes, save output files
    #-------------------------------------------------------------------------
    shapeFlow.add_nodes([AreaNode, DepthNode, CurvNode])
    if do_input_vtk:
        mbFlow.connect([(Surf, shapeFlow,
                         [('surface_files','Area.surface_file')])])
        mbFlow.connect([(Surf, shapeFlow,
                         [('surface_files','Depth.surface_file')])])
        mbFlow.connect([(Surf, shapeFlow,
                         [('surface_files','Curvature.surface_file')])])
    else:
        mbFlow.connect([(ConvertSurf, shapeFlow,
                         [('vtk_file', 'Area.surface_file')])])
        mbFlow.connect([(ConvertSurf, shapeFlow,
                         [('vtk_file', 'Depth.surface_file')])])
        mbFlow.connect([(ConvertSurf, shapeFlow,
                         [('vtk_file', 'Curvature.surface_file')])])
    mbFlow.connect([(shapeFlow, Sink,
                     [('Area.area_file', 'shapes.@area')])])
    mbFlow.connect([(shapeFlow, Sink,
                     [('Depth.depth_file', 'shapes.@depth')])])
    mbFlow.connect([(shapeFlow, Sink,
     [('Curvature.mean_curvature_file', 'shapes.@mean_curvature'),
      ('Curvature.gauss_curvature_file', 'shapes.@gauss_curvature'),
      ('Curvature.max_curvature_file', 'shapes.@max_curvature'),
      ('Curvature.min_curvature_file', 'shapes.@min_curvature'),
      ('Curvature.min_curvature_vector_file', 'shapes.@min_curvature_vectors')])])

##############################################################################
#
#   Feature extraction workflow
#
##############################################################################
if run_featureFlow:

    featureFlow = Workflow(name='Features')
    min_fold_size = 50

    #=========================================================================
    # Folds
    #=========================================================================
    FoldsNode = Node(name='Folds',
                     interface = Fn(function = extract_folds,
                                    input_names = ['depth_file',
                                                   'min_fold_size',
                                                   'tiny_depth',
                                                   'save_file'],
                                    output_names = ['folds',
                                                    'n_folds',
                                                    'depth_threshold',
                                                    'bins',
                                                    'bin_edges',
                                                    'folds_file']))
    featureFlow.add_nodes([FoldsNode])
    mbFlow.connect([(shapeFlow, featureFlow,
                     [('Depth.depth_file','Folds.depth_file')])])
    FoldsNode.inputs.min_fold_size = min_fold_size
    FoldsNode.inputs.tiny_depth = 0.001
    FoldsNode.inputs.save_file = True
    # Save folds
    mbFlow.connect([(featureFlow, Sink,
                     [('Folds.folds_file','features.@folds')])])

    """
    # Subfolds
    SubfoldsNode = Node(name='Subfolds',
                        interface = Fn(function = extract_subfolds,
                                       input_names = ['depth_file',
                                                      'folds',
                                                      'depth_factor',
                                                      'depth_ratio',
                                                      'tolerance',
                                                      'save_file'],
                                       output_names = ['subfolds',
                                                       'n_subfolds',
                                                       'subfolds_file']))
    featureFlow.add_nodes([SubfoldsNode])
    mbFlow.connect([(shapeFlow, featureFlow,
                     [('Depth.depth_file','Subfolds.depth_file')])])
    featureFlow.connect([(FoldsNode, SubfoldsNode, [('folds','folds')])])
    SubfoldsNode.inputs.depth_factor = 0.25
    SubfoldsNode.inputs.depth_ratio = 0.1
    SubfoldsNode.inputs.tolerance = 0.01
    SubfoldsNode.inputs.save_file = True
    # Save subfolds
    mbFlow.connect([(featureFlow, Sink,
                     [('Subfolds.subfolds_file','features.@subfolds')])])
    """

    #=========================================================================
    # Rescaled travel depth
    #=========================================================================
    if run_shapeFlow:

        #=====================================================================
        # Rescale surface depth
        #=====================================================================
        RescaleDepth = Node(name='Rescale_depth',
                            interface = Fn(function = rescale_by_neighborhood,
                                           input_names = ['input_vtk',
                                                          'indices',
                                                          'nedges',
                                                          'p',
                                                          'set_max_to_1',
                                                          'save_file',
                                                          'output_filestring'],
                                           output_names = ['rescaled_scalars',
                                                           'rescaled_scalars_file']))
        shapeFlow.add_nodes([RescaleDepth])
        mbFlow.connect([(DepthNode, RescaleDepth, [('depth_file','input_vtk')])])
        RescaleDepth.inputs.indices = []
        RescaleDepth.inputs.nedges = 10
        RescaleDepth.inputs.p = 99
        RescaleDepth.inputs.set_max_to_1 = True
        RescaleDepth.inputs.save_file = True
        RescaleDepth.inputs.output_filestring = 'depth_rescaled'
        # Save rescaled depth
        mbFlow.connect([(shapeFlow, Sink,
                         [('Rescale_depth.rescaled_scalars_file','shapes.@depth_rescaled')])])

    #=========================================================================
    # Sulci
    #=========================================================================
    if do_sulci:
        LabelPairs = Node(name='Label_pairs',
                          interface = Fn(function = sulcus_boundaries,
                                         input_names = [],
                                         output_names = ['label_pair_lists']))
        featureFlow.add_nodes([LabelPairs])

        SulciNode = Node(name='Sulci',
                         interface = Fn(function = extract_sulci,
                                        input_names = ['labels_file',
                                                       'folds_or_file',
                                                       'label_pair_lists',
                                                       'min_boundary',
                                                       'sulcus_names',
                                                       'save_file'],
                                        output_names = ['sulci',
                                                        'n_sulci',
                                                        'sulci_file']))
        featureFlow.add_nodes([SulciNode])
        mbFlow.connect([(labelFlow, featureFlow,
                         [(init_labels_plug, 'Sulci.labels_file')])])
        featureFlow.connect([(FoldsNode, SulciNode, [('folds','folds_or_file')])])
        featureFlow.connect([(LabelPairs, SulciNode,
                              [('label_pair_lists','label_pair_lists')])])
        SulciNode.inputs.min_boundary = 1
        sulcus_names_file = os.path.join(protocol_path, 'sulci.names.DKT25.txt')
        fid = open(sulcus_names_file, 'r')
        sulcus_names = fid.readlines()
        sulcus_names = [x.strip('\n') for x in sulcus_names]
        SulciNode.inputs.sulcus_names = sulcus_names
        SulciNode.inputs.save_file = True
        # Save sulci
        mbFlow.connect([(featureFlow, Sink,
                         [('Sulci.sulci_file','features.@sulci')])])

    #=========================================================================
    # Fundi (curves at the bottoms of folds/sulci)
    #=========================================================================
    if do_fundi:
        LikelihoodNode = Node(name='Likelihood',
                              interface = Fn(function = compute_likelihood,
                                             input_names = ['trained_file',
                                                            'depth_file',
                                                            'curvature_file',
                                                            'folds'],
                                             output_names = ['likelihoods']))

        featureFlow.add_nodes([LikelihoodNode])
        LikelihoodNode.inputs.trained_file = os.path.join(data_path, 'atlases',
            'depth_curv_border_nonborder_parameters.pkl')
        mbFlow.connect([(shapeFlow, featureFlow,
                         [('Rescale_depth.rescaled_scalars_file',
                           'Likelihood.depth_file'),
                          ('Curvature.mean_curvature_file',
                           'Likelihood.curvature_file')])])
        featureFlow.connect([(FoldsNode, LikelihoodNode, [('folds','folds')])])

        FundiNode = Node(name='Fundi',
                         interface = Fn(function = extract_fundi,
                                        input_names = ['folds',
                                                       'sulci',
                                                       'likelihoods',
                                                       'rescaled_depth_file',
                                                       'depth_file',
                                                       'smooth_skeleton',
                                                       'save_file'],
                                        output_names = ['fundi',
                                                        'n_fundi',
                                                        'fundi_file']))
        featureFlow.connect([(FoldsNode, FundiNode, [('folds','folds')])])
        featureFlow.connect([(SulciNode, FundiNode, [('sulci','sulci')])])
        featureFlow.connect([(LikelihoodNode, FundiNode,
                              [('likelihoods', 'likelihoods')])])
        mbFlow.connect([(shapeFlow, featureFlow,
                         [('Rescale_depth.rescaled_scalars_file',
                           'Fundi.rescaled_depth_file'),
                          ('Depth.depth_file','Fundi.depth_file')])])
        FundiNode.inputs.smooth_skeleton = False
        FundiNode.inputs.save_file = True
        # Save VTK file with fundi:
        mbFlow.connect([(featureFlow, Sink,
                         [('Fundi.fundi_file','features.@fundi')])])

##############################################################################
#
#   Shape measurement workflow (continued, for features)
#
##############################################################################
if run_shapeFlow:

    if do_measure_spectra:
        #=====================================================================
        # Measure Laplace-Beltrami spectra of labeled regions
        #=====================================================================
        SpectraLabels = Node(name='Spectra_labels',
                             interface = Fn(function = fem_laplacian_from_labels,
                                            input_names = ['vtk_file',
                                                           'n_eigenvalues',
                                                           'normalization'],
                                            output_names = ['spectrum_lists',
                                                            'label_list']))
        shapeFlow.add_nodes([SpectraLabels])
        mbFlow.connect([(labelFlow, shapeFlow,
                         [(init_labels_plug, 'Spectra_labels.vtk_file')])])
        SpectraLabels.inputs.n_eigenvalues = 6
        SpectraLabels.inputs.normalization = "area"

        #=====================================================================
        # Measure Laplace-Beltrami spectra of sulci
        #=====================================================================
        if do_sulci:
            SpectraSulci = SpectraLabels.clone('Spectra_sulci')
            shapeFlow.add_nodes([SpectraSulci])
            mbFlow.connect([(SulciNode, SpectraSulci,
                             [('sulci_file', 'vtk_file')])])

##############################################################################
#
#   Table construction workflow
#
##############################################################################
if run_tableFlow:

    tableFlow = Workflow(name='Tables')

    #=========================================================================
    # Shape tables of surface: labels, fundi, and sulci
    #=========================================================================
    ShapeTables = Node(name='Shape_tables',
                       interface = Fn(function = write_mean_shapes_tables,
                                      input_names = ['labels_or_file',
                                                     'sulci',
                                                     'fundi',
                                                     'area_file',
                                                     'depth_file',
                                                     'mean_curvature_file',
                                                     'gauss_curvature_file',
                                                     'max_curvature_file',
                                                     'min_curvature_file',
                                                     'thickness_file',
                                                     'convexity_file',
                                                     'labels_spectra',
                                                     'labels_spectra_norm',
                                                     'labels_spectra_IDs',
                                                     'sulci_spectra',
                                                     'sulci_spectra_norm',
                                                     'sulci_spectra_IDs',
                                                     'exclude_labels',
                                                     'delimiter'],
                                      output_names = ['label_table',
                                                      'fundus_table',
                                                      'sulcus_table',
                                                      'norm_label_table',
                                                      'norm_fundus_table',
                                                      'norm_sulcus_table']))
    tableFlow.add_nodes([ShapeTables])
    mbFlow.connect([(labelFlow, tableFlow,
                     [(init_labels_plug, 'Shape_tables.labels_or_file')])])
    mbFlow.connect([(featureFlow, tableFlow,
                     [('Sulci.sulci', 'Shape_tables.sulci')])])
    if do_fundi:
        mbFlow.connect([(featureFlow, tableFlow,
                         [('Fundi.fundi', 'Shape_tables.fundi')])])
    else:
        ShapeTables.inputs.fundi = []
    #-------------------------------------------------------------------------
    mbFlow.connect([(shapeFlow, tableFlow,
                     [('Area.area_file','Shape_tables.area_file')])])
    mbFlow.connect([(shapeFlow, tableFlow,
                     [('Depth.depth_file','Shape_tables.depth_file')])])
    mbFlow.connect([(shapeFlow, tableFlow,
                     [('Curvature.mean_curvature_file',
                       'Shape_tables.mean_curvature_file')])])
    mbFlow.connect([(shapeFlow, tableFlow,
                     [('Curvature.gauss_curvature_file',
                       'Shape_tables.gauss_curvature_file')])])
    mbFlow.connect([(shapeFlow, tableFlow,
                     [('Curvature.max_curvature_file',
                       'Shape_tables.max_curvature_file')])])
    mbFlow.connect([(shapeFlow, tableFlow,
                     [('Curvature.min_curvature_file',
                       'Shape_tables.min_curvature_file')])])
    if do_thickness:
        mbFlow.connect([(shapeFlow, tableFlow,
                         [('Thickness_to_VTK.output_vtk',
                           'Shape_tables.thickness_file')])])
    else:
        ShapeTables.inputs.thickness_file = ''
    if do_convexity:
        mbFlow.connect([(shapeFlow, tableFlow,
                         [('Convexity_to_VTK.output_vtk',
                           'Shape_tables.convexity_file')])])
    else:
        ShapeTables.inputs.convexity_file = ''
    if do_measure_spectra:
        mbFlow.connect([(shapeFlow, tableFlow,
                         [('Spectra_labels.spectrum_lists',
                           'Shape_tables.labels_spectra')])])
        if do_sulci:
            mbFlow.connect([(shapeFlow, tableFlow,
                             [('Spectra_sulci.spectrum_lists',
                               'Shape_tables.sulci_spectra')])])
        else:
            ShapeTables.inputs.sulci_spectra = []
    else:
        ShapeTables.inputs.labels_spectra = []
        ShapeTables.inputs.sulci_spectra = []

    #-------------------------------------------------------------------------
    ShapeTables.inputs.exclude_labels = [-1]
    ShapeTables.inputs.delimiter = ","
    # Save results
    mbFlow.connect([(tableFlow, Sink,
                     [('Shape_tables.label_table', 'tables.@labels'),
                      ('Shape_tables.sulcus_table', 'tables.@sulci'),
                      ('Shape_tables.fundus_table', 'tables.@fundi'),
                      ('Shape_tables.norm_label_table', 'tables.@labels_norm'),
                      ('Shape_tables.norm_fundus_table', 'tables.@fundi_norm'),
                      ('Shape_tables.norm_sulcus_table', 'tables.@sulci_norm')])])

    #=========================================================================
    # Per-vertex shapes
    #=========================================================================
    if do_vertex_tables:

        VertexTable = Node(name='Vertex_table',
                           interface = Fn(function = write_vertex_shapes_table,
                                          input_names = ['table_file',
                                                         'labels_or_file',
                                                         'sulci',
                                                         'fundi',
                                                         'area_file',
                                                         'depth_file',
                                                         'depth_rescaled_file',
                                                         'mean_curvature_file',
                                                         'gauss_curvature_file',
                                                         'max_curvature_file',
                                                         'min_curvature_file',
                                                         'thickness_file',
                                                         'convexity_file',
                                                         'delimiter'],
                                          output_names = ['shapes_table']))
        tableFlow.add_nodes([VertexTable])
        VertexTable.inputs.table_file = 'vertex_shapes.csv'
        mbFlow.connect([(labelFlow, tableFlow,
                         [(init_labels_plug, 'Vertex_table.labels_or_file')])])
        mbFlow.connect([(featureFlow, tableFlow,
                         [('Sulci.sulci', 'Vertex_table.sulci')])])
        if do_fundi:
            mbFlow.connect([(featureFlow, tableFlow,
                             [('Fundi.fundi', 'Vertex_table.fundi')])])
        else:
            ShapeTables.inputs.fundi = []
        #---------------------------------------------------------------------
        mbFlow.connect([(shapeFlow, tableFlow,
                         [('Area.area_file','Vertex_table.area_file')])])
        mbFlow.connect([(shapeFlow, tableFlow,
                         [('Depth.depth_file','Vertex_table.depth_file')])])
        mbFlow.connect([(shapeFlow, tableFlow,
                         [('Rescale_depth.rescaled_scalars_file',
                           'Vertex_table.depth_rescaled_file')])])
        mbFlow.connect([(shapeFlow, tableFlow,
                         [('Curvature.mean_curvature_file',
                           'Vertex_table.mean_curvature_file')])])
        mbFlow.connect([(shapeFlow, tableFlow,
                         [('Curvature.gauss_curvature_file',
                           'Vertex_table.gauss_curvature_file')])])
        mbFlow.connect([(shapeFlow, tableFlow,
                         [('Curvature.max_curvature_file',
                           'Vertex_table.max_curvature_file')])])
        mbFlow.connect([(shapeFlow, tableFlow,
                         [('Curvature.min_curvature_file',
                           'Vertex_table.min_curvature_file')])])
        if do_thickness:
            mbFlow.connect([(shapeFlow, tableFlow,
                             [('Thickness_to_VTK.output_vtk',
                               'Vertex_table.thickness_file')])])
        else:
            VertexTable.inputs.thickness_file = ''
        if do_convexity:
            mbFlow.connect([(shapeFlow, tableFlow,
                             [('Convexity_to_VTK.output_vtk',
                               'Vertex_table.convexity_file')])])
        else:
            VertexTable.inputs.convexity_file = ''
        #---------------------------------------------------------------------
        VertexTable.inputs.delimiter = ","
        mbFlow.connect([(tableFlow, Sink,
                         [('Vertex_table.shapes_table',
                           'tables.@vertex_table')])])

#=============================================================================
# Surface label evaluation
#=============================================================================
if do_evaluate_surface:

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
    mbFlow.connect([(labelFlow, EvalSurfLabels,
                     [(init_labels_plug, 'labels_file2')])])
    mbFlow.connect([(EvalSurfLabels, Sink,
                     [('overlap_file', 'evaluate_labels')])])

##############################################################################
#
#   Fill volume prep workflow:
#   Convert labels from VTK to .annot format
#
##############################################################################
if run_volumeFlow and do_fill:

    annotflow = Workflow(name='Fill_volume_prep')

    #=========================================================================
    # Convert VTK labels to .annot format.
    #=========================================================================
    #-------------------------------------------------------------------------
    # Write .label files for surface vertices
    #-------------------------------------------------------------------------
    WriteLabels = Node(name='Write_label_files',
                       interface = Fn(function = vtk_to_labels,
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
    mbFlow.connect([(labelFlow, annotflow,
                     [(init_labels_plug, 'Write_label_files.surface_file')])])
    if init_labels == 'max':
        WriteLabels.inputs.scalar_name = 'Max_(majority_labels)'
    else:
        WriteLabels.inputs.scalar_name = 'Labels'
    #-------------------------------------------------------------------------
    # Write .annot file from .label files
    # NOTE:  incorrect labels to be corrected below!
    #-------------------------------------------------------------------------
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
    WriteAnnot.inputs.annot_name = 'labels.' + protocol + '.' + init_labels
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

##############################################################################
#
#   Label volumes workflow:
#   * Fill volume
#   * Measure label volumes
#   * Evaluate volume labels
#
##############################################################################
if run_volumeFlow:

    mbFlow2 = Workflow(name='Label_volumes')
    mbFlow2.base_dir = temp_path

    #=========================================================================
    # Fill (gray matter) volume using FreeSurfer
    #=========================================================================
    #-------------------------------------------------------------------------
    # Iterate inputs over subjects
    #-------------------------------------------------------------------------
    Info2 = Info.clone('Inputs2')
    Info2.iterables = ([('subject', subjects)])
    Sink2 = Sink.clone('Results2')

    #-------------------------------------------------------------------------
    # Fill volume mask with surface vertex labels from .annot file.
    # Convert label volume from FreeSurfer 'unconformed' to original space.
    #-------------------------------------------------------------------------
    if do_fill:

        FillVolume = Node(name='Fill_volume',
                          interface = Fn(function = labels_to_volume,
                                         input_names = ['subject',
                                                        'annot_name',
                                                        'original_space',
                                                        'reference'],
                                         output_names = ['output_file']))
        mbFlow2.add_nodes([FillVolume])
        mbFlow2.connect([(Info2, FillVolume, [('subject', 'subject')])])
        FillVolume.inputs.annot_name = 'labels.' + protocol + '.' + init_labels
        FillVolume.inputs.original_space = True
        mbFlow2.connect([(Info2, Vol, [('subject','subject')])])
        mbFlow2.connect([(Vol, FillVolume, [('original_volume', 'reference')])])
        #---------------------------------------------------------------------
        # Relabel file, replacing colortable labels with real labels
        #---------------------------------------------------------------------
        Relabel = Node(name='Correct_labels',
                       interface = Fn(function = relabel_volume,
                                      input_names = ['input_file',
                                                     'old_labels',
                                                     'new_labels'],
                                      output_names = ['output_file']))
        mbFlow2.add_nodes([Relabel])
        mbFlow2.connect([(FillVolume, Relabel, [('output_file', 'input_file')])])
        relabel_file = os.path.join(protocol_path,
                            'labels.volume.annot_errors.' + protocol + '.txt')
        old_labels, new_labels = read_columns(relabel_file, 2)
        Relabel.inputs.old_labels = old_labels
        Relabel.inputs.new_labels = new_labels
        mbFlow2.connect([(Relabel, Sink2, [('output_file', 'labels_volume')])])

    #=========================================================================
    # Compute volume per label.
    #=========================================================================
    if do_measure_volume:

        #---------------------------------------------------------------------
        # Measure volume of each region of a labeled image file.
        #---------------------------------------------------------------------
        MeasureVolumes = Node(name='Measure_volumes',
                              interface = Fn(function = volume_per_label,
                                             input_names = ['labels',
                                                            'input_file'],
                                             output_names = ['volumes',
                                                             'labels']))
        mbFlow2.add_nodes([MeasureVolumes])
        volume_labels_list_file = os.path.join(protocol_path,
                                               'labels.volume.'+protocol+'.txt')
        volume_labels_list = read_columns(volume_labels_list_file, 1)[0]
        volume_labels_list = [int(x) for x in volume_labels_list]
        MeasureVolumes.inputs.labels = volume_labels_list
        if do_fill:
            mbFlow2.connect([(Relabel, MeasureVolumes,
                              [('output_file', 'input_file')])])
        else:
            sys.exit('No alternative set of label volumes provided...')

        #---------------------------------------------------------------------
        # Create a table to save the volume measures.
        #---------------------------------------------------------------------
        InitVolTable = Node(name='Initialize_Volume_label_table',
                            interface = Fn(function = write_columns,
                                           input_names = ['columns',
                                                          'column_names',
                                                          'output_table',
                                                          'delimiter',
                                                          'input_table'],
                                           output_names = ['output_table']))
        mbFlow2.add_nodes([InitVolTable])
        mbFlow2.connect([(MeasureVolumes, InitVolTable,
                          [('labels', 'columns')])])
        InitVolTable.inputs.column_names = ['label']
        InitVolTable.inputs.output_table = 'volume_labels.csv'
        InitVolTable.inputs.delimiter = ','
        InitVolTable.inputs.input_table = ''

        VolumeLabelTable = Node(name='Volume_label_table',
                                interface = Fn(function = write_columns,
                                               input_names = ['columns',
                                                              'column_names',
                                                              'output_table',
                                                              'delimiter',
                                                              'input_table'],
                                               output_names = ['output_table']))
        mbFlow2.connect([(MeasureVolumes, VolumeLabelTable,
                          [('volumes', 'columns')])])
        VolumeLabelTable.inputs.column_names = ['volume']
        VolumeLabelTable.inputs.output_table = 'label_volume_shapes.csv'
        VolumeLabelTable.inputs.delimiter = ','
        mbFlow2.connect([(InitVolTable, VolumeLabelTable,
                          [('output_table', 'input_table')])])
        # Save table of label volumes
        mbFlow2.connect([(VolumeLabelTable, Sink2,
                          [('output_table', 'tables.@volume_labels')])])

    #=========================================================================
    # Evaluate label volume overlaps.
    #=========================================================================
    if do_evaluate_volume:

        #---------------------------------------------------------------------
        # Evaluation inputs: location and structure of atlas volumes
        #---------------------------------------------------------------------
        AtlasVol = Node(name = 'Atlas_volume',
                         interface = DataGrabber(infields=['subject'],
                         outfields=['atlas_vol_file']))
        AtlasVol.inputs.base_directory = atlases_path
        AtlasVol.inputs.template = '%s/mri/labels.' + protocol + '.manual.nii.gz'
        AtlasVol.inputs.template_args['atlas_vol_file'] = [['subject']]
        mbFlow2.connect([(Info2, AtlasVol, [('subject','subject')])])
        #---------------------------------------------------------------------
        # Evaluate volume labels
        #---------------------------------------------------------------------
        EvalVolLabels = Node(name='Evaluate_volume_labels',
                               interface = Fn(function = measure_volume_overlap,
                                              input_names = ['labels',
                                                             'file2',
                                                             'file1'],
                                              output_names = ['overlaps',
                                                              'out_file']))
        labels_file = os.path.join(protocol_path, 'labels.volume.' + protocol + '.txt')
        labels = read_columns(labels_file, 1)[0]
        EvalVolLabels.inputs.labels = labels
        mbFlow2.add_nodes([EvalVolLabels])
        mbFlow2.connect([(AtlasVol, EvalVolLabels,
                          [('atlas_vol_file','file2')])])
        mbFlow2.connect([(Relabel, EvalVolLabels,
                          [('output_file', 'file1')])])
        mbFlow2.connect([(EvalVolLabels, Sink2,
                          [('out_file', 'evaluate_labels_volume')])])

##############################################################################
#
#    Run workflows
#
##############################################################################
if __name__== '__main__':

    #from nipype import config, logging
    #config.set('logging', 'interface_level', 'DEBUG')
    #config.set('logging', 'workflow_level', 'DEBUG')
    #logging.update_logging(config)

    run_flow1 = True
    run_flow2 = True
    generate_graphs = 0#True
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

"""
import os
from mindboggle.utils.io_file import read_columns

out_path = '/Users/arno/Data/Mindboggle-101/'
x_path = os.path.join(os.environ['MINDBOGGLE'], 'x')
atlas_list_file = os.path.join(x_path, 'mindboggle101_atlases.txt')
atlas_list = read_columns(atlas_list_file, 1)[0]

for atlas in atlas_list:
    if 'MMRR-21' in atlas:
        cmd = ' '.join(['python pipeline.py', out_path, atlas])
        print(cmd); os.system(cmd)
"""