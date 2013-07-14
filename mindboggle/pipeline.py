#!/usr/bin/env python
"""
This is Mindboggle's nipype software pipeline!

Examples
--------
$ python pipeline.py <output path> <1 or more subject names>
$ python pipeline.py output HLN-12-1 HLN-12-2

..SURFACE workflows ::

      ======================
    * Surface label workflow *
      ======================
      - DKT classifier (default option)
          This option gives the same cortical labels as FreeSurfer
          (version 5.2 or greater) if the DKT-40 classifier is used
          instead of the DKT-100 classifier, and the DKT31 protocol labels
          are not converted to DKT25 protocol labels.
      - FreeSurfer labels
          Version 5.2 or greater recommended (see above).
      - Multi-atlas labeling
          If surface template and labeled surface atlases are supplied,
          this option uses FreeSurfer registration and
          majority vote rule on multiple label assignments.
      - Manual labels
          If labeled surfaces are supplied (such as the Mindboggle-101 set)
          these labels may also be used to evaluate any of the above labels.

      ==================================
    * Surface shape measurement workflow *
      ==================================
      - Surface area
      - Travel depth
      - Geodesic depth
      - Mean curvature
      - Convexity (FreeSurfer)
      - Thickness (FreeSurfer)

      ===================================
    * Surface feature extraction workflow *
      ===================================
      - Folds
      - Sulci
      - Fundi

      ==============================
    * Surface feature shape workflow *
      ==============================
      - Laplace-Beltrami spectra

      ==========================
    * Label volume prep workflow *
      ==========================
      - Label format conversion

..VOLUME workflows ::

      =====================
    * Volume label workflow *
      =====================
      - Fill cortical gray matter with labels
      - Register image volume to template in MNI152 space
      - Label subcortical volumes

      =============================
    * Volume feature shape workflow *
      =============================
      - Find positions (native and MNI152 spaces)
      - Measure label volumes
      - (Evaluate label volumes vs. manual labels)

      ==============
    * Table workflow *
      ==============
      - Volume shape tables:
          - Label volumes
      - Surface feature shape tables:
          - Label shapes
          - Sulcus shapes
          - Fundus shapes
      - Vertex measures table


.. Note::
      Mindboggle currently uses FreeSurfer for label initialization
      (its label output or its surface registration algorithm),
      and assumes that input files reside within a directory structure
      like that created by FreeSurfer (autorecon -all).
      For example, each scan is assigned a unique subject name that is also
      the name of a folder within FreeSurfer's subjects directory; within this
      folder are the subfolders surf, label, and mri.

For more information about Mindboggle:
http://mindboggle.info/software/documentation.html

For information on Nipype (http://www.nipy.org/nipype/):
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3159964/


Authors:
    - Arno Klein, 2011-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Each file lists Mindboggle team members who contributed to its content.

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
# Workflow options
#=============================================================================
#-----------------------------------------------------------------------------
run_RegFlows = True  # Run Mindboggle's registration workflows
#-----------------------------------------------------------------------------
do_register_standard = True  # Register volume to standard template
vol_reg_method = 'ANTS' # Volume registration: 'antsRegister', 'ANTS', 'flirt'
standard_template = 'OASIS-TRT-20_template_to_MNI152.nii.gz'
    # 'MNI152_T1_1mm_brain.nii.gz'

#-----------------------------------------------------------------------------
run_SurfFlows = True  # Run Mindboggle's surface workflows
run_VolFlows = True  # Run Mindboggle's volume workflows
#-----------------------------------------------------------------------------
do_input_vtk = False  # Load VTK surfaces directly (not FreeSurfer surfaces)
do_input_nifti = False  # Load nifti directly (not FreeSurfer mgh file)
do_surf_table = 0#True  # Store surface feature shape measures in a table
do_vertex_table = 0#True  # Create per-vertex shape table
do_vol_table = 0#True  # Store volume feature shape measures in a table

#-----------------------------------------------------------------------------
run_SurfLabelFlow = True
#-----------------------------------------------------------------------------
# Initialize labels with:
# * 'DKT_atlas': FreeSurfer-style classifier atlas trained on the DKT protocol
# * 'FreeSurfer': FreeSurfer (with atlas trained on the DK or DKT protocol)
# * 'max_prob': majority vote labels from multiple atlases (DISABLED)
# * 'manual': process manual labels (individual atlas)
init_labels = 'manual'
classifier_atlas = 'DKTatlas40.gcs'  # DKT_atlas: 'DKTatlas[40,100].gcs'
#free_template = 'OASIS-TRT-20'  # max_prob (FreeSurfer .tif) surface template
#atlas_list = read_columns('mindboggle101_atlases.txt', 1)[0]
#
# Labeling protocol used by Mindboggle:
# * 'DKT31': 'Desikan-Killiany-Tourville (DKT) protocol with 31 label regions
# * 'DKT25': 'fundus-friendly' version of the DKT protocol following fundi
protocol = 'DKT31'
#
# Type of atlas labels:
# * 'manual': manual edits
# * FUTURE: <'adjusted': manual edits after automated alignment to fundi>
atlas_label_type = 'manual'
#
do_evaluate_surf_labels = False  # Surface overlap: auto vs. manual labels

#-----------------------------------------------------------------------------
run_WholeSurfShapeFlow = True
#-----------------------------------------------------------------------------
do_thickness = True  # Include FreeSurfer's thickness measure
do_convexity = True  # Include FreeSurfer's convexity measure (sulc.pial)
do_measure_spectra = True  # Measure Laplace-Beltrami spectra for features

#-----------------------------------------------------------------------------
run_SurfFeatureFlow = True
#-----------------------------------------------------------------------------
do_sulci = True  # Extract sulci
do_fundi = False  # Extract fundi

#-----------------------------------------------------------------------------
run_SurfShapeFlow = True
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
run_VolLabelFlow = True
#-----------------------------------------------------------------------------
do_fill_cortex = True  # Fill cortical gray matter with surface labels
do_input_mask = False
do_label_whole_volume = 0#True  # Label whole brain (not just cortex)
atlas_volume = 'OASIS-TRT-20_atlas_to_MNI152.nii.gz'
do_evaluate_vol_labels = False  # Volume overlap: auto vs. manual labels

#-----------------------------------------------------------------------------
run_VolShapeFlow = 0#True
#-----------------------------------------------------------------------------

#=============================================================================
# Setup: import libraries, set file paths
#=============================================================================
#-----------------------------------------------------------------------------
# Import system and Nipype Python libraries
#-----------------------------------------------------------------------------
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import Function as Fn
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import DataGrabber, DataSink
from nipype.interfaces.freesurfer import MRIConvert
#from nipype.interfaces.fsl import FLIRT
#from nipype.interfaces.ants import Registration
#-----------------------------------------------------------------------------
# Import Mindboggle Python libraries
#-----------------------------------------------------------------------------
from mindboggle.utils.io_vtk import rewrite_scalars, read_vtk, \
    apply_affine_transform
from mindboggle.utils.io_table import read_columns, write_columns, \
    write_shape_stats, write_vertex_measures
from mindboggle.utils.io_free import surface_to_vtk, curvature_to_vtk, \
    annot_to_vtk
from mindboggle.utils.ants import ANTS, WarpImageMultiTransform, \
    fill_volume_with_surface_labels
from mindboggle.labels.protocol import dkt_protocol
from mindboggle.labels.label_free import label_with_classifier
from mindboggle.labels.relabel import relabel_volume, relabel_surface, \
    overwrite_volume_labels
from mindboggle.shapes.measure import area, travel_depth, geodesic_depth, \
    curvature, volume_per_label, rescale_by_neighborhood
from mindboggle.shapes.laplace_beltrami import fem_laplacian_from_labels
from mindboggle.shapes.likelihood import compute_likelihood
from mindboggle.features.folds import extract_folds
from mindboggle.features.fundi import extract_fundi
from mindboggle.features.sulci import extract_sulci
from mindboggle.evaluate.evaluate_labels import measure_surface_overlap, \
    measure_volume_overlap
#-----------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------
subjects_path = os.environ['SUBJECTS_DIR']  # FreeSurfer subjects directory
data_path = os.environ['MINDBOGGLE_DATA']  # Mindboggle data directory
ccode_path = os.environ['MINDBOGGLE_TOOLS']  # Mindboggle C++ code directory

# Name of volume template for transforming data to a standard MNI152 space:
atlas_path = os.path.join(data_path, 'atlases')
volume_template = os.path.join(atlas_path, standard_template)

# Create output directories:
temp_path = os.path.join(output_path, 'workspace')  # Where to save temp files
if not os.path.isdir(output_path):
    os.makedirs(output_path)
if not os.path.isdir(temp_path):
    os.makedirs(temp_path)
#-----------------------------------------------------------------------------
# Protocol information
#-----------------------------------------------------------------------------
sulcus_names, sulcus_label_pair_lists, unique_sulcus_label_pairs, \
    label_names, label_numbers, cortex_names, cortex_numbers, \
    noncortex_names, noncortex_numbers = dkt_protocol(protocol)

#=============================================================================
# Initialize all workflow inputs and outputs
#=============================================================================
mbFlow = Workflow(name='Mindboggle_workflow')
mbFlow.base_dir = temp_path
#-----------------------------------------------------------------------------
# Iterate inputs over subjects, hemispheres
# (surfaces are assumed to take the form: lh.pial or lh.pial.vtk)
#-----------------------------------------------------------------------------
Info = Node(name='Inputs',
            interface=IdentityInterface(fields=['subject', 'hemi']))
Info.iterables = ([('subject', subjects), ('hemi', ['lh','rh'])])
#-------------------------------------------------------------------------
# Outputs
#-------------------------------------------------------------------------
Sink = Node(DataSink(), name='Results')
Sink.inputs.base_directory = output_path
Sink.inputs.container = 'results'

if run_SurfFlows:
    #-------------------------------------------------------------------------
    # Location and structure of the surface inputs
    #-------------------------------------------------------------------------
    Surf = Node(name='Surfaces',
                interface=DataGrabber(infields=['subject', 'hemi'],
                                      outfields=['surface_files',
                                                 'sphere_files'],
                                      sort_filelist=False))
    Surf.inputs.base_directory = subjects_path
    Surf.inputs.template = '%s/surf/%s.%s'
    Surf.inputs.template_args['surface_files'] = [['subject','hemi','pial']]
    Surf.inputs.template_args['sphere_files'] = [['subject','hemi','sphere']]
    if do_thickness:
        Surf.inputs.template_args['thickness_files'] = \
            [['subject','hemi','thickness']]
    if do_convexity:
        Surf.inputs.template_args['convexity_files'] = \
            [['subject','hemi','sulc']]
    #-------------------------------------------------------------------------
    # Location and structure of the FreeSurfer label inputs
    #-------------------------------------------------------------------------
    if init_labels == 'FreeSurfer':
        Annot = Node(name='Annots',
                     interface=DataGrabber(infields=['subject', 'hemi'],
                                             outfields=['annot_files'],
                                             sort_filelist=False))
        Annot.inputs.base_directory = subjects_path
        Annot.inputs.template = '%s/label/%s.aparc.annot'
        Annot.inputs.template_args['annot_files'] = [['subject','hemi']]

if run_VolFlows or run_RegFlows:
    #-------------------------------------------------------------------------
    # Location and structure of the volume inputs
    #-------------------------------------------------------------------------
    if do_input_nifti:
        niftiBrain = Node(name='nifti',
                          interface=DataGrabber(infields=['subject'],
                                                outfields=['nifti'],
                                                sort_filelist=False))
        niftiBrain.inputs.base_directory = subjects_path
        niftiBrain.inputs.template = '%s/mri/brain.nii.gz'
        niftiBrain.inputs.template_args['nifti'] = [['subject']]
        mbFlow.connect(Info, 'subject', niftiBrain, 'subject')
    else:
        mghBrain = Node(name='mgh',
                        interface=DataGrabber(infields=['subject'],
                                              outfields=['mgh'],
                                              sort_filelist=False))
        mghBrain.inputs.base_directory = subjects_path
        mghBrain.inputs.template = '%s/mri/brain.mgz'
        mghBrain.inputs.template_args['mgh'] = [['subject']]
        mbFlow.connect(Info, 'subject', mghBrain, 'subject')

        # FreeSurfer mgh of original for converting from conformal (below):
        mghOrig = Node(name='mgh_orig',
                       interface=DataGrabber(infields=['subject'],
                                             outfields=['mgh_orig'],
                                             sort_filelist=False))
        mghOrig.inputs.base_directory = subjects_path
        mghOrig.inputs.template = '%s/mri/orig/001.mgz'
        mghOrig.inputs.template_args['mgh_orig'] = [['subject']]
        mbFlow.connect(Info, 'subject', mghOrig, 'subject')

        # Convert FreeSurfer mgh conformal brain volume to nifti format:        #---------------------------------------------------------------------
        mgh2nifti = Node(name='mgh_to_nifti', interface=MRIConvert())
        mbFlow.add_nodes([mgh2nifti])
        mbFlow.connect(mghBrain, 'mgh', mgh2nifti, 'in_file')
        mbFlow.connect(mghOrig, 'mgh_orig', mgh2nifti, 'reslice_like')
        mgh2nifti.inputs.resample_type = 'interpolate'
        mgh2nifti.inputs.out_type = 'niigz'
        mgh2nifti.inputs.out_file = 'brain.nii.gz'
        mbFlow.connect(mgh2nifti, 'out_file', Sink, 'nifti')

if run_VolFlows and run_VolLabelFlow and do_fill_cortex:
    if do_input_mask:
        niftiMask = Node(name='nifti_mask',
                         interface=DataGrabber(infields=['subject','hemi'],
                                               outfields=['nifti_mask'],
                                               sort_filelist=False))
        niftiMask.inputs.base_directory = subjects_path
        niftiMask.inputs.template = '%s/mri/%s.ribbon.nii.gz'
        niftiMask.inputs.template_args['nifti_mask'] = [['subject','hemi']]
    else:
        mghMask = Node(name='mgh_mask',
                       interface=DataGrabber(infields=['subject','hemi'],
                                             outfields=['mgh_mask'],
                                             sort_filelist=False))
        mghMask.inputs.base_directory = subjects_path
        mghMask.inputs.template = '%s/mri/%s.ribbon.mgz'
        mghMask.inputs.template_args['mgh_mask'] = [['subject','hemi']]

        # Convert FreeSurfer mgh conformal gray matter mask to nifti format:
        mgh_mask2nifti = Node(name='mgh_mask_to_nifti', interface=MRIConvert())
        mbFlow.add_nodes([mgh_mask2nifti])
        mbFlow.connect(mghMask, 'mgh_mask', mgh_mask2nifti, 'in_file')
        mbFlow.connect(mghOrig, 'mgh_orig', mgh_mask2nifti, 'reslice_like')
        mgh_mask2nifti.inputs.resample_type = 'nearest'
        mgh_mask2nifti.inputs.out_type = 'niigz'
        mgh_mask2nifti.inputs.out_file = 'mask.nii.gz'
        mbFlow.connect(mgh_mask2nifti, 'out_file', Sink, 'nifti_mask')


#=============================================================================
##############################################################################
#=============================================================================
#
#
#   Registration workflows
#
#       - Register image volume to template in MNI152 space
#
#=============================================================================
##############################################################################
#=============================================================================

if run_RegFlows:

    #=========================================================================
    # Register image volume to template in MNI152 space
    #=========================================================================
    if do_register_standard:

#        #---------------------------------------------------------------------
#        # Register image volume to template in MNI152 space using antsRegister:
#        #---------------------------------------------------------------------
#         if vol_reg_method == 'antsRegister':
#             regAnts = Node(name='antsRegister_standard', interface=Registration())
#             mbFlow.add_nodes([regAnts])
#             if do_input_nifti:
#                 mbFlow.connect(niftiBrain, 'nifti',
#                                 regAnts, 'moving_image')
#             else:
#                 mbFlow.connect(mgh2nifti, 'out_file',
#                                 regAnts, 'moving_image')
#             regAnts.inputs.fixed_image = [volume_template]
#             regAnts.inputs.num_threads = 2
#             regAnts.inputs.winsorize_lower_quantile = 0.01
#             regAnts.inputs.winsorize_upper_quantile = 0.99
#             regAnts.inputs.output_warped_image = False
#             regAnts.inputs.dimension = 3
#             regAnts.inputs.write_composite_transform = True
#             regAnts.inputs.collapse_output_transforms = True
#             regAnts.inputs.write_composite_transform = True
#             regAnts.inputs.output_transform_prefix = 'standard_'
#             if do_label_whole_volume:
#                 regAnts.inputs.transforms = ['Rigid', 'Affine', 'SyN']
#                 regAnts.inputs.transform_parameters = [(0.1,), (0.1,), (0.1, 3.0, 0.0)]
#                 regAnts.inputs.number_of_iterations = [[1000,500,250,100]]*2 + [[100,100,70,20]]
# #                regAnts.inputs.number_of_iterations = [[10,5,2,1]]*2 + [[1,1,7,2]]
#                 regAnts.inputs.metric = ['MI']*2 + ['CC']
#                 regAnts.inputs.metric_weight = [1]*3
#                 regAnts.inputs.radius_or_number_of_bins = [32]*2 + [4]
#                 regAnts.inputs.sampling_strategy = ['Regular']*2 + [None]
#                 regAnts.inputs.sampling_percentage = [0.25]*2 + [None]
#                 regAnts.inputs.convergence_threshold = [1.e-8]*2 + [1e-9]
#                 regAnts.inputs.convergence_window_size = [10]*2 + [15]
#                 regAnts.inputs.smoothing_sigmas = [[3,2,1,0]]*3
#                 regAnts.inputs.sigma_units = ['mm']*3
#                 regAnts.inputs.shrink_factors = [[8,4,2,1]]*2 + [[6,4,2,1]]
#                 regAnts.inputs.use_estimate_learning_rate_once = [True, True, True]
#                 regAnts.inputs.use_histogram_matching = [False]*2 + [True]
#                 regAnts.inputs.initial_moving_transform_com = True
#             else:
#                 regAnts.inputs.transforms = ['Rigid', 'Affine']
#                 regAnts.inputs.transform_parameters = [(0.1,), (0.1,)]
# #                regAnts.inputs.number_of_iterations = [[10,5,2,1]]*2
#                 regAnts.inputs.number_of_iterations = [[1000,500,250,100]]*2
#                 regAnts.inputs.metric = ['MI']*2
#                 regAnts.inputs.metric_weight = [1]*2
#                 regAnts.inputs.radius_or_number_of_bins = [32]*2
#                 regAnts.inputs.sampling_strategy = ['Regular']*2
#                 regAnts.inputs.sampling_percentage = [0.25]*2
#                 regAnts.inputs.convergence_threshold = [1.e-8]*2
#                 regAnts.inputs.convergence_window_size = [10]*2
#                 regAnts.inputs.smoothing_sigmas = [[3,2,1,0]]*2
#                 regAnts.inputs.sigma_units = ['mm']*2
#                 regAnts.inputs.shrink_factors = [[8,4,2,1]]*2
#                 regAnts.inputs.use_estimate_learning_rate_once = [True, True]
#                 regAnts.inputs.use_histogram_matching = [False]*2
#             mbFlow.connect(regAnts, 'composite_transform',
#                            Sink, 'transforms.@affine_antsRegistration')
        #---------------------------------------------------------------------
        # Register image volume to template in MNI152 space using ANTS:
        #---------------------------------------------------------------------
        if vol_reg_method == 'ANTS':
            regANTS = Node(name='ANTS',
                           interface=Fn(function = ANTS,
                                        input_names=['source',
                                                     'target',
                                                     'iterations',
                                                     'output_stem'],
                                        output_names=['affine_transform',
                                                      'nonlinear_transform',
                                                      'nonlinear_inverse_transform',
                                                      'output_stem']))
            mbFlow.add_nodes([regANTS])
            if do_input_nifti:
                mbFlow.connect(niftiBrain, 'nifti',
                               regANTS, 'source')
            else:
                mbFlow.connect(mgh2nifti, 'out_file',
                               regANTS, 'source')
            regANTS.inputs.target = volume_template
            if do_label_whole_volume:
                regANTS.inputs.iterations = '1' #'33x99x11'
            else:
                regANTS.inputs.iterations = '0'
            regANTS.inputs.output_stem = ''
            mbFlow.connect(regANTS, 'affine_transform',
                           Sink, 'transforms.@affine_ANTS')
            mbFlow.connect(regANTS, 'nonlinear_transform',
                           Sink, 'transforms.@nonlinear_ANTS')
            mbFlow.connect(regANTS, 'nonlinear_inverse_transform',
                           Sink, 'transforms.@nonlinearinverse_ANTS')
        #---------------------------------------------------------------------
        # Register image volume to template in MNI152 space using FSL's flirt:
        #---------------------------------------------------------------------
        elif vol_reg_method == 'flirt':
            regFlirt = Node(name='FLIRT_standard', interface=FLIRT())
            mbFlow.add_nodes([regFlirt])
            if do_input_nifti:
                mbFlow.connect(niftiBrain, 'nifti',
                               regFlirt, 'in_file')
            else:
                mbFlow.connect(mgh2nifti, 'out_file',
                               regFlirt, 'in_file')
            regFlirt.inputs.bins = 640
            regFlirt.inputs.cost_func = 'mutualinfo'
            regFlirt.inputs.dof = 12
            regFlirt.inputs.reference = volume_template
            regFlirt.inputs.out_matrix_file = 'affine_to_template.mat'
            regFlirt.inputs.out_file = 'affine_to_template.nii.gz'
            mbFlow.connect(regFlirt, 'out_matrix_file',
                           Sink, 'transforms.@affine_flirt')
            mbFlow.connect(regFlirt, 'out_file',
                           Sink, 'transforms.@affine_volume')


#=============================================================================
##############################################################################
#=============================================================================
#
#
#   Surface workflows
#
#
#=============================================================================
##############################################################################
#=============================================================================

if run_SurfFlows:
    mbFlow.connect([(Info, Surf, [('subject','subject'), ('hemi','hemi')])])

    #-------------------------------------------------------------------------
    # Convert surfaces to VTK
    #-------------------------------------------------------------------------
    if not do_input_vtk:
        ConvertSurf = Node(name='Surface_to_vtk',
                           interface=Fn(function = surface_to_vtk,
                                        input_names=['surface_file'],
                                        output_names=['vtk_file']))
        mbFlow.connect(Surf, 'surface_files', ConvertSurf, 'surface_file')
    #-------------------------------------------------------------------------
    # Evaluation inputs: location and structure of atlas surfaces
    #-------------------------------------------------------------------------
    if do_evaluate_surf_labels or init_labels == 'manual':
        Atlas = Node(name='Atlases',
                     interface=DataGrabber(infields=['subject','hemi'],
                                           outfields=['atlas_file'],
                                           sort_filelist=False))
        Atlas.inputs.base_directory = subjects_path
        Atlas.inputs.template = '%s/label/%s.labels.' +\
                                protocol + '.' + atlas_label_type + '.vtk'
        Atlas.inputs.template_args['atlas_file'] = [['subject','hemi']]
    
        mbFlow.connect([(Info, Atlas, [('subject','subject'), ('hemi','hemi')])])

##############################################################################
#
#   Surface label workflow
#
#      - DKT classifier (default option)
#      - FreeSurfer labels
#      - Multi-atlas labeling
#      - Manual labels
#
##############################################################################
if run_SurfLabelFlow and run_SurfFlows:

    SurfLabelFlow = Workflow(name='Surface_labels')

    #=========================================================================
    # Initialize labels with the DKT classifier atlas
    #=========================================================================
    if init_labels == 'DKT_atlas':
        #---------------------------------------------------------------------
        # Label a brain with the DKT atlas using FreeSurfer's mris_ca_label
        #---------------------------------------------------------------------
        Classifier = Node(name='Label_with_DKT_atlas',
                          interface=Fn(function = label_with_classifier,
                                       input_names=['hemi',
                                                    'subject',
                                                    'subjects_path',
                                                    'sphere_file',
                                                    'classifier_path',
                                                    'classifier_atlas'],
                                       output_names=['annot_name',
                                                     'annot_file']))
        SurfLabelFlow.add_nodes([Classifier])
        mbFlow.connect([(Info, SurfLabelFlow,
                           [('hemi', 'Label_with_DKT_atlas.hemi'),
                            ('subject', 'Label_with_DKT_atlas.subject')])])
        Classifier.inputs.subjects_path = subjects_path
        mbFlow.connect(Surf, 'sphere_files',
                         SurfLabelFlow, 'Label_with_DKT_atlas.sphere_file')
        Classifier.inputs.classifier_path = atlas_path
        Classifier.inputs.classifier_atlas = classifier_atlas

        #---------------------------------------------------------------------
        # Convert .annot file to .vtk format
        #---------------------------------------------------------------------
        Classifier2vtk = Node(name='DKT_annot_to_vtk',
                              interface=Fn(function = annot_to_vtk,
                                           input_names=['annot_file',
                                                        'vtk_file'],
                                           output_names=['labels',
                                                         'output_vtk']))
        SurfLabelFlow.add_nodes([Classifier2vtk])
        SurfLabelFlow.connect(Classifier, 'annot_file',
                              Classifier2vtk, 'annot_file')
        if do_input_vtk:
            mbFlow.connect(Surf, 'surface_files',
                           SurfLabelFlow, 'DKT_annot_to_vtk.vtk_file')
        else:
            mbFlow.connect(ConvertSurf, 'vtk_file',
                           SurfLabelFlow, 'DKT_annot_to_vtk.vtk_file')
        mbFlow.connect(SurfLabelFlow, 'DKT_annot_to_vtk.output_vtk',
                       Sink, 'labels.@DKT_surface')
        plug = 'DKT_annot_to_vtk.output_vtk'
        plug1 = Classifier2vtk
        plug2 = 'output_vtk'

    #=========================================================================
    # Initialize labels with FreeSurfer's standard DK classifier atlas
    #=========================================================================
    elif init_labels == 'FreeSurfer':
        FreeLabels = Node(name='DK_annot_to_vtk',
                          interface=Fn(function = annot_to_vtk,
                                       input_names=['annot_file',
                                                    'vtk_file'],
                                       output_names=['labels',
                                                     'output_vtk']))
        SurfLabelFlow.add_nodes([FreeLabels])
        mbFlow.connect(Annot, 'annot_files',
                       SurfLabelFlow, 'DK_annot_to_vtk.annot_file')
        if do_input_vtk:
            mbFlow.connect(Surf, 'surface_files',
                           SurfLabelFlow, 'DK_annot_to_vtk.vtk_file')
        else:
            mbFlow.connect(ConvertSurf, 'vtk_file',
                           SurfLabelFlow, 'DK_annot_to_vtk.vtk_file')
        mbFlow.connect(SurfLabelFlow, 'DK_annot_to_vtk.output_vtk',
                       Sink, 'labels.@Free_surface')
        plug = 'DK_annot_to_vtk.output_vtk'
        plug1 = FreeLabels
        plug2 = 'output_vtk'

        # #=========================================================================
        # # Initialize labels using multi-atlas registration
        # #=========================================================================
        # elif init_labels == 'max_prob':
        #     #---------------------------------------------------------------------
        #     # Register surfaces to average template
        #     #---------------------------------------------------------------------
        #     Register = Node(name='Register_template',
        #                     interface=Fn(function = register_template,
        #                                  input_names=['hemi',
        #                                               'sphere_file',
        #                                               'transform',
        #                                               'atlas_path',
        #                                               'template'],
        #                                  output_names=['transform']))
        #     SurfLabelFlow.add_nodes([Register])
        #     mbFlow.connect(Info, 'hemi', SurfLabelFlow, 'Register_template.hemi')
        #     mbFlow.connect(Surf, 'sphere_files',
        #                    SurfLabelFlow, 'Register_template.sphere_file')
        #     Register.inputs.transform = 'sphere_to_' + free_template + '.reg'
        #     Register.inputs.atlas_path = atlas_path
        #     Register.inputs.template = free_template + '.tif'
        #     #---------------------------------------------------------------------
        #     # Register atlases to subject via template
        #     #---------------------------------------------------------------------
        #     Transform = MapNode(name='Transform_labels',
        #                         iterfield = ['atlas'],
        #                         interface=Fn(function = transform_atlas_labels,
        #                                      input_names=['hemi',
        #                                                   'subject',
        #                                                   'transform',
        #                                                   'subjects_path',
        #                                                   'atlas',
        #                                                   'atlas_string'],
        #                                      output_names=['output_file']))
        #     SurfLabelFlow.add_nodes([Transform])
        #     mbFlow.connect([(Info, SurfLabelFlow,
        #                        [('hemi', 'Transform_labels.hemi'),
        #                         ('subject', 'Transform_labels.subject')])])
        #     SurfLabelFlow.connect(Register, 'transform', Transform, 'transform')
        #     #Transform.inputs.transform = 'sphere_to_' + template + '_template.reg'
        #     Transform.inputs.subjects_path = subjects_path
        #     Transform.inputs.atlas = atlas_list
        #     Transform.inputs.atlas_string = 'labels.'+protocol+'.'+atlas_label_type
        #     #---------------------------------------------------------------------
        #     # Majority vote label
        #     #---------------------------------------------------------------------
        #     Vote = Node(name='Label_vote',
        #                 interface=Fn(function = majority_vote_label,
        #                              input_names=['surface_file',
        #                                           'annot_files'],
        #                              output_names=['labels_max',
        #                                            'label_counts',
        #                                            'label_votes',
        #                                            'consensus_vertices',
        #                                            'maxlabel_file',
        #                                            'labelcounts_file',
        #                                            'labelvotes_file']))
        #     SurfLabelFlow.add_nodes([Vote])
        #     if do_input_vtk:
        #         mbFlow.connect(Surf, 'surface_files',
        #                          SurfLabelFlow, 'Label_vote.surface_file')
        #     else:
        #         mbFlow.connect(ConvertSurf, 'vtk_file',
        #                          SurfLabelFlow, 'Label_vote.surface_file')
        #     SurfLabelFlow.connect(Transform, 'output_file', Vote, 'annot_files')
        #     mbFlow.connect([(SurfLabelFlow, Sink,
        #                        [('Label_vote.maxlabel_file', 'labels.@max'),
        #                         ('Label_vote.labelcounts_file', 'labels.@counts'),
        #                         ('Label_vote.labelvotes_file', 'labels.@votes')])])
        #     plug = 'Label_vote.maxlabel_file'
        #     plug1 = Vote
        #     plug2 = 'maxlabel_file'

    #=========================================================================
    # Skip label initialization and process manual (atlas) labels
    #=========================================================================
    elif init_labels == 'manual':
        AtlasLabels = Node(name='Atlas_labels',
                           interface=Fn(function = read_vtk,
                                        input_names=['input_vtk',
                                                     'return_first',
                                                     'return_array'],
                                        output_names=['faces',
                                                      'lines',
                                                      'indices',
                                                      'points',
                                                      'npoints',
                                                      'scalars',
                                                      'scalar_names',
                                                      'input_vtk']))
        SurfLabelFlow.add_nodes([AtlasLabels])
        mbFlow.connect(Atlas, 'atlas_file',
                       SurfLabelFlow, 'Atlas_labels.input_vtk')
        AtlasLabels.inputs.return_first = 'True'
        AtlasLabels.inputs.return_array = 'False'
        plug = 'Atlas_labels.input_vtk'
        plug1 = AtlasLabels
        plug2 = 'input_vtk'

    #=========================================================================
    # Surface label evaluation against manual labels
    #=========================================================================
    if do_evaluate_surf_labels:

        EvalSurfLabels = Node(name='Evaluate_surface_labels',
                              interface=Fn(function = measure_surface_overlap,
                                           input_names=['command',
                                                        'labels_file1',
                                                        'labels_file2'],
                                           output_names=['overlap_file']))
        mbFlow.add_nodes([EvalSurfLabels])
        surface_overlap_command = os.path.join(ccode_path,
            'surface_overlap', 'SurfaceOverlapMain')
        EvalSurfLabels.inputs.command = surface_overlap_command
        mbFlow.connect(Atlas, 'atlas_file', EvalSurfLabels, 'labels_file1')
        mbFlow.connect(SurfLabelFlow, plug,
                       'EvalSurfLabels.labels_file2')
        mbFlow.connect(EvalSurfLabels, 'overlap_file', Sink, 'evaluate_labels')

    #=========================================================================
    # Convert surface label numbers to volume label numbers
    #=========================================================================
    RelabelSurface = Node(name='Relabel_surface',
                          interface=Fn(function = relabel_surface,
                                       input_names=['vtk_file',
                                                    'hemi',
                                                    'old_labels',
                                                    'new_labels',
                                                    'output_file'],
                                       output_names=['output_file']))
    SurfLabelFlow.add_nodes([RelabelSurface])
    SurfLabelFlow.connect(plug1, plug2, RelabelSurface, 'vtk_file')
    mbFlow.connect(Info, 'hemi', SurfLabelFlow, 'Relabel_surface.hemi')
    RelabelSurface.inputs.old_labels = ''
    RelabelSurface.inputs.new_labels = ''
    RelabelSurface.inputs.output_file = ''
    mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                   Sink, 'labels.@surface')

##############################################################################
#
#   Surface shape measurement workflow
#
#      - Surface area
#      - Travel depth
#      - Geodesic depth
#      - Mean curvature
#      - Convexity (from FreeSurfer)
#      - Thickness (from FreeSurfer)
#
##############################################################################
if run_WholeSurfShapeFlow and run_SurfFlows:

    WholeSurfShapeFlow = Workflow(name='Surface_shapes')

    #=========================================================================
    # Measure surface area
    #=========================================================================
    SurfaceArea = Node(name='Surface_area',
                interface=Fn(function = area,
                             input_names=['command',
                                          'surface_file'],
                             output_names=['area_file']))
    area_command = os.path.join(ccode_path, 'area', 'PointAreaMain')
    SurfaceArea.inputs.command = area_command

    #=========================================================================
    # Measure surface travel depth
    #=========================================================================
    TravelDepth = Node(name='Travel_depth',
                       interface=Fn(function = travel_depth,
                                    input_names=['command',
                                                 'surface_file'],
                                    output_names=['depth_file']))
    WholeSurfShapeFlow.add_nodes([TravelDepth])
    TravelDepth.inputs.command = os.path.join(ccode_path,
                                              'travel_depth',
                                              'TravelDepthMain')

    #=========================================================================
    # Rescale surface travel depth
    #=========================================================================
    if do_fundi:
        RescaleTravelDepth = Node(name='Rescale_travel_depth',
                            interface=Fn(function = rescale_by_neighborhood,
                                         input_names=['input_vtk',
                                                      'indices',
                                                      'nedges',
                                                      'p',
                                                      'set_max_to_1',
                                                      'save_file',
                                                      'output_filestring'],
                                         output_names=['rescaled_scalars',
                                                       'rescaled_scalars_file']))
        WholeSurfShapeFlow.add_nodes([RescaleTravelDepth])
        WholeSurfShapeFlow.connect(TravelDepth, 'depth_file',
                                   RescaleTravelDepth, 'input_vtk')
        RescaleTravelDepth.inputs.indices = []
        RescaleTravelDepth.inputs.nedges = 10
        RescaleTravelDepth.inputs.p = 99
        RescaleTravelDepth.inputs.set_max_to_1 = True
        RescaleTravelDepth.inputs.save_file = True
        RescaleTravelDepth.inputs.output_filestring = 'travel_depth_rescaled'
        #mbFlow.connect(WholeSurfShapeFlow, 'Rescale_travel_depth.rescaled_scalars_file',
        #                 Sink, 'shapes.@travel_depth_rescaled')

    #=========================================================================
    # Measure surface geodesic depth
    #=========================================================================
    GeodesicDepth = Node(name='Geodesic_depth',
                         interface=Fn(function = geodesic_depth,
                                      input_names=['command',
                                                   'surface_file'],
                                      output_names=['depth_file']))
    GeodesicDepth.inputs.command = os.path.join(ccode_path,
                                                'geodesic_depth',
                                                'GeodesicDepthMain')

    #=========================================================================
    # Measure surface curvature
    #=========================================================================
    CurvNode = Node(name='Curvature',
                     interface=Fn(function = curvature,
                                  input_names=['command',
                                               'method',
                                               'arguments',
                                               'surface_file'],
                                  output_names=['mean_curvature_file',
                                                'gauss_curvature_file',
                                                'max_curvature_file',
                                                'min_curvature_file',
                                                'min_curvature_vector_file']))
    CurvNode.inputs.command = os.path.join(ccode_path,
                                           'curvature',
                                           'CurvatureMain')
    CurvNode.inputs.method = 2
    CurvNode.inputs.arguments = '-n 0.7'

    #=========================================================================
    # Convert FreeSurfer surface measures to VTK
    #=========================================================================
    if do_convexity:
        ConvexNode = Node(name='Convexity_to_vtk',
                          interface=Fn(function = curvature_to_vtk,
                                       input_names=['surface_file',
                                                    'vtk_file'],
                                       output_names=['output_vtk']))
        WholeSurfShapeFlow.add_nodes([ConvexNode])
        mbFlow.connect(Surf, 'convexity_files',
                         WholeSurfShapeFlow, 'Convexity_to_vtk.surface_file')
        mbFlow.connect(ConvertSurf, 'vtk_file',
                         WholeSurfShapeFlow, 'Convexity_to_vtk.vtk_file')
        mbFlow.connect(WholeSurfShapeFlow, 'Convexity_to_vtk.output_vtk',
                         Sink, 'shapes.@convexity')
    if do_thickness:
        ThickNode = Node(name='Thickness_to_vtk',
                         interface=Fn(function = curvature_to_vtk,
                                      input_names=['surface_file',
                                                   'vtk_file'],
                                      output_names=['output_vtk']))
        WholeSurfShapeFlow.add_nodes([ThickNode])
        mbFlow.connect(Surf, 'thickness_files',
                         WholeSurfShapeFlow, 'Thickness_to_vtk.surface_file')
        mbFlow.connect(ConvertSurf, 'vtk_file',
                         WholeSurfShapeFlow, 'Thickness_to_vtk.vtk_file')
        mbFlow.connect(WholeSurfShapeFlow, 'Thickness_to_vtk.output_vtk',
                         Sink, 'shapes.@thickness')
    #-------------------------------------------------------------------------
    # Add and connect nodes, save output files
    #-------------------------------------------------------------------------
    WholeSurfShapeFlow.add_nodes([SurfaceArea, GeodesicDepth, CurvNode])
    if do_input_vtk:
        mbFlow.connect([(Surf, WholeSurfShapeFlow,
                           [('surface_files','Surface_area.surface_file'),
                            ('surface_files','Travel_depth.surface_file'),
                            ('surface_files','Geodesic_depth.surface_file'),
                            ('surface_files','Curvature.surface_file')])])
    else:
        mbFlow.connect([(ConvertSurf, WholeSurfShapeFlow,
                           [('vtk_file', 'Surface_area.surface_file'),
                            ('vtk_file', 'Travel_depth.surface_file'),
                            ('vtk_file', 'Geodesic_depth.surface_file'),
                            ('vtk_file', 'Curvature.surface_file')])])
    mbFlow.connect([(WholeSurfShapeFlow, Sink,
                       [('Surface_area.area_file', 'shapes.@surface_area'),
                        ('Travel_depth.depth_file', 'shapes.@travel_depth'),
                        ('Geodesic_depth.depth_file', 'shapes.@geodesic_depth'),
                        ('Curvature.mean_curvature_file', 'shapes.@mean_curvature')])])

##############################################################################
#
#   Surface feature extraction workflow
#
#      - Folds
#      - Sulci
#      - Fundi (curves at the bottoms of folds/sulci)
#
##############################################################################
if run_SurfFeatureFlow and run_SurfFlows:

    SurfFeatureFlow = Workflow(name='Surface_features')
    min_fold_size = 50  # Minimum number or vertices constituting a fold

    #=========================================================================
    # Folds
    #=========================================================================
    FoldsNode = Node(name='Folds',
                     interface=Fn(function = extract_folds,
                                  input_names=['depth_file',
                                               'min_fold_size',
                                               'tiny_depth',
                                               'save_file'],
                                  output_names=['folds',
                                                'n_folds',
                                                'depth_threshold',
                                                'bins',
                                                'bin_edges',
                                                'folds_file']))
    SurfFeatureFlow.add_nodes([FoldsNode])
    mbFlow.connect(WholeSurfShapeFlow, 'Travel_depth.depth_file',
                     SurfFeatureFlow, 'Folds.depth_file')
    FoldsNode.inputs.min_fold_size = min_fold_size
    FoldsNode.inputs.tiny_depth = 0.001
    FoldsNode.inputs.save_file = True
    mbFlow.connect(SurfFeatureFlow, 'Folds.folds_file', Sink, 'features.@folds')

    # # Subfolds
    # SubfoldsNode = Node(name='Subfolds',
    #                     interface=Fn(function = extract_subfolds,
    #                                  input_names=['depth_file',
    #                                               'folds',
    #                                               'depth_factor',
    #                                               'depth_ratio',
    #                                               'tolerance',
    #                                               'save_file'],
    #                                  output_names=['subfolds',
    #                                                'n_subfolds',
    #                                                'subfolds_file']))
    # SurfFeatureFlow.add_nodes([SubfoldsNode])
    # mbFlow.connect(WholeSurfShapeFlow, 'Travel_depth.depth_file',
    #               SurfFeatureFlow, 'Subfolds.depth_file')
    # SurfFeatureFlow.connect(FoldsNode, 'folds', SubfoldsNode, 'folds')
    # SubfoldsNode.inputs.depth_factor = 0.25
    # SubfoldsNode.inputs.depth_ratio = 0.1
    # SubfoldsNode.inputs.tolerance = 0.01
    # SubfoldsNode.inputs.save_file = True
    # # Save subfolds
    # mbFlow.connect(SurfFeatureFlow, 'Subfolds.subfolds_file',
    #               Sink, 'features.@subfolds')

    #=========================================================================
    # Sulci
    #=========================================================================
    if do_sulci:
        SulciNode = Node(name='Sulci',
                         interface=Fn(function = extract_sulci,
                                      input_names=['labels_file',
                                                   'folds_or_file',
                                                   'hemi',
                                                   'sulcus_label_pair_lists',
                                                   'unique_sulcus_label_pairs',
                                                   'min_boundary',
                                                   'sulcus_names'],
                                      output_names=['sulci',
                                                    'n_sulci',
                                                    'sulci_file']))
        SurfFeatureFlow.add_nodes([SulciNode])
        mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                       SurfFeatureFlow, 'Sulci.labels_file')
        SurfFeatureFlow.connect(FoldsNode, 'folds', SulciNode, 'folds_or_file')
        mbFlow.connect(Info, 'hemi', SurfFeatureFlow, 'Sulci.hemi')
        SulciNode.inputs.sulcus_label_pair_lists = sulcus_label_pair_lists
        SulciNode.inputs.unique_sulcus_label_pairs = unique_sulcus_label_pairs
        SulciNode.inputs.min_boundary = 1
        SulciNode.inputs.sulcus_names = sulcus_names
        mbFlow.connect(SurfFeatureFlow, 'Sulci.sulci_file',
                       Sink, 'features.@sulci')

    #=========================================================================
    # Fundi
    #=========================================================================
    if do_fundi:
        LikelihoodNode = Node(name='Likelihood',
                              interface=Fn(function = compute_likelihood,
                                           input_names=['trained_file',
                                                        'depth_file',
                                                        'curvature_file',
                                                        'folds'
                                                        'save_file'],
                                           output_names=['likelihoods',
                                                         'likelihoods_file']))

        SurfFeatureFlow.add_nodes([LikelihoodNode])
        LikelihoodNode.inputs.trained_file = os.path.join(atlas_path,
            'depth_curv_border_nonborder_parameters.pkl')
        mbFlow.connect([(WholeSurfShapeFlow, SurfFeatureFlow,
                           [('Rescale_travel_depth.rescaled_scalars_file',
                             'Likelihood.depth_file'),
                            ('Curvature.mean_curvature_file',
                             'Likelihood.curvature_file')])])
        SurfFeatureFlow.connect(FoldsNode, 'folds', LikelihoodNode, 'folds')
        LikelihoodNode.inputs.save_file = True
        mbFlow.connect(SurfFeatureFlow, 'Likelihood.likelihoods_file',
                         Sink, 'features.@likelihoods')

        FundiNode = Node(name='Fundi',
                         interface=Fn(function = extract_fundi,
                                      input_names=['folds',
                                                   'sulci',
                                                   'likelihoods',
                                                   'rescaled_depth_file',
                                                   'depth_file',
                                                   'min_edges',
                                                   'erosion_ratio',
                                                   'normalize_likelihoods',
                                                   'smooth_skeleton',
                                                   'save_file'],
                                      output_names=['fundi',
                                                    'n_fundi',
                                                    'fundi_file']))
        SurfFeatureFlow.connect(FoldsNode, 'folds', FundiNode, 'folds')
        SurfFeatureFlow.connect(SulciNode, 'sulci', FundiNode, 'sulci')
        SurfFeatureFlow.connect(LikelihoodNode, 'likelihoods', 
                            FundiNode, 'likelihoods')
        mbFlow.connect([(WholeSurfShapeFlow, SurfFeatureFlow,
                       [('Rescale_travel_depth.rescaled_scalars_file',
                         'Fundi.rescaled_depth_file'),
                        ('Travel_depth.depth_file','Fundi.depth_file')])])
        FundiNode.inputs.min_edges = 10
        FundiNode.inputs.erosion_ratio = 0.25
        FundiNode.inputs.normalize_likelihoods = True
        FundiNode.inputs.smooth_skeleton = False
        FundiNode.inputs.save_file = True
        mbFlow.connect(SurfFeatureFlow, 'Fundi.fundi_file',
                         Sink, 'features.@fundi')

##############################################################################
#
#   Surface feature shape workflow
#
#       - Laplace-Beltrami spectra
#
##############################################################################
if run_SurfShapeFlow and run_SurfFlows:
    SurfFeatureShapeFlow = Workflow(name='Surface_feature_shapes')

    if do_measure_spectra:
        #=====================================================================
        # Measure Laplace-Beltrami spectra of labeled regions
        #=====================================================================
        SpectraLabels = Node(name='Spectra_labels',
                             interface=Fn(function = fem_laplacian_from_labels,
                                          input_names=['vtk_file',
                                                       'n_eigenvalues',
                                                       'normalization'],
                                          output_names=['spectrum_lists',
                                                        'label_list']))
        SurfFeatureShapeFlow.add_nodes([SpectraLabels])
        mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                       SurfFeatureShapeFlow, 'Spectra_labels.vtk_file')
        SpectraLabels.inputs.n_eigenvalues = 6
        SpectraLabels.inputs.normalization = "area"

        #=====================================================================
        # Measure Laplace-Beltrami spectra of sulci
        #=====================================================================
        if do_sulci:
            SpectraSulci = SpectraLabels.clone('Spectra_sulci')
            SurfFeatureShapeFlow.add_nodes([SpectraSulci])
            mbFlow.connect(SulciNode, 'sulci_file', SpectraSulci, 'vtk_file')


##############################################################################
#
#   Surface feature shape table workflow
#
#       - Surface feature shape tables:
#           - Label shapes
#           - Sulcus shapes
#           - Fundus shapes
#       - Vertex measures table
#
##############################################################################
if run_SurfFlows:

    #=========================================================================
    # Surface feature shape tables: labels, sulci, fundi
    #=========================================================================
    if do_surf_table:
        ShapeTables = Node(name='Shape_tables',
                           interface=Fn(function = write_shape_stats,
                                        input_names=['labels_or_file',
                                                     'sulci',
                                                     'fundi',
                                                     'affine_transform_file',
                                                     'transform_format',
                                                     'area_file',
                                                     'mean_curvature_file',
                                                     'travel_depth_file',
                                                     'geodesic_depth_file',
                                                     'convexity_file',
                                                     'thickness_file',
                                                     'labels_spectra',
                                                     'labels_spectra_IDs',
                                                     'sulci_spectra',
                                                     'sulci_spectra_IDs',
                                                     'exclude_labels',
                                                     'delimiter'],
                                        output_names=['label_table',
                                                      'sulcus_table',
                                                      'fundus_table']))
        mbFlow.add_nodes([ShapeTables])
        mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                       ShapeTables, 'labels_or_file')
        if do_sulci:
            mbFlow.connect(SurfFeatureFlow, 'Sulci.sulci',
                           ShapeTables, 'sulci')
        else:
            ShapeTables.inputs.sulci = []
        if do_fundi:
            mbFlow.connect(SurfFeatureFlow, 'Fundi.fundi',
                           ShapeTables, 'fundi')
        else:
            ShapeTables.inputs.fundi = []
        if run_VolFlows and run_RegFlows and do_register_standard:
            #if vol_reg_method == 'antsRegister':
            #    mbFlow.connect(regAnts, 'composite_transform',
            #                   ShapeTables, 'affine_transform_file')
            #    # Apply the affine part of a complex transform:
            #    #pickfirst = lambda x: x[:1]
            #    #def pickfirst(x):
            #    #  return lambda x: x[:1]
            #    #mbFlow.connect(regAnts, ('forward_transforms', pickfirst),
            #    #                 ShapeTables, 'affine_transform_file')
            #    ShapeTables.inputs.transform_format = 'mat'
            if vol_reg_method == 'ANTS':
                mbFlow.connect(regANTS, 'affine_transform',
                               ShapeTables, 'affine_transform_file')
                ShapeTables.inputs.transform_format = 'itk'
            elif vol_reg_method == 'flirt':
                mbFlow.connect(regFlirt, 'out_matrix_file',
                               ShapeTables, 'affine_transform_file')
                ShapeTables.inputs.transform_format = 'txt'
        else:
            ShapeTables.inputs.affine_transform_file = None
            ShapeTables.inputs.transform_format = None

        #---------------------------------------------------------------------
        mbFlow.connect([(WholeSurfShapeFlow, ShapeTables,
                           [('Surface_area.area_file',
                             'area_file'),
                            ('Curvature.mean_curvature_file',
                             'mean_curvature_file'),
                            ('Travel_depth.depth_file',
                             'travel_depth_file'),
                            ('Geodesic_depth.depth_file',
                             'geodesic_depth_file')])])
        if do_convexity:
            mbFlow.connect(WholeSurfShapeFlow, 'Convexity_to_vtk.output_vtk',
                           ShapeTables, 'convexity_file')
        else:
            ShapeTables.inputs.convexity_file = ''
        if do_thickness:
            mbFlow.connect(WholeSurfShapeFlow, 'Thickness_to_vtk.output_vtk',
                           ShapeTables, 'thickness_file')
        else:
            ShapeTables.inputs.thickness_file = ''
        if do_measure_spectra:
            mbFlow.connect(SurfFeatureShapeFlow, 'Spectra_labels.spectrum_lists',
                           ShapeTables, 'labels_spectra')
            mbFlow.connect(SurfFeatureShapeFlow, 'Spectra_labels.label_list',
                           ShapeTables, 'labels_spectra_IDs')
            if do_sulci:
                mbFlow.connect(SurfFeatureShapeFlow, 'Spectra_sulci.spectrum_lists',
                               ShapeTables, 'sulci_spectra')
                mbFlow.connect(SurfFeatureShapeFlow, 'Spectra_sulci.label_list',
                               ShapeTables, 'sulci_spectra_IDs')
            else:
                ShapeTables.inputs.sulci_spectra = []
                ShapeTables.inputs.sulci_spectra_IDs = []
        else:
            ShapeTables.inputs.labels_spectra = []
            ShapeTables.inputs.sulci_spectra = []
            ShapeTables.inputs.labels_spectra_IDs = []
            ShapeTables.inputs.sulci_spectra_IDs = []

        #---------------------------------------------------------------------
        ShapeTables.inputs.exclude_labels = [-1]
        ShapeTables.inputs.delimiter = ","
        mbFlow.connect(ShapeTables, 'label_table', Sink, 'tables.@labels')
        if do_sulci:
            mbFlow.connect(ShapeTables, 'sulcus_table', Sink, 'tables.@sulci')
        if do_fundi:
            mbFlow.connect(ShapeTables, 'fundus_table', Sink, 'tables.@fundi')

    #=========================================================================
    # Vertex measures table
    #=========================================================================
    if do_vertex_table:

        VertexTable = Node(name='Vertex_table',
                           interface=Fn(function = write_vertex_measures,
                                        input_names=['table_file',
                                                     'labels_or_file',
                                                     'sulci',
                                                     'fundi',
                                                     'affine_transform_file',
                                                     'transform_format',
                                                     'area_file',
                                                     'mean_curvature_file',
                                                     'travel_depth_file',
                                                     'geodesic_depth_file',
                                                     'convexity_file',
                                                     'thickness_file',
                                                     'delimiter'],
                                        output_names=['shapes_table']))
        mbFlow.add_nodes([VertexTable])
        VertexTable.inputs.table_file = 'vertex_shapes.csv'
        mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                       VertexTable, 'labels_or_file')
        if do_sulci:
            mbFlow.connect(SurfFeatureFlow, 'Sulci.sulci',
                           VertexTable, 'sulci')
        else:
            ShapeTables.inputs.sulci = []
        if do_fundi:
            mbFlow.connect(SurfFeatureFlow, 'Fundi.fundi',
                           VertexTable, 'fundi')
        else:
            ShapeTables.inputs.fundi = []

        if run_VolFlows and run_RegFlows and do_register_standard:
            # if vol_reg_method == 'antsRegister':
            #     mbFlow.connect(regAnts, 'composite_transform',
            #                    VertexTable, 'affine_transform_file')
            #     # Apply the affine part of a complex transform:
            #     #pickfirst = lambda x: x[:1]
            #     #mbFlow.connect(regAnts, ('forward_transforms', pickfirst),
            #     #                 VertexTable, 'affine_transform_file')
            #     VertexTable.inputs.transform_format = 'mat'
            if vol_reg_method == 'ANTS':
                mbFlow.connect(regANTS, 'affine_transform',
                               VertexTable, 'affine_transform_file')
            elif vol_reg_method == 'flirt':
                mbFlow.connect(regFlirt, 'out_matrix_file',
                               VertexTable, 'affine_transform_file')
                VertexTable.inputs.transform_format = 'txt'
        else:
            VertexTable.inputs.affine_transform_file = None
            VertexTable.inputs.transform_format = None
        #---------------------------------------------------------------------
        mbFlow.connect([(WholeSurfShapeFlow, VertexTable,
                           [('Surface_area.area_file','area_file'),
                            ('Travel_depth.depth_file',
                             'travel_depth_file'),
                            ('Geodesic_depth.depth_file',
                             'geodesic_depth_file'),
                            ('Curvature.mean_curvature_file',
                             'mean_curvature_file')])])
        if do_thickness:
            mbFlow.connect(WholeSurfShapeFlow, 'Thickness_to_vtk.output_vtk',
                           VertexTable, 'thickness_file')
        else:
            VertexTable.inputs.thickness_file = ''
        if do_convexity:
            mbFlow.connect(WholeSurfShapeFlow, 'Convexity_to_vtk.output_vtk',
                           VertexTable, 'convexity_file')
        else:
            VertexTable.inputs.convexity_file = ''
        #---------------------------------------------------------------------
        VertexTable.inputs.delimiter = ","
        mbFlow.connect(VertexTable, 'shapes_table',
                         Sink, 'tables.@vertex_table')

    # #---------------------------------------------------------------------
    # # Apply RegFlows's affine transform to surface coordinates:
    # #---------------------------------------------------------------------
    # TransformPoints = Node(name='Transform_points',
    #                        interface=Fn(function = apply_affine_transform,
    #                                     input_names=['transform_file',
    #                                                  'vtk_or_points',
    #                                                  'transform_format',
    #                                                  'save_file'],
    #                                     output_names=['affine_points',
    #                                                   'output_file']))
    # VolLabelFlow.add_nodes([TransformPoints])
    # if vol_reg_method == 'antsRegister':
    #     TransformPoints.inputs.transform_format = 'mat'
    #     VolLabelFlow.connect(regAnts, 'output_transform_prefix',
    #                          TransformPoints, 'transform_file')
    # elif vol_reg_method == 'ANTS':
    #     TransformPoints.inputs.transform_format = 'itk'
    #     VolLabelFlow.connect(regANTS, 'affine_transform',
    #                          TransformPoints, 'transform_file')
    # elif vol_reg_method == 'flirt':
    #     TransformPoints.inputs.transform_format = 'txt'
    #     VolLabelFlow.connect(regFlirt, 'out_matrix_file',
    #                          TransformPoints, 'transform_file')
    # SurfShapeFlow.connect(TravelDepth, 'depth_file',
    #                       TransformPoints, 'vtk_or_points')
    # TransformPoints.inputs.save_file = True
    # mbFlow.connect(SurfShapeFlow, 'Transform_points.output_file',
    #                Sink, 'transforms.@points_to_template')


#=============================================================================
##############################################################################
#=============================================================================
#
#
#   Volume workflows
#
#
#=============================================================================
##############################################################################
#=============================================================================

##############################################################################
#
#   Volume label workflow
#
#       - Fill cortical gray matter with labels
#       - Label subcortical volumes
#       - Evaluate volume labels
#
##############################################################################
if run_VolLabelFlow and run_VolFlows:
    VolLabelFlow = Workflow(name='Volume_labels')

    #=========================================================================
    # Fill cortical gray matter volume with surface vertex labels
    #=========================================================================
    if do_fill_cortex:
        use_ants_to_fill = True
        if use_ants_to_fill:

            FillCortex = Node(name='Fill_cortex',
                              interface=Fn(function = fill_volume_with_surface_labels,
                                           input_names=['volume_mask',
                                                        'surface_files',
                                                        'output_file',
                                                        'binarize'],
                                           output_names=['output_file']))
            VolLabelFlow.add_nodes([FillCortex])
            if do_input_mask:
                mbFlow.connect([(Info, niftiMask,
                                 [('hemi', 'hemi'), ('subject', 'subject')])])
                mbFlow.connect(niftiMask, 'nifti_mask',
                               VolLabelFlow, 'Fill_cortex.volume_mask')
            else:
                mbFlow.connect([(Info, mghMask,
                                 [('hemi', 'hemi'), ('subject', 'subject')])])

                mbFlow.connect(mgh_mask2nifti, 'out_file',
                               VolLabelFlow, 'Fill_cortex.volume_mask')
            if run_SurfFlows and run_SurfLabelFlow:
                mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                               VolLabelFlow, 'Fill_cortex.surface_files')
            else:
                sys.exit('No input surface file specified for Fill_cortex.')
            FillCortex.inputs.output_file = os.path.join(os.getcwd(),
                                                         'filled_cortex.nii.gz')
            FillCortex.inputs.binarize = False
            mbFlow.connect(VolLabelFlow, 'Fill_cortex.output_file',
                           Sink, 'labels.@cortex_volume')

    #=========================================================================
    # Label non-cortical volumes
    #=========================================================================
    if run_RegFlows and do_register_standard and do_label_whole_volume:

        # Inverse transform subcortical label volumes to subject via template
        LabelVolume = Node(name='Label_volume',
                           interface=Fn(function = WarpImageMultiTransform,
                                        input_names=['source',
                                                     'target',
                                                     'output'
                                                     'interp',
                                                     'xfm_stem',
                                                     'affine_transform',
                                                     'nonlinear_transform',
                                                     'inverse',
                                                     'affine_only'],
                                        output_names=['output']))
        VolLabelFlow.add_nodes([LabelVolume])
        LabelVolume.inputs.source = os.path.join(atlas_path, atlas_volume)
        if do_input_nifti:
            mbFlow.connect(niftiBrain, 'nifti',
                           VolLabelFlow, 'Label_volume.target')
        else:
            mbFlow.connect(mgh2nifti, 'out_file',
                           VolLabelFlow, 'Label_volume.target')
        LabelVolume.inputs.output = os.path.join(os.getcwd(),
                                                 'label_volume.nii.gz')
        LabelVolume.inputs.interp = '--use-NN'
        if vol_reg_method == 'ANTS':
            mbFlow.connect(regANTS, 'output_stem',
                           VolLabelFlow, 'Label_volume.xfm_stem')
            LabelVolume.inputs.affine_transform = ''
            LabelVolume.inputs.nonlinear_transform = ''
        else:
            sys.exit('No other vol_reg_method set up.')
        LabelVolume.inputs.inverse = True
        LabelVolume.inputs.affine_only = False
        mbFlow.connect(VolLabelFlow, 'Label_volume.output',
                       Sink, 'labels.@noncortex_volume')

        #=====================================================================
        # Combine cortical and noncortical volume labels
        #=====================================================================
        if do_fill_cortex:
            CombineLabels = Node(name='Combine_labels',
                                 interface=Fn(function = overwrite_volume_labels,
                                              input_names=['source',
                                                           'target',
                                                           'output_file',
                                                           'ignore_labels'],
                                              output_names=['output_file']))
            VolLabelFlow.add_nodes([CombineLabels])
            VolLabelFlow.connect(FillCortex, 'output_file',
                                 CombineLabels, 'source')
            VolLabelFlow.connect(LabelVolume, 'output',
                                 CombineLabels, 'target')
            CombineLabels.inputs.output_file = os.path.join(os.getcwd(),
                                                            'combined_labels.nii.gz')
            CombineLabels.inputs.ignore_labels = [0]
            mbFlow.connect(VolLabelFlow, 'Combine_labels.output_file',
                           Sink, 'labels.@brain')

    #=========================================================================
    # Evaluate label volume overlaps
    #=========================================================================
    if do_evaluate_vol_labels:

        #---------------------------------------------------------------------
        # Evaluation inputs: location and structure of atlas volumes
        #---------------------------------------------------------------------
        AtlasVol = Node(name='Atlas_volume',
                        interface=DataGrabber(infields=['subject'],
                                              outfields=['atlas_vol_file'],
                                              sort_filelist=False))
        VolLabelFlow.add_nodes([AtlasVol])
        AtlasVol.inputs.base_directory = subjects_path
        AtlasVol.inputs.template = '%s/mri/labels.' + protocol + '.manual.nii.gz'
        AtlasVol.inputs.template_args['atlas_vol_file'] = [['subject']]
        mbFlow.connect(Info, 'subject', VolLabelFlow, 'Atlas_volume.subject')
        #---------------------------------------------------------------------
        # Evaluate volume labels
        #---------------------------------------------------------------------
        EvalVolLabels = Node(name='Evaluate_volume_labels',
                             interface=Fn(function = measure_volume_overlap,
                                          input_names=['labels',
                                                       'file2',
                                                       'file1'],
                                          output_names=['overlaps',
                                                        'out_file']))
        VolLabelFlow.add_nodes([EvalVolLabels])
        if do_label_whole_volume and do_fill_cortex:
            EvalVolLabels.inputs.labels = label_numbers
        elif do_fill_cortex:
            EvalVolLabels.inputs.labels = cortex_numbers
        VolLabelFlow.connect(AtlasVol, 'atlas_vol_file', EvalVolLabels, 'file2')
        VolLabelFlow.connect(FillCortex, 'output_file', EvalVolLabels, 'file1')
        mbFlow.connect(VolLabelFlow, 'Evaluate_volume_labels.out_file',
                        Sink, 'evaluate_labels_volume')


##############################################################################
#
#   Volume feature shape workflow
#
#       - Volumes of labeled regions
#
##############################################################################
if run_VolShapeFlow and run_VolLabelFlow and run_VolFlows:
    VolShapeFlow = Workflow(name='Volume_feature_shapes')

    #=========================================================================
    # Measure volume of each region of a labeled image file
    #=========================================================================
    MeasureVolumes = Node(name='Measure_volumes',
                          interface=Fn(function = volume_per_label,
                                       input_names=['labels',
                                                    'input_file'],
                                       output_names=['labels_volumes']))
    VolShapeFlow.add_nodes([MeasureVolumes])
    if do_label_whole_volume and do_fill_cortex:
        MeasureVolumes.inputs.labels = label_numbers
    elif do_fill_cortex:
        MeasureVolumes.inputs.labels = cortex_numbers
    if do_label_whole_volume and do_fill_cortex:
        mbFlow.connect(VolLabelFlow, 'Combine_labels.output_file',
                       VolShapeFlow, 'Measure_volumes.input_file')
    elif do_fill_cortex:
        mbFlow.connect(VolLabelFlow, 'Fill_cortex.output_file',
                       VolShapeFlow, 'Measure_volumes.input_file')
    else:
        sys.exit('No alternative set of label volumes provided...')

    #=========================================================================
    # Create a table to save the volume measures
    #=========================================================================
    LabelVolTable = Node(name='Label_volume_table',
                         interface=Fn(function = write_columns,
                                      input_names=['columns',
                                                   'column_names',
                                                   'output_table',
                                                   'delimiter',
                                                   'quote',
                                                   'input_table'],
                                      output_names=['output_table']))
    mbFlow.add_nodes([LabelVolTable])
    mbFlow.connect(VolShapeFlow, 'Measure_volumes.labels_volumes',
                   LabelVolTable, 'columns')
    LabelVolTable.inputs.column_names = ['label', 'volume']
    LabelVolTable.inputs.output_table = 'label_volumes.csv'
    LabelVolTable.inputs.delimiter = ','
    LabelVolTable.inputs.quote = True
    LabelVolTable.inputs.input_table = ''
    mbFlow.connect(LabelVolTable, 'output_table',
                    Sink, 'tables.@volume_labels')


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

    generate_graphs = True
    if generate_graphs:
        mbFlow.write_graph(graph2use='flat')
        mbFlow.write_graph(graph2use='hierarchical')
    mbFlow.run()

"""
# Script for running Mindboggle on Mindboggle-101 set
import os
from mindboggle.utils.io_table import read_columns

out_path = '/homedir/Data/Mindboggle-101/'
x_path = os.path.join(os.environ['MINDBOGGLE'], 'x')
atlas_list_file = '/home/arno/Data/Brains/Mindboggle101/code/mindboggle101_atlases.txt'
atlas_list = read_columns(atlas_list_file, 1)[0]

for atlas in atlas_list:
    if 'HLN' in atlas or 'Twins' in atlas or \
       'Colin' in atlas or 'After' in atlas or \
       'MMRR-3T7T' in atlas:
    #if 'MMRR-21' in atlas:
    #if 'OASIS-TRT' in atlas:
    #if 'NKI-TRT' in atlas:
    #if 'NKI-RS' in atlas:
        cmd = ' '.join(['python pipeline.py', out_path, atlas])
        print(cmd); os.system(cmd)
"""
