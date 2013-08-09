#!/usr/bin/env python
"""
This is the main program to run Mindboggle.

For help in using Mindboggle ::

    - Online `documentation <http://mindboggle.info/documentation.html>`_
    - README file
    - Help on the command line::

        >>> python mindboggler.py --help

This file uses Nipype (http://www.nipy.org/nipype/) to create a workflow
environment to enable Mindboggle to run in a flexible, modular manner
while storing provenance information.

Authors:
    - Arno Klein, 2010-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Satrajit S. Ghosh, 2013  (satra@mit.edu)  http://www.mit.edu/~satra/
    - Each file lists Mindboggle team members who contributed to its content.

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

#=============================================================================
#
#   Import libraries
#
#=============================================================================
import os
import sys
import argparse
#-----------------------------------------------------------------------------
# Nipype libraries
#-----------------------------------------------------------------------------
from nipype.pipeline.engine import Workflow, Node
from nipype.interfaces.utility import Function as Fn
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import DataGrabber, DataSink
from nipype.interfaces.freesurfer import MRIConvert
#from nipype.interfaces.fsl import FLIRT
#from nipype.interfaces.ants import Registration
#-----------------------------------------------------------------------------
# Mindboggle libraries
#-----------------------------------------------------------------------------
from mindboggle.utils.io_vtk import read_vtk
from mindboggle.utils.io_table import write_columns, \
    write_shape_stats, write_vertex_measures
from mindboggle.data import hashes_url
from mindboggle.utils.io_uri import retrieve_data
from mindboggle.utils.io_free import surface_to_vtk, curvature_to_vtk, \
    annot_to_vtk
from mindboggle.utils.ants import ANTS, WarpImageMultiTransform, ImageMath, \
    fill_volume_with_surface_labels
from mindboggle.labels.protocol import dkt_protocol
from mindboggle.labels.label_free import label_with_classifier
from mindboggle.labels.relabel import relabel_surface, overwrite_volume_labels
from mindboggle.shapes.measure import area, travel_depth, geodesic_depth, \
    curvature, volume_per_label, rescale_by_neighborhood
from mindboggle.utils.mesh import decimate_file
from mindboggle.shapes.laplace_beltrami import spectrum_per_label
from mindboggle.shapes.zernike.zernike import zernike_moments_per_label
from mindboggle.shapes.likelihood import compute_likelihood
from mindboggle.features.folds import extract_folds
from mindboggle.features.fundi import extract_fundi
from mindboggle.utils.paths import smooth_skeleton
from mindboggle.features.sulci import extract_sulci
from mindboggle.evaluate.evaluate_labels import measure_surface_overlap, \
    measure_volume_overlap

#=============================================================================
#
#   Command-line arguments
#
#=============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("subjects",
                    help=("Example: \"python %(prog)s sub1 sub2 sub3\" "
                          "\"sub1\",... are subject names corresponding to "
                          "subject directories within $SUBJECTS_DIR"),
                    nargs='+')
parser.add_argument("-o", help="output directory: \"-o $HOME/mindboggled\" "
                               "(default)",
                    default=os.path.join(os.environ['HOME'], 'mindboggled'),
                    metavar='')
parser.add_argument("-n",
                    help=("number of processors: \"-n 1\" (default)"),
                    type=int,
                    default=1, metavar='')
parser.add_argument("-c", action='store_true',
                    help="Use HTCondor cluster")
# Labels:
parser.add_argument("--only_label", action='store_true',
                    help="only output labels (no features or tables)")
parser.add_argument("--no_labels", action='store_true',
                    help="do not label surfaces or volumes")
parser.add_argument("--surface_labels",
                    help=("surface labels: "
                    "\"atlas\" (default), \"freesurfer\", \"manual\""),
                    choices=['atlas', 'freesurfer', 'manual'],
                    default='atlas', metavar='')
parser.add_argument("--atlases", help=("label with additional volume "
                                       "atlas(es) in MNI152 space"),
                    nargs='+', metavar='')
# Registration:
parser.add_argument("--no_register", help="do not register to volume template (ANTS)",
                    action='store_true')
parser.add_argument("--register_iters", help="ANTS iterations: "
                    "\"--register_iters 33x99x11\" (default)",
                    default='33x99x11', metavar='')
# Features:
parser.add_argument("--no_features", action='store_true',
                    help="do not extract sulci or fundi")
#parser.add_argument("--no_fundi", action='store_true',
#                    help="do not extract fundi")
parser.add_argument("--fundi", action='store_true',
                    help="extract fundi (default is not to)")
# Tables:
parser.add_argument("--no_tables", action='store_true',
                    help="do not generate shape tables")
parser.add_argument("--vertex_table", action='store_true',
                    help=("make table of per-vertex surface shape measures"))
# Shapes:
parser.add_argument("--spectra_values",
                    help="number of eigenvalues: \"--spectra_values 10\" "
                         "(default)",
                    default=10, type=int, metavar='')
#parser.add_argument("--no_zernike", action='store_true',
#                    help="do not compute Zernike moments")
#parser.add_argument("--zernike", action='store_true',
#                    help="compute Zernike moments (default is not to)")
#parser.add_argument("--order", help="order of Zernike moments: "
#                                    "\"--order 10\" (default)",
#                    default=10, type=int, metavar='')
#parser.add_argument("--reduction", help="mesh decimation fraction: "
#                                        "\"--reduction 0.75\" (default)",
#                    default=0.75, type=float, metavar='')
# Top settings:
parser.add_argument("--no_volumes", action='store_true',
                    help="do not process volumes")
parser.add_argument("--no_surfaces", action='store_true',
                    help="do not process surfaces")
parser.add_argument("--no_freesurfer", action='store_true',
                    help="do not use FreeSurfer or its outputs (UNTESTED)")
parser.add_argument("--visual", help=("generate py/graphviz workflow visual"),
                    choices=['hier', 'flat', 'exec'], metavar='')
parser.add_argument("--version", help="version number",
                    action='version', version='%(prog)s 0.1')
args = parser.parse_args()
#-----------------------------------------------------------------------------
# Arguments:
#-----------------------------------------------------------------------------
subjects = args.subjects
output_path = args.o
n_processes = args.n
cluster = args.c

# Labels:
only_label = args.only_label
no_labels = args.no_labels
if no_labels:
    do_label = False
else:
    do_label = True
surface_labels = args.surface_labels
atlases = args.atlases

# Registration:
no_register = args.no_register
register_iters = args.register_iters

# Features:
no_features = args.no_features
#no_fundi = args.no_fundi
do_fundi = args.fundi

# Tables:
no_tables = args.no_tables
vertex_table = args.vertex_table

# Shapes:
no_spectra = False #args.no_spectra
n_eigenvalues = args.spectra_values
#no_zernike = args.no_zernike
do_zernike = False #args.zernike
zernike_order = 10 #args.order
reduction = 0.75 #args.reduction

# General:
no_volumes = args.no_volumes
no_surfaces = args.no_surfaces
no_freesurfer = args.no_freesurfer
graph_vis = args.visual
if graph_vis == 'hier':
    graph_vis = 'hierarchical'
#-----------------------------------------------------------------------------
# Non-FreeSurfer input:
#-----------------------------------------------------------------------------
do_input_vtk = False  # Load VTK surfaces directly (not FreeSurfer surfaces)
do_input_nifti = False  # Load nifti directly (not FreeSurfer mgh file)
do_input_mask = False  # Load nifti directly (not FreeSurfer mgh file)
if no_freesurfer:
    do_input_vtk = True
    do_input_nifti = True
    do_input_mask = True
#-----------------------------------------------------------------------------
# Volume workflows:
#-----------------------------------------------------------------------------
run_VolFlows = False
run_VolShapeFlow = False
run_VolLabelFlow = False
if not no_volumes:
    run_VolFlows = True
    if not only_label and not no_tables:
        run_VolShapeFlow = True
    if do_label:
        run_VolLabelFlow = True
#-----------------------------------------------------------------------------
# Registration to template:
#-----------------------------------------------------------------------------
vol_reg_method = 'ANTS'
save_transforms = True  # NOTE: must save transforms!
if no_register or no_volumes:
    do_register = False
else:
    do_register = True
#-----------------------------------------------------------------------------
# Surface workflows:
#-----------------------------------------------------------------------------
run_SurfFlows = False
run_WholeSurfShapeFlow = False
run_SurfFeatureFlow = False
run_SurfLabelFlow = False
if not no_surfaces:
    run_SurfFlows = True
    if not only_label:
        run_WholeSurfShapeFlow = True
    if do_label:
        run_SurfLabelFlow = True
        if not only_label and not no_features:
            run_SurfFeatureFlow = True
#-----------------------------------------------------------------------------
# Surface features:
#-----------------------------------------------------------------------------
do_folds = False  # Extract folds
do_sulci = False  # Extract sulci
#do_fundi = False  # Extract fundi
do_smooth_fundi = False
if run_SurfFeatureFlow:
    do_folds = True
    do_sulci = True
#        if not no_fundi:
#            do_fundi = True
#            do_smooth_fundi = True
if do_fundi:
    do_smooth_fundi = True
#-----------------------------------------------------------------------------
# Surface shapes:
#-----------------------------------------------------------------------------
do_decimate = False
do_spectra = False  # Measure Laplace-Beltrami spectra for features
#do_zernike = False  # Measure Zernike moments for features
do_thickness = False  # Include FreeSurfer's thickness measure
do_convexity = False  # Include FreeSurfer's convexity measure (sulc.pial)
if run_WholeSurfShapeFlow:
    if not no_freesurfer:
        do_thickness = True
        do_convexity = True
if run_SurfFeatureFlow and not no_tables:
    if not no_spectra:
        do_spectra = True
#    if not no_zernike:
#        do_zernike = True
#        do_decimate = True
    if do_zernike:
        do_decimate = True
else:
    do_zernike = False

#-----------------------------------------------------------------------------
# Labels:
#-----------------------------------------------------------------------------
classifier_name = 'DKTatlas40'
do_label_surf = False
do_register_labels = False  # Label whole brain via volume registration
do_mask_cortex = False  # Mask whole brain labels with gray matter mask
do_fill_cortex = False  # Fill cortical gray matter with surface labels
if do_label:
    if do_register and not no_volumes:
        do_register_labels = True
        do_mask_cortex = True
    if not no_surfaces:
        do_label_surf = True
        if not no_volumes:
            do_fill_cortex = True
#-----------------------------------------------------------------------------
# Tables:
#-----------------------------------------------------------------------------
do_vol_table = False
do_surf_table = False
# Surface/volume feature shape measures:
if not no_tables and not only_label:
    if run_VolLabelFlow:
        do_vol_table = True
    do_surf_table = True
# Per-vertex surface shape measures:
do_vertex_table = False
if vertex_table and not only_label:
    do_vertex_table = True

#=============================================================================
#
#   Hidden arguments: paths, label and template data
#
#=============================================================================
#-----------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------
subjects_path = os.environ['SUBJECTS_DIR']  # FreeSurfer subjects directory
ccode_path = os.environ['MINDBOGGLE_TOOLS']  # Mindboggle C++ code directory
#-----------------------------------------------------------------------------
# Hashes to verify retrieved data
#-----------------------------------------------------------------------------
hashes, url, cache_env, cache = hashes_url()
#-----------------------------------------------------------------------------
# Cache and output directories
#-----------------------------------------------------------------------------
if cache_env in os.environ.keys():
    cache = os.environ[cache_env]
if not os.path.exists(cache):
    print("Create missing cache directory: {0}".format(cache))
    os.mkdir(cache)
temp_path = os.path.join(cache, 'workspace')  # Where to save workflow files
if not os.path.isdir(temp_path):
    os.makedirs(temp_path)
if not os.path.isdir(output_path):
    os.makedirs(output_path)
#-----------------------------------------------------------------------------
# Protocol information
#-----------------------------------------------------------------------------
sulcus_names, sulcus_label_pair_lists, unique_sulcus_label_pairs, \
    label_names, label_numbers, cortex_names, cortex_names_DKT25, \
    cortex_numbers, cortex_numbers_DKT25, noncortex_names, \
    noncortex_numbers = dkt_protocol()
#-----------------------------------------------------------------------------
# Volume template and atlas:
#-----------------------------------------------------------------------------
template_volume = 'OASIS-TRT-20_template_to_MNI152.nii.gz'
atlas_volumes = ['OASIS-TRT-20_atlas_to_MNI152.nii.gz']
if atlases:
    atlas_volumes.extend(atlases)
#-----------------------------------------------------------------------------
# Surface atlas labels:
# - 'manual': manual edits
# - FUTURE: <'adjusted': manual edits after automated alignment to fundi>
#-----------------------------------------------------------------------------
atlas_label_type = 'manual'
#-----------------------------------------------------------------------------
# Evaluation
#-----------------------------------------------------------------------------
do_evaluate_surf_labels = False  # Surface overlap: auto vs. manual labels
do_evaluate_vol_labels = False  # Volume overlap: auto vs. manual labels

#=============================================================================
#
#   Initialize workflow inputs and outputs
#
#=============================================================================
mbFlow = Workflow(name='Mindboggle')
mbFlow.base_dir = temp_path
#-----------------------------------------------------------------------------
# Iterate inputs over subjects, hemispheres, and atlases
# (surfaces are assumed to take the form: lh.pial or lh.pial.vtk)
#-----------------------------------------------------------------------------
if isinstance(atlas_volumes, str):
    atlas_volumes = list(atlas_volumes)
InputVolAtlases = Node(name='Input_volume_atlases',
                    interface=IdentityInterface(fields=['atlas']))
InputVolAtlases.iterables = ('atlas', atlas_volumes)
InputSubjects = Node(name='Input_subjects',
                     interface=IdentityInterface(fields=['subject']))
InputSubjects.iterables = ('subject', subjects)
InputHemis = Node(name='Input_hemispheres',
                  interface=IdentityInterface(fields=['hemi']))
InputHemis.iterables = ('hemi', ['lh','rh'])
#-------------------------------------------------------------------------
# Outputs
#-------------------------------------------------------------------------
Sink = Node(DataSink(), name='Results')
Sink.inputs.base_directory = output_path
Sink.inputs.container = ''
Sink.inputs.substitutions = [('_hemi_lh', 'left'),
    ('_hemi_rh', 'right'),
    ('_subject_', ''),
    ('_atlas_', ''),
    ('smooth_skeletons.vtk', 'smooth_fundi.vtk'),
    ('propagated_labels.nii.gz', 'cortical_surface_labels.nii.gz'),
    ('multiplied_volumes.nii.gz', 'cortical_surface_masked_volume_labels.nii.gz'),
    ('combined_labels.nii.gz',
     'cortical_surface_and_noncortical_volume_labels.nii.gz'),
    ('transformed.nii.gz', 'registered_volume_labels.nii.gz'),
    ('brain_to_OASIS-TRT-20_template_to_MNI152', 'brain_to_template_'),
    ('OASIS-TRT-20to_MNI152.nii.gz', 'template')]

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
    Surf.inputs.template_args['surface_files'] = [['subject', 'hemi', 'pial']]
    Surf.inputs.template_args['sphere_files'] = [['subject', 'hemi','sphere']]
    if do_thickness:
        Surf.inputs.template_args['thickness_files'] = \
            [['subject', 'hemi', 'thickness']]
    if do_convexity:
        Surf.inputs.template_args['convexity_files'] = \
            [['subject', 'hemi', 'sulc']]
    #-------------------------------------------------------------------------
    # Location and structure of the FreeSurfer label inputs
    #-------------------------------------------------------------------------
    if not no_freesurfer and do_label and surface_labels == 'freesurfer':
        Annot = Node(name='Annots',
                     interface=DataGrabber(infields=['subject', 'hemi'],
                                           outfields=['annot_files'],
                                           sort_filelist=False))
        Annot.inputs.base_directory = subjects_path
        Annot.inputs.template = '%s/label/%s.aparc.annot'
        Annot.inputs.template_args['annot_files'] = [['subject', 'hemi']]

if run_VolFlows or do_register:
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
        mbFlow.connect(InputSubjects, 'subject', niftiBrain, 'subject')
    else:
        mghBrain = Node(name='mgh',
                        interface=DataGrabber(infields=['subject'],
                                              outfields=['mgh'],
                                              sort_filelist=False))
        mghBrain.inputs.base_directory = subjects_path
        mghBrain.inputs.template = '%s/mri/brain.mgz'
        mghBrain.inputs.template_args['mgh'] = [['subject']]
        mbFlow.connect(InputSubjects, 'subject', mghBrain, 'subject')

        # FreeSurfer mgh of original for converting from conformal (below):
        mghOrig = Node(name='mgh_orig',
                       interface=DataGrabber(infields=['subject'],
                                             outfields=['mgh_orig'],
                                             sort_filelist=False))
        mghOrig.inputs.base_directory = subjects_path
        mghOrig.inputs.template = '%s/mri/orig/001.mgz'
        mghOrig.inputs.template_args['mgh_orig'] = [['subject']]
        mbFlow.connect(InputSubjects, 'subject', mghOrig, 'subject')

        # Convert FreeSurfer mgh conformal brain volume to nifti format:
        mgh2nifti = Node(name='mgh_to_nifti', interface=MRIConvert())
        mbFlow.add_nodes([mgh2nifti])
        mbFlow.connect(mghBrain, 'mgh', mgh2nifti, 'in_file')
        mbFlow.connect(mghOrig, 'mgh_orig', mgh2nifti, 'reslice_like')
        mgh2nifti.inputs.resample_type = 'interpolate'
        mgh2nifti.inputs.out_type = 'niigz'
        mgh2nifti.inputs.out_file = 'brain.nii.gz'
        #mbFlow.connect(mgh2nifti, 'out_file', Sink, 'brain')

if run_VolLabelFlow and (do_fill_cortex or do_mask_cortex):
    if do_input_mask:
        niftiMask = Node(name='nifti_mask',
                         interface=DataGrabber(infields=['subject', 'hemi'],
                                               outfields=['mask'],
                                               sort_filelist=False))
        niftiMask.inputs.base_directory = subjects_path
        niftiMask.inputs.template = '%s/mri/%s.ribbon.nii.gz'
        niftiMask.inputs.template_args['mask'] = [['subject', 'hemi']]
    else:
        mghMask = Node(name='mgh_mask',
                       interface=DataGrabber(infields=['subject', 'hemi'],
                                             outfields=['mask'],
                                             sort_filelist=False))
        mghMask.inputs.base_directory = subjects_path
        mghMask.inputs.template = '%s/mri/%s.ribbon.mgz'
        mghMask.inputs.template_args['mask'] = [['subject', 'hemi']]
        mbFlow.connect(InputSubjects, 'subject', mghMask, 'subject')
        mbFlow.connect(InputHemis, 'hemi', mghMask, 'hemi')

        # Convert FreeSurfer mgh conformal gray matter mask to nifti format:
        mgh_mask2nifti = Node(name='mgh_mask_to_nifti',
                              interface=MRIConvert())
        mbFlow.add_nodes([mgh_mask2nifti])
        mbFlow.connect(mghMask, 'mask', mgh_mask2nifti, 'in_file')
        mbFlow.connect(mghOrig, 'mgh_orig', mgh_mask2nifti, 'reslice_like')
        mgh_mask2nifti.inputs.resample_type = 'nearest'
        mgh_mask2nifti.inputs.out_type = 'niigz'
        mgh_mask2nifti.inputs.out_file = 'mask.nii.gz'
        #mbFlow.connect(mgh_mask2nifti, 'out_file', Sink, 'brain.@mask')


#=============================================================================
#-----------------------------------------------------------------------------
#
#   Register image volume to template in MNI152 space
#
#-----------------------------------------------------------------------------
#=============================================================================
if do_register:

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
        volume_template_file = retrieve_data(template_volume, url,
                                             hashes, cache_env, cache)
        regANTS.inputs.target = volume_template_file
        if do_register_labels:
            regANTS.inputs.iterations = register_iters
        else:
            regANTS.inputs.iterations = '0'
        regANTS.inputs.output_stem = ''
        if save_transforms:
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
            mbFlow.connect(niftiBrain, 'nifti.@brain',
                           regFlirt, 'in_file')
        else:
            mbFlow.connect(mgh2nifti, 'out_file',
                           regFlirt, 'in_file')
        regFlirt.inputs.bins = 640
        regFlirt.inputs.cost_func = 'mutualinfo'
        regFlirt.inputs.dof = 12
        volume_template_file = retrieve_data(template_volume, url,
                                      hashes, cache_env, cache)
        regFlirt.inputs.reference = volume_template_file
        regFlirt.inputs.out_matrix_file = 'affine_to_template.mat'
        regFlirt.inputs.out_file = 'affine_to_template.nii.gz'
        if save_transforms:
            mbFlow.connect(regFlirt, 'out_matrix_file',
                           Sink, 'transforms.@affine_flirt')
            mbFlow.connect(regFlirt, 'out_file',
                           Sink, 'transforms.@affine_volume')
    # #---------------------------------------------------------------------
    # # Register image volume to template in MNI152 space using antsRegister:
    # #---------------------------------------------------------------------
    # elif vol_reg_method == 'antsRegister':
    #     regAnts = Node(name='antsRegister_standard', interface=Registration())
    #     mbFlow.add_nodes([regAnts])
    #     if do_input_nifti:
    #         mbFlow.connect(niftiBrain, 'nifti',
    #                         regAnts, 'moving_image')
    #     else:
    #         mbFlow.connect(mgh2nifti, 'out_file',
    #                         regAnts, 'moving_image')
    #     regAnts.inputs.fixed_image = [volume_template]
    #     regAnts.inputs.num_threads = 2
    #     regAnts.inputs.winsorize_lower_quantile = 0.01
    #     regAnts.inputs.winsorize_upper_quantile = 0.99
    #     regAnts.inputs.output_warped_image = False
    #     regAnts.inputs.dimension = 3
    #     regAnts.inputs.write_composite_transform = True
    #     regAnts.inputs.collapse_output_transforms = True
    #     regAnts.inputs.write_composite_transform = True
    #     regAnts.inputs.output_transform_prefix = 'standard_'
    #     if do_register_labels:
    #         regAnts.inputs.transforms = ['Rigid', 'Affine', 'SyN']
    #         regAnts.inputs.transform_parameters = [(0.1,), (0.1,), (0.1, 3.0, 0.0)]
    #         regAnts.inputs.number_of_iterations = [[1000,500,250,100]]*2 + [[100,100,70,20]]
    # #                regAnts.inputs.number_of_iterations = [[10,5,2,1]]*2 + [[1,1,7,2]]
    #         regAnts.inputs.metric = ['MI']*2 + ['CC']
    #         regAnts.inputs.metric_weight = [1]*3
    #         regAnts.inputs.radius_or_number_of_bins = [32]*2 + [4]
    #         regAnts.inputs.sampling_strategy = ['Regular']*2 + [None]
    #         regAnts.inputs.sampling_percentage = [0.25]*2 + [None]
    #         regAnts.inputs.convergence_threshold = [1.e-8]*2 + [1e-9]
    #         regAnts.inputs.convergence_window_size = [10]*2 + [15]
    #         regAnts.inputs.smoothing_sigmas = [[3,2,1,0]]*3
    #         regAnts.inputs.sigma_units = ['mm']*3
    #         regAnts.inputs.shrink_factors = [[8,4,2,1]]*2 + [[6,4,2,1]]
    #         regAnts.inputs.use_estimate_learning_rate_once = [True, True, True]
    #         regAnts.inputs.use_histogram_matching = [False]*2 + [True]
    #         regAnts.inputs.initial_moving_transform_com = True
    #     else:
    #         regAnts.inputs.transforms = ['Rigid', 'Affine']
    #         regAnts.inputs.transform_parameters = [(0.1,), (0.1,)]
    # #                regAnts.inputs.number_of_iterations = [[10,5,2,1]]*2
    #         regAnts.inputs.number_of_iterations = [[1000,500,250,100]]*2
    #         regAnts.inputs.metric = ['MI']*2
    #         regAnts.inputs.metric_weight = [1]*2
    #         regAnts.inputs.radius_or_number_of_bins = [32]*2
    #         regAnts.inputs.sampling_strategy = ['Regular']*2
    #         regAnts.inputs.sampling_percentage = [0.25]*2
    #         regAnts.inputs.convergence_threshold = [1.e-8]*2
    #         regAnts.inputs.convergence_window_size = [10]*2
    #         regAnts.inputs.smoothing_sigmas = [[3,2,1,0]]*2
    #         regAnts.inputs.sigma_units = ['mm']*2
    #         regAnts.inputs.shrink_factors = [[8,4,2,1]]*2
    #         regAnts.inputs.use_estimate_learning_rate_once = [True, True]
    #         regAnts.inputs.use_histogram_matching = [False]*2
    #     if save_transforms:
    #         mbFlow.connect(regAnts, 'composite_transform',
    #                        Sink, 'transforms.@affine_antsRegistration')


#=============================================================================
#-----------------------------------------------------------------------------
#
#   Surface workflows
#
#-----------------------------------------------------------------------------
#=============================================================================
if run_SurfFlows:

    mbFlow.connect(InputSubjects, 'subject', Surf, 'subject')
    mbFlow.connect(InputHemis, 'hemi', Surf, 'hemi')

    #-------------------------------------------------------------------------
    # Convert surfaces to VTK
    #-------------------------------------------------------------------------
    if not do_input_vtk:
        ConvertSurf = Node(name='Surface_to_vtk',
                           interface=Fn(function = surface_to_vtk,
                                        input_names=['surface_file',
                                                     'output_vtk'],
                                        output_names=['output_vtk']))
        mbFlow.connect(Surf, 'surface_files', ConvertSurf, 'surface_file')
        ConvertSurf.inputs.output_vtk = ''
    #-------------------------------------------------------------------------
    # Evaluation inputs: location and structure of atlas surfaces
    #-------------------------------------------------------------------------
    if (do_evaluate_surf_labels or surface_labels == 'manual') and do_label:
        SurfaceAtlas = Node(name='Surface_atlas',
                            interface=DataGrabber(infields=['subject','hemi'],
                                                  outfields=['atlas_file'],
                                                  sort_filelist=False))
        SurfaceAtlas.inputs.base_directory = subjects_path
        SurfaceAtlas.inputs.template = '%s/label/%s.labels.DKT31.' +\
                                       atlas_label_type + '.vtk'
        SurfaceAtlas.inputs.template_args['atlas_file'] = [['subject','hemi']]
    
        mbFlow.connect(InputSubjects, 'subject', SurfaceAtlas, 'subject')
        mbFlow.connect(InputHemis, 'hemi', SurfaceAtlas, 'hemi')

#=============================================================================
#
#   Surface labels
#
#=============================================================================
if run_SurfLabelFlow:

    SurfLabelFlow = Workflow(name='Surface_labels')

    #=========================================================================
    # Initialize labels with the DKT classifier atlas
    #=========================================================================
    if surface_labels == 'atlas' and not no_freesurfer:
        #---------------------------------------------------------------------
        # Label a brain with the DKT atlas using FreeSurfer's mris_ca_label
        #---------------------------------------------------------------------
        Classifier = Node(name='Label_surface_with_classifier',
                          interface=Fn(function = label_with_classifier,
                                       input_names=['subject',
                                                    'hemi',
                                                    'left_classifier',
                                                    'right_classifier',
                                                    'annot_file'],
                                       output_names=['annot_file']))
        SurfLabelFlow.add_nodes([Classifier])
        mbFlow.connect(InputSubjects, 'subject',
                       SurfLabelFlow, 'Label_surface_with_classifier.subject')
        mbFlow.connect(InputHemis, 'hemi',
                       SurfLabelFlow, 'Label_surface_with_classifier.hemi')
        left_classifier_file = 'lh.' + classifier_name + '.gcs'
        right_classifier_file = 'rh.' + classifier_name + '.gcs'
        left_classifier = retrieve_data(left_classifier_file, url,
                                        hashes, cache_env, cache)
        right_classifier = retrieve_data(right_classifier_file, url,
                                         hashes, cache_env, cache)
        Classifier.inputs.left_classifier = left_classifier
        Classifier.inputs.right_classifier = right_classifier
        Classifier.inputs.annot_file = ''

        #---------------------------------------------------------------------
        # Convert .annot file to .vtk format
        #---------------------------------------------------------------------
        Classifier2vtk = Node(name='annot_to_vtk',
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
                           SurfLabelFlow, 'annot_to_vtk.vtk_file')
        else:
            mbFlow.connect(ConvertSurf, 'output_vtk',
                           SurfLabelFlow, 'annot_to_vtk.vtk_file')
        mbFlow.connect(SurfLabelFlow, 'annot_to_vtk.output_vtk',
                       Sink, 'labels.@DKT_surface')
        plug = 'annot_to_vtk.output_vtk'
        plug1 = Classifier2vtk
        plug2 = 'output_vtk'

    #=========================================================================
    # Initialize labels with FreeSurfer's standard DK classifier atlas
    #=========================================================================
    elif surface_labels == 'freesurfer' and not no_freesurfer:
        FreeLabels = Node(name='FreeSurfer_annot_to_vtk',
                          interface=Fn(function = annot_to_vtk,
                                       input_names=['annot_file',
                                                    'vtk_file'],
                                       output_names=['labels',
                                                     'output_vtk']))
        SurfLabelFlow.add_nodes([FreeLabels])
        mbFlow.connect(Annot, 'annot_files',
                       SurfLabelFlow, 'FreeSurfer_annot_to_vtk.annot_file')
        if do_input_vtk:
            mbFlow.connect(Surf, 'surface_files',
                           SurfLabelFlow, 'FreeSurfer_annot_to_vtk.vtk_file')
        else:
            mbFlow.connect(ConvertSurf, 'output_vtk',
                           SurfLabelFlow, 'FreeSurfer_annot_to_vtk.vtk_file')
        mbFlow.connect(SurfLabelFlow, 'FreeSurfer_annot_to_vtk.output_vtk',
                       Sink, 'labels.@Free_surface')
        plug = 'FreeSurfer_annot_to_vtk.output_vtk'
        plug1 = FreeLabels
        plug2 = 'output_vtk'

        # #=========================================================================
        # # Initialize labels using multi-atlas registration
        # #=========================================================================
        # elif surface_labels == 'max_prob':
        #     #---------------------------------------------------------------------
        #     # Register surfaces to average template
        #     #---------------------------------------------------------------------
        #     Register = Node(name='Register_surface_template',
        #                     interface=Fn(function = register_template,
        #                                  input_names=['hemi',
        #                                               'sphere_file',
        #                                               'transform',
        #                                               'atlas_path',
        #                                               'template'],
        #                                  output_names=['transform']))
        #     SurfLabelFlow.add_nodes([Register])
        #     mbFlow.connect(InputHemis, 'hemi', SurfLabelFlow, 'Register_surface_template.hemi')
        #     mbFlow.connect(Surf, 'sphere_files',
        #                    SurfLabelFlow, 'Register_surface_template.sphere_file')
        #     Register.inputs.transform = 'sphere_to_' + free_template + '.reg'
        #     Register.inputs.atlas_path = atlas_path
        #     Register.inputs.template = free_template + '.tif'
        #     #---------------------------------------------------------------------
        #     # Register atlases to subject via template
        #     #---------------------------------------------------------------------
        #     Transform = MapNode(name='Label_surface',
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
        #     mbFlow.connect(InputSubjects, 'subject',
        #                    SurfLabelFlow, 'Label_surface.subject')
        #     mbFlow.connect(InputHemis, 'hemi',
        #                    SurfLabelFlow, 'Label_surface.hemi')
        #     SurfLabelFlow.connect(Register, 'transform', Transform, 'transform')
        #     #Transform.inputs.transform = 'sphere_to_' + template + '_template.reg'
        #     Transform.inputs.subjects_path = subjects_path
        #     Transform.inputs.atlas = atlas_list
        #     Transform.inputs.atlas_string = 'labels.DKT31.'+atlas_label_type
        #     #---------------------------------------------------------------------
        #     # Majority vote label
        #     #---------------------------------------------------------------------
        #     Vote = Node(name='Vote_surface_labels',
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
        #                          SurfLabelFlow, 'Vote_surface_labels.surface_file')
        #     else:
        #         mbFlow.connect(ConvertSurf, 'output_vtk',
        #                          SurfLabelFlow, 'Vote_surface_labels.surface_file')
        #     SurfLabelFlow.connect(Transform, 'output_file', Vote, 'annot_files')
        #     mbFlow.connect([(SurfLabelFlow, Sink,
        #                        [('Vote_surface_labels.maxlabel_file', 'labels.@max'),
        #                         ('Vote_surface_labels.labelcounts_file', 'labels.@counts'),
        #                         ('Vote_surface_labels.labelvotes_file', 'labels.@votes')])])
        #     plug = 'Vote_surface_labels.maxlabel_file'
        #     plug1 = Vote
        #     plug2 = 'maxlabel_file'

    #=========================================================================
    # Skip label initialization and process manual (atlas) labels
    #=========================================================================
    elif surface_labels == 'manual':
        ManualSurfLabels = Node(name='Manual_surface_labels',
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
        SurfLabelFlow.add_nodes([ManualSurfLabels])
        mbFlow.connect(SurfaceAtlas, 'atlas_file',
                       SurfLabelFlow, 'Manual_surface_labels.input_vtk')
        ManualSurfLabels.inputs.return_first = 'True'
        ManualSurfLabels.inputs.return_array = 'False'
        plug = 'Manual_surface_labels.input_vtk'
        plug1 = ManualSurfLabels
        plug2 = 'input_vtk'

    else:
        sys.exit('Label initialization improperly set.')

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
        mbFlow.connect(SurfaceAtlas, 'atlas_file',
                       EvalSurfLabels, 'labels_file1')
        mbFlow.connect(SurfLabelFlow, plug,
                       'EvalSurfLabels.labels_file2')
        #mbFlow.connect(EvalSurfLabels, 'overlap_file', Sink, 'evaluate_labels')

    #=========================================================================
    # Convert surface label numbers to volume label numbers
    #=========================================================================
    if not no_freesurfer:
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
        mbFlow.connect(InputHemis, 'hemi', SurfLabelFlow, 'Relabel_surface.hemi')
        RelabelSurface.inputs.old_labels = ''
        RelabelSurface.inputs.new_labels = ''
        RelabelSurface.inputs.output_file = ''
        mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                       Sink, 'labels.@surface')

#=============================================================================
#
#   Surface shape measurements
#
#=============================================================================
if run_WholeSurfShapeFlow:

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
                                                    'vtk_file',
                                                    'output_vtk'],
                                       output_names=['output_vtk']))
        WholeSurfShapeFlow.add_nodes([ConvexNode])
        mbFlow.connect(Surf, 'convexity_files',
                       WholeSurfShapeFlow, 'Convexity_to_vtk.surface_file')
        mbFlow.connect(ConvertSurf, 'output_vtk',
                       WholeSurfShapeFlow, 'Convexity_to_vtk.vtk_file')
        ConvexNode.inputs.output_vtk = ''
        mbFlow.connect(WholeSurfShapeFlow, 'Convexity_to_vtk.output_vtk',
                       Sink, 'shapes.@convexity')
    if do_thickness:
        ThickNode = Node(name='Thickness_to_vtk',
                         interface=Fn(function = curvature_to_vtk,
                                      input_names=['surface_file',
                                                   'vtk_file',
                                                   'output_vtk'],
                                      output_names=['output_vtk']))
        WholeSurfShapeFlow.add_nodes([ThickNode])
        mbFlow.connect(Surf, 'thickness_files',
                       WholeSurfShapeFlow, 'Thickness_to_vtk.surface_file')
        mbFlow.connect(ConvertSurf, 'output_vtk',
                       WholeSurfShapeFlow, 'Thickness_to_vtk.vtk_file')
        ThickNode.inputs.output_vtk = ''
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
                           [('output_vtk', 'Surface_area.surface_file'),
                            ('output_vtk', 'Travel_depth.surface_file'),
                            ('output_vtk', 'Geodesic_depth.surface_file'),
                            ('output_vtk', 'Curvature.surface_file')])])
    mbFlow.connect([(WholeSurfShapeFlow, Sink,
                       [('Surface_area.area_file', 'shapes.@surface_area'),
                        ('Travel_depth.depth_file', 'shapes.@travel_depth'),
                        ('Geodesic_depth.depth_file', 'shapes.@geodesic_depth'),
                        ('Curvature.mean_curvature_file', 'shapes.@mean_curvature')])])

#=============================================================================
#
#   Surface feature extraction
#
#=============================================================================
if run_SurfFeatureFlow:
    SurfFeatureFlow = Workflow(name='Surface_features')

    #=========================================================================
    # Folds
    #=========================================================================
    if do_folds:
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
        FoldsNode.inputs.min_fold_size = 50
        FoldsNode.inputs.tiny_depth = 0.001
        FoldsNode.inputs.save_file = True
        mbFlow.connect(SurfFeatureFlow, 'Folds.folds_file',
                       Sink, 'features.@folds')

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
        mbFlow.connect(InputHemis, 'hemi', SurfFeatureFlow, 'Sulci.hemi')
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
        FundiNode = Node(name='Fundi',
                         interface=Fn(function = extract_fundi,
                                      input_names=['folds',
                                                   'curv_file',
                                                   'depth_file',
                                                   'min_separation',
                                                   'erode_ratio',
                                                   'erode_min_size',
                                                   'sulci',
                                                   'save_file'],
                                      output_names=['fundi',
                                                    'n_fundi',
                                                    'fundi_file',
                                                    'n_fundi_all_folds',
                                                    'fundi_all_folds',
                                                    'fundi_all_folds_file']))
        SurfFeatureFlow.connect(FoldsNode, 'folds', FundiNode, 'folds')
        mbFlow.connect([(WholeSurfShapeFlow, SurfFeatureFlow,
                       [('Curvature.mean_curvature_file', 'Fundi.curv_file'),
                        ('Rescale_travel_depth.rescaled_scalars_file',
                         'Fundi.depth_file')])])
        FundiNode.inputs.min_separation = 10
        FundiNode.inputs.erode_ratio = 0.10
        FundiNode.inputs.erode_min_size = 10
        SurfFeatureFlow.connect(SulciNode, 'sulci', FundiNode, 'sulci')
        FundiNode.inputs.save_file = True
        mbFlow.connect(SurfFeatureFlow, 'Fundi.fundi_file',
                       Sink, 'features.@fundi')
        mbFlow.connect(SurfFeatureFlow, 'Fundi.fundi_all_folds_file',
                       Sink, 'features.@fundi_all_folds')

        if do_smooth_fundi:

            #-----------------------------------------------------------------
            # Compute likelihoods for smoothing fundi:
            #-----------------------------------------------------------------
            LikelihoodNode = Node(name='Likelihood',
                                  interface=Fn(function = compute_likelihood,
                                               input_names=['trained_file',
                                                            'depth_file',
                                                            'curvature_file',
                                                            'folds',
                                                            'save_file'],
                                               output_names=['likelihoods',
                                                             'likelihoods_file']))
            SurfFeatureFlow.add_nodes([LikelihoodNode])
            border_params_file = 'depth_curv_border_nonborder_parameters.pkl'
            border_params_path = retrieve_data(border_params_file, url,
                                               hashes, cache_env, cache)
            LikelihoodNode.inputs.trained_file = border_params_path
            mbFlow.connect([(WholeSurfShapeFlow, SurfFeatureFlow,
                               [('Rescale_travel_depth.rescaled_scalars_file',
                                 'Likelihood.depth_file'),
                                ('Curvature.mean_curvature_file',
                                 'Likelihood.curvature_file')])])
            SurfFeatureFlow.connect(FoldsNode, 'folds',
                                    LikelihoodNode, 'folds')
            LikelihoodNode.inputs.save_file = True
            #mbFlow.connect(SurfFeatureFlow, 'Likelihood.likelihoods_file',
            #               Sink, 'features.@likelihoods')

            #-----------------------------------------------------------------
            # Smooth fundi:
            #-----------------------------------------------------------------
            SmoothFundi = Node(name='Smooth_fundi',
                             interface=Fn(function = smooth_skeleton,
                                          input_names=['skeletons',
                                                       'bounds',
                                                       'vtk_file',
                                                       'likelihoods',
                                                       'wN_max',
                                                       'erode_again',
                                                       'save_file'],
                                          output_names=['smooth_skeletons',
                                                        'n_skeletons',
                                                        'skeletons_file']))
            SurfFeatureFlow.connect(FundiNode, 'fundi',
                                    SmoothFundi, 'skeletons')
            SurfFeatureFlow.connect(FoldsNode, 'folds', SmoothFundi, 'bounds')
            mbFlow.connect(WholeSurfShapeFlow, 'Curvature.mean_curvature_file',
                           SurfFeatureFlow, 'Smooth_fundi.vtk_file')
            SurfFeatureFlow.connect(LikelihoodNode, 'likelihoods',
                                    SmoothFundi, 'likelihoods')
            SmoothFundi.inputs.wN_max = 1.0
            SmoothFundi.inputs.erode_again = False
            SmoothFundi.inputs.save_file = True
            mbFlow.connect(SurfFeatureFlow, 'Smooth_fundi.skeletons_file',
                           Sink, 'features.@smooth_fundi')


#=============================================================================
#
#   Surface feature shapes
#
#=============================================================================
if run_SurfFlows:
    SurfFeatureShapeFlow = Workflow(name='Surface_feature_shapes')

    #=========================================================================
    # Decimate patches of surface mesh
    #=========================================================================
    if do_decimate:

        #---------------------------------------------------------------------
        # Decimate labeled regions
        #---------------------------------------------------------------------
        DecimateLabels = Node(name='Decimate_labels',
                              interface=Fn(function = decimate_file,
                                           input_names=['input_vtk',
                                                        'reduction',
                                                        'smooth_steps',
                                                        'save_vtk',
                                                        'output_vtk'],
                                           output_names=['output_vtk']))
        SurfFeatureShapeFlow.add_nodes([DecimateLabels])
        mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                       SurfFeatureShapeFlow, 'Decimate_labels.input_vtk')
        DecimateLabels.inputs.reduction = reduction
        DecimateLabels.inputs.smooth_steps = 100
        DecimateLabels.inputs.save_vtk = True
        DecimateLabels.inputs.output_vtk = ''

        #---------------------------------------------------------------------
        # Decimate sulci
        #---------------------------------------------------------------------
        if do_sulci:
            DecimateSulci = DecimateLabels.clone('Decimate_sulci')
            SurfFeatureShapeFlow.add_nodes([DecimateSulci])
            mbFlow.connect(SurfFeatureFlow, 'Sulci.sulci_file',
                           SurfFeatureShapeFlow, 'Decimate_sulci.input_vtk')

    #=========================================================================
    # Measure Laplace-Beltrami spectra
    #=========================================================================
    if do_spectra:

        #---------------------------------------------------------------------
        # Measure spectra of labeled regions
        #---------------------------------------------------------------------
        SpectraLabels = Node(name='Spectra_labels',
                             interface=Fn(function = spectrum_per_label,
                                          input_names=['vtk_file',
                                                       'n_eigenvalues',
                                                       'exclude_labels',
                                                       'normalization',
                                                       'area_file',
                                                       'largest_segment'],
                                          output_names=['spectrum_lists',
                                                        'label_list']))
        SurfFeatureShapeFlow.add_nodes([SpectraLabels])
        mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                       SurfFeatureShapeFlow, 'Spectra_labels.vtk_file')
        SpectraLabels.inputs.n_eigenvalues = n_eigenvalues
        SpectraLabels.inputs.exclude_labels = [0]
        SpectraLabels.inputs.normalization = "area"
        SpectraLabels.inputs.area_file = ""
        SpectraLabels.inputs.largest_segment = True
        mbFlow.connect(WholeSurfShapeFlow, 'Surface_area.area_file',
                       SurfFeatureShapeFlow, 'Spectra_labels.area_file')

        #---------------------------------------------------------------------
        # Measure spectra of sulci
        #---------------------------------------------------------------------
        if do_sulci:
            SpectraSulci = SpectraLabels.clone('Spectra_sulci')
            SurfFeatureShapeFlow.add_nodes([SpectraSulci])
            mbFlow.connect(SurfFeatureFlow, 'Sulci.sulci_file',
                           SurfFeatureShapeFlow, 'Spectra_sulci.vtk_file')
            SpectraSulci.inputs.exclude_labels = [-1]

    #=========================================================================
    # Measure Zernike moments
    #=========================================================================
    if do_zernike:

        #---------------------------------------------------------------------
        # Measure Zernike moments of labeled regions
        #---------------------------------------------------------------------
        ZernikeLabels = Node(name='Zernike_labels',
                             interface=Fn(function = zernike_moments_per_label,
                                          input_names=['vtk_file',
                                                       'order',
                                                       'exclude_labels',
                                                       'area_file',
                                                       'largest_segment'],
                                          output_names=['descriptors_lists',
                                                        'label_list']))
        SurfFeatureShapeFlow.add_nodes([ZernikeLabels])
        if do_decimate:
            SurfFeatureShapeFlow.connect(DecimateLabels, 'output_vtk',
                                         ZernikeLabels, 'vtk_file')
        else:
            mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                           SurfFeatureShapeFlow, 'Zernike_labels.vtk_file')
        ZernikeLabels.inputs.order = zernike_order
        ZernikeLabels.inputs.exclude_labels = [0]
        ZernikeLabels.inputs.largest_segment = True
        mbFlow.connect(WholeSurfShapeFlow, 'Surface_area.area_file',
                       SurfFeatureShapeFlow, 'Zernike_labels.area_file')

        #---------------------------------------------------------------------
        # Measure Zernike moments of sulci
        #---------------------------------------------------------------------
        if do_sulci:
            ZernikeSulci = ZernikeLabels.clone('Zernike_sulci')
            SurfFeatureShapeFlow.add_nodes([ZernikeSulci])
            ZernikeSulci.inputs.exclude_labels = [-1]
            if do_decimate:
                SurfFeatureShapeFlow.connect(DecimateSulci, 'output_vtk',
                                             ZernikeSulci, 'vtk_file')
            else:
                mbFlow.connect(SulciNode, 'sulci_file',
                               ZernikeSulci, 'vtk_file')


#=============================================================================
#
#   Surface feature shape tables
#
#=============================================================================
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
                                                     'labels_zernike',
                                                     'labels_zernike_IDs',
                                                     'sulci_zernike',
                                                     'sulci_zernike_IDs',
                                                     'exclude_labels',
                                                     'delimiter'],
                                        output_names=['label_table',
                                                      'sulcus_table',
                                                      'fundus_table']))
        mbFlow.add_nodes([ShapeTables])
        if do_label:
            mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                           ShapeTables, 'labels_or_file')
        else:
            ShapeTables.inputs.labels_or_file = []
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
        if do_register:
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

        # Laplace-Beltrami spectra:
        if do_spectra:
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

        # Zernike moments:
        if do_zernike:
            mbFlow.connect(SurfFeatureShapeFlow, 'Zernike_labels.descriptors_lists',
                           ShapeTables, 'labels_zernike')
            mbFlow.connect(SurfFeatureShapeFlow, 'Zernike_labels.label_list',
                           ShapeTables, 'labels_zernike_IDs')
            if do_sulci:
                mbFlow.connect(SurfFeatureShapeFlow, 'Zernike_sulci.descriptors_lists',
                               ShapeTables, 'sulci_zernike')
                mbFlow.connect(SurfFeatureShapeFlow, 'Zernike_sulci.label_list',
                               ShapeTables, 'sulci_zernike_IDs')
            else:
                ShapeTables.inputs.sulci_zernike = []
                ShapeTables.inputs.sulci_zernike_IDs = []
        else:
            ShapeTables.inputs.labels_zernike = []
            ShapeTables.inputs.sulci_zernike = []
            ShapeTables.inputs.labels_zernike_IDs = []
            ShapeTables.inputs.sulci_zernike_IDs = []

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
                                        input_names=['output_table',
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
                                        output_names=['output_table']))
        mbFlow.add_nodes([VertexTable])
        VertexTable.inputs.output_table = ''
        if do_label:
            mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                           VertexTable, 'labels_or_file')
        else:
            VertexTable.inputs.labels_or_file = []
        if do_sulci:
            mbFlow.connect(SurfFeatureFlow, 'Sulci.sulci',
                           VertexTable, 'sulci')
        else:
            VertexTable.inputs.sulci = []
        if do_fundi:
            mbFlow.connect(SurfFeatureFlow, 'Fundi.fundi',
                           VertexTable, 'fundi')
        else:
            VertexTable.inputs.fundi = []

        if run_VolFlows and do_register:
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
        mbFlow.connect(VertexTable, 'output_table', Sink, 'tables.@vertices')

    # #---------------------------------------------------------------------
    # # Apply RegFlows's affine transform to surface coordinates:
    # #---------------------------------------------------------------------
    # TransformPoints = Node(name='Transform_surface_points',
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
    # mbFlow.connect(SurfShapeFlow, 'Transform_surface_points.output_file',
    #                Sink, 'transforms.@points_to_template')


#=============================================================================
#-----------------------------------------------------------------------------
#
#   Volume workflows
#
#-----------------------------------------------------------------------------
#=============================================================================

#=============================================================================
#
#   Volume labels
#
#=============================================================================
if run_VolLabelFlow:
    VolLabelFlow = Workflow(name='Volume_labels')

    #=========================================================================
    # Fill cortical gray matter volume with surface labels
    #=========================================================================
    if do_fill_cortex:
        FillCortex = Node(name='Fill_cortex',
                          interface=Fn(function = fill_volume_with_surface_labels,
                                       input_names=['volume_mask',
                                                    'surface_files',
                                                    'output_file',
                                                    'binarize'],
                                       output_names=['output_file']))
        VolLabelFlow.add_nodes([FillCortex])
        if do_input_mask:
            mbFlow.connect(niftiMask, 'mask',
                           VolLabelFlow, 'Fill_cortex.volume_mask')
        else:
            mbFlow.connect(mgh_mask2nifti, 'out_file',
                           VolLabelFlow, 'Fill_cortex.volume_mask')
        if run_SurfFlows and run_SurfLabelFlow:
            mbFlow.connect(SurfLabelFlow, 'Relabel_surface.output_file',
                           VolLabelFlow, 'Fill_cortex.surface_files')
        else:
            sys.exit('No input surface file specified for Fill_cortex.')
        FillCortex.inputs.output_file = ''
        FillCortex.inputs.binarize = False
        mbFlow.connect(VolLabelFlow, 'Fill_cortex.output_file',
                       Sink, 'labels.@filled_cortex')

    #=========================================================================
    # Label using registered volume
    #=========================================================================
    if do_register_labels:

        # Retrieve full atlas path(s):
        RetrieveVolAtlas = Node(name='Retrieve_volume_atlas',
                                interface=Fn(function = retrieve_data,
                                             input_names=['data_file',
                                                          'url',
                                                          'hashes',
                                                          'cache_env',
                                                          'cache',
                                                          'return_missing'],
                                             output_names=['data_path']))
        VolLabelFlow.add_nodes([RetrieveVolAtlas])
        mbFlow.connect(InputVolAtlases, 'atlas',
                       VolLabelFlow, 'Retrieve_volume_atlas.data_file')
        RetrieveVolAtlas.inputs.url = url
        RetrieveVolAtlas.inputs.hashes = hashes
        RetrieveVolAtlas.inputs.cache_env = cache_env
        RetrieveVolAtlas.inputs.cache = cache
        RetrieveVolAtlas.inputs.return_missing = True


        # Inverse transform atlas label volume to subject via template
        LabelVolume = Node(name='Label_volume',
                           interface=Fn(function = WarpImageMultiTransform,
                                        input_names=['source',
                                                     'target',
                                                     'output',
                                                     'interp',
                                                     'xfm_stem',
                                                     'affine_transform',
                                                     'nonlinear_transform',
                                                     'inverse',
                                                     'affine_only'],
                                        output_names=['output']))
        VolLabelFlow.add_nodes([LabelVolume])
        VolLabelFlow.connect(RetrieveVolAtlas, 'data_path',
                             LabelVolume, 'source')
        if do_input_nifti:
            mbFlow.connect(niftiBrain, 'nifti',
                           VolLabelFlow, 'Label_volume.target')
        else:
            mbFlow.connect(mgh2nifti, 'out_file',
                           VolLabelFlow, 'Label_volume.target')
        LabelVolume.inputs.output = ''
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
                       Sink, 'labels.@registered')

        #=====================================================================
        # Mask cortical volume labels
        #=====================================================================
        if do_mask_cortex:
            MaskVolume = Node(name='Mask_volume',
                              interface=Fn(function = ImageMath,
                                           input_names=['volume1',
                                                        'volume2',
                                                        'operator',
                                                        'output_file'],
                                           output_names=['output_file']))
            VolLabelFlow.add_nodes([MaskVolume])
            if do_input_mask:
                mbFlow.connect(niftiMask, 'nifti.@mask',
                               VolLabelFlow, 'Mask_volume.volume1')
            else:
                mbFlow.connect(mgh_mask2nifti, 'out_file',
                               VolLabelFlow, 'Mask_volume.volume1')
            VolLabelFlow.connect(LabelVolume, 'output',
                                 MaskVolume, 'volume2')
            MaskVolume.inputs.operator = 'm'
            MaskVolume.inputs.output_file = ''
            mbFlow.connect(VolLabelFlow, 'Mask_volume.output_file',
                           Sink, 'labels.@masked_registered')

        #=====================================================================
        # Combine cortical and noncortical volume labels
        #=====================================================================
        if do_fill_cortex:
            CombineLabels = Node(name='Combine_labels',
                                 interface=Fn(function = overwrite_volume_labels,
                                              input_names=['source',
                                                           'target',
                                                           'output_file',
                                                           'ignore_labels',
                                                           'replace'],
                                              output_names=['output_file']))
            VolLabelFlow.connect(FillCortex, 'output_file',
                                 CombineLabels, 'source')
            VolLabelFlow.connect(LabelVolume, 'output',
                                 CombineLabels, 'target')
            CombineLabels.inputs.output_file = ''
            CombineLabels.inputs.ignore_labels = [0]
            CombineLabels.inputs.replace = True
            mbFlow.connect(VolLabelFlow, 'Combine_labels.output_file',
                           Sink, 'labels.@combined')

    #=========================================================================
    # Evaluate label volume overlaps
    #=========================================================================
    if do_evaluate_vol_labels:

        #---------------------------------------------------------------------
        # Evaluation inputs: location and structure of atlas volumes
        #---------------------------------------------------------------------
        VolAtlas = Node(name='Volume_atlas',
                        interface=DataGrabber(infields=['subject'],
                                              outfields=['atlas_vol_file'],
                                              sort_filelist=False))
        VolLabelFlow.add_nodes([VolAtlas])
        VolAtlas.inputs.base_directory = subjects_path
        VolAtlas.inputs.template = '%s/mri/labels.DKT31.manual.nii.gz'
        VolAtlas.inputs.template_args['atlas_vol_file'] = [['subject']]
        mbFlow.connect(InputSubjects, 'subject',
                       VolLabelFlow, 'Volume_atlas.subject')
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
        if do_register_labels and do_fill_cortex:
            EvalVolLabels.inputs.labels = label_numbers
        elif do_fill_cortex:
            EvalVolLabels.inputs.labels = cortex_numbers
        VolLabelFlow.connect(VolAtlas, 'atlas_vol_file', EvalVolLabels, 'file2')
        VolLabelFlow.connect(FillCortex, 'output_file', EvalVolLabels, 'file1')
        #mbFlow.connect(VolLabelFlow, 'Evaluate_volume_labels.out_file',
        #               Sink, 'evaluate_labels_volume')


#=============================================================================
#
#   Volume feature shapes
#
#=============================================================================
if run_VolShapeFlow:
    VolShapeFlow = Workflow(name='Volume_feature_shapes')

    #=========================================================================
    # Measure volume of each region of a labeled image file
    #=========================================================================
    LabelVolumes = Node(name='Label_volumes',
                          interface=Fn(function = volume_per_label,
                                       input_names=['labels',
                                                    'input_file'],
                                       output_names=['labels_volumes']))
    LabelVolumes.inputs.labels = label_numbers
    #-------------------------------------------------------------------------
    # Volumes of the registered labels:
    #-------------------------------------------------------------------------
    if do_register_labels:
        VolShapeFlow.add_nodes([LabelVolumes])
        mbFlow.connect(VolLabelFlow, 'Label_volume.output',
                       VolShapeFlow, 'Label_volumes.input_file')
        #---------------------------------------------------------------------
        # Volumes of the cortex-masked registered labels:
        #---------------------------------------------------------------------
        if do_mask_cortex:
            MaskedVolumes = LabelVolumes.clone('Masked_volumes')
            VolShapeFlow.add_nodes([MaskedVolumes])
            mbFlow.connect(VolLabelFlow, 'Mask_volume.output_file',
                           VolShapeFlow, 'Masked_volumes.input_file')
        #---------------------------------------------------------------------
        # Volumes of the combined registered and cortex-filled surface labels:
        #---------------------------------------------------------------------
        if do_fill_cortex:
            CombinedVolumes = LabelVolumes.clone('Combined_volumes')
            VolShapeFlow.add_nodes([CombinedVolumes])
            mbFlow.connect(VolLabelFlow, 'Combine_labels.output_file',
                           VolShapeFlow, 'Combined_volumes.input_file')
    #-------------------------------------------------------------------------
    # Volumes of the cortex-filled surface labels:
    #-------------------------------------------------------------------------
    if do_fill_cortex:
        FilledVolumes = LabelVolumes.clone('Filled_volumes')
        VolShapeFlow.add_nodes([FilledVolumes])
        mbFlow.connect(VolLabelFlow, 'Fill_cortex.output_file',
                       VolShapeFlow, 'Filled_volumes.input_file')

    #=========================================================================
    # Create tables to save the volume measures
    #=========================================================================
    if do_vol_table:
        VolumesTable = Node(name='Label_volumes_table',
                            interface=Fn(function = write_columns,
                                         input_names=['columns',
                                                      'column_names',
                                                      'delimiter',
                                                      'quote',
                                                      'input_table',
                                                      'output_table'],
                                         output_names=['output_table']))
        VolumesTable.inputs.column_names = ['label', 'volume']
        VolumesTable.inputs.delimiter = ','
        VolumesTable.inputs.quote = True
        VolumesTable.inputs.input_table = ''
        #---------------------------------------------------------------------
        # Volumes of the registered labels:
        #---------------------------------------------------------------------
        if do_register_labels:
            mbFlow.add_nodes([VolumesTable])
            VolShapeFlow.connect(LabelVolumes, 'labels_volumes',
                                 VolumesTable, 'columns')
            s = 'volumes_of_registered_volume_labels.csv'
            VolumesTable.inputs.output_table = s
            mbFlow.connect(VolumesTable, 'output_table',
                           Sink, 'tables.@registered_label_volumes')
            #-----------------------------------------------------------------
            # Volumes of the masked registered labels:
            #-----------------------------------------------------------------
            if do_mask_cortex:
                MaskedTable = VolumesTable.clone('Masked_label_volumes_table')
                VolShapeFlow.add_nodes([MaskedTable])
                VolShapeFlow.connect(MaskedVolumes, 'labels_volumes',
                                     MaskedTable, 'columns')
                s = 'volumes_of_cortical_surface_masked_volume_labels.csv'
                MaskedTable.inputs.output_table = s
                mbFlow.connect(MaskedTable, 'output_table',
                               Sink, 'tables.@masked_label_volumes')
            #-----------------------------------------------------------------
            # Volumes of the combined registered and filled labels:
            #-----------------------------------------------------------------
            if do_fill_cortex:
                CombinedTable = VolumesTable.clone('Combined_label_volumes_table')
                VolShapeFlow.add_nodes([CombinedTable])
                VolShapeFlow.connect(CombinedVolumes, 'labels_volumes',
                                     CombinedTable, 'columns')
                s = 'volumes_of_cortical_surface_and_noncortical_volume_labels.csv'
                CombinedTable.inputs.output_table = s
                mbFlow.connect(CombinedTable, 'output_table',
                               Sink, 'tables.@combined_label_volumes')
        #---------------------------------------------------------------------
        # Volumes of the filled labels:
        #---------------------------------------------------------------------
        if do_fill_cortex:
            FilledTable = VolumesTable.clone('Filled_label_volumes_table')
            VolShapeFlow.add_nodes([FilledTable])
            VolShapeFlow.connect(FilledVolumes, 'labels_volumes',
                                 FilledTable, 'columns')
            s = 'volumes_of_cortical_surface_labels.csv'
            FilledTable.inputs.output_table = s
            mbFlow.connect(FilledTable, 'output_table',
                           Sink, 'tables.@filled_label_volumes')


#=============================================================================
#-----------------------------------------------------------------------------
#
#   Run workflows
#
#-----------------------------------------------------------------------------
#=============================================================================
if __name__== '__main__':

    #-------------------------------------------------------------------------
    # Generate a visual graph:
    #-------------------------------------------------------------------------
    if graph_vis:
        if graph_vis == 'exec':
            mbFlow.write_graph(graph2use=graph_vis, simple_form=False)
        else:
            mbFlow.write_graph(graph2use=graph_vis)

    #-------------------------------------------------------------------------
    # Run (HTCondor) cluster processes, such as on the Mindboggler cluster:
    #-------------------------------------------------------------------------
    if cluster:
        mbFlow.run(plugin='CondorDAGMan')
    #-------------------------------------------------------------------------
    # Run multiple processes or not:
    #-------------------------------------------------------------------------
    else:
        if n_processes:
            if n_processes > 1:
                mbFlow.run(plugin='MultiProc',
                           plugin_args={'n_procs': n_processes})
            else:
                mbFlow.run()
        else:
            mbFlow.run()
        # # Default is to use all processors:
        #else:
        #    mbFlow.run(plugin='MultiProc')

"""
# Script for running Mindboggle on the Mindboggle-101 set:
import os
from mindboggle.utils.io_table import read_columns
from mindboggle.utils.utils import execute

out_path = '/homedir/Mindboggle101/mindboggled/'
atlas_list_file = '/homedir/Mindboggle101/code/mindboggle101_atlases.txt'
atlas_list = read_columns(atlas_list_file, 1)[0]

for atlas in atlas_list:
    if 'MMRR-21-' in atlas:
        cmd = ' '.join(['python mindboggler.py', '-o', out_path, atlas])
        execute(cmd)

for atlas in atlas_list:
    if 'NKI-RS-' in atlas:
        cmd = ' '.join(['python mindboggler.py', '-o', out_path, atlas])
        execute(cmd)

for atlas in atlas_list:
    if 'NKI-TRT-' in atlas:
        cmd = ' '.join(['python mindboggler.py', '-o', out_path, atlas])
        execute(cmd)

for atlas in atlas_list:
    if 'OASIS-TRT-20' in atlas:
        cmd = ' '.join(['python mindboggler.py', '-o', out_path, atlas])
        execute(cmd)

for atlas in atlas_list:
    if 'HLN-' in atlas or 'Twins-' in atlas or 'Colin' in atlas or 'After' in atlas or '3T7T' in atlas:
        cmd = ' '.join(['python mindboggler.py', '-o', out_path, atlas])
        execute(cmd)

"""
