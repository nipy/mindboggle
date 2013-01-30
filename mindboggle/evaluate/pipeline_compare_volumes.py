#!/usr/bin/env python
"""
This is a Nipype pipeline for comparing nibabel-readable image volumes.


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""
import os

#=============================================================================
# Setup: import libraries, set file paths, and initialize main workflow
#=============================================================================
#-----------------------------------------------------------------------------
# Paths
#-----------------------------------------------------------------------------
# images_path is the beginning of the path not in the text of image_list file
data_path = '/Users/arno/Data'
results_path = '/desk/Results'
temp_path = '/desk'
#-----------------------------------------------------------------------------
# Data to run
#-----------------------------------------------------------------------------
process_phantoms = False  # phantom or human data?
process_DTI = True  # structural or diffusion data?
#-----------------------------------------------------------------------------
# Steps to run
#-----------------------------------------------------------------------------
do_register_images_to_ref_images = True
do_threshold_images = True
do_compute_image_similarities = True
do_compare_image_histograms = True
do_compute_image_overlaps = False  # Not as useful
#-----------------------------------------------------------------------------
# Settings
#-----------------------------------------------------------------------------
run_bet = False
if process_phantoms:
    max_angle = 15
    if process_DTI:
        temp_path = os.path.join(temp_path, 'workspace_DTI_phantoms')
    else:
        temp_path = os.path.join(temp_path, 'workspace_phantoms')
else:
    max_angle = 45
    if process_DTI:
        temp_path = os.path.join(temp_path, 'workspace_DTI_humans')
    else:
        run_bet = True
        temp_path = os.path.join(temp_path, 'workspace_humans')
if process_DTI:
    threshold_value = 0.3
    interp = 'nearestneighbour'
else:
    threshold_value = 0.1
    interp = 'trilinear'
#-----------------------------------------------------------------------------
# Lists of images to process
#-----------------------------------------------------------------------------
if process_phantoms:
    if process_DTI:
        output_path = 'DTI_phantoms'
        images_path = os.path.join(data_path, 'Phantom_DTI_Updated20121214')
        reg_image_list = 'Phantom_DTI_1stvol_rotated.txt'
        image_list = 'Phantom_DTI_FA_rotated.txt'
    else:
        output_path = 'ADNI_phantoms'
        images_path = os.path.join(data_path, 'Phantom_ADNI_Updated20121213')
        reg_image_list = 'Phantom_ADNI.txt'
        image_list = 'Phantom_ADNI.txt'
else:
    if process_DTI:
        output_path = 'Human_DTI_retests'
        images_path = os.path.join(data_path, 'Human_DTI_Updated20130123')
        reg_image_list = 'Human_retest_1stvol_files.txt'
        image_list = 'Human_retest_FA_files.txt'
    else:
        output_path = 'Human_retests'
        images_path = os.path.join(data_path, 'Human_STR_Updated20130123')
        reg_image_list = 'Human_retest_files.txt'
        image_list = 'Human_retest_files.txt'
image_list = os.path.join(images_path, image_list)
reg_image_list = os.path.join(images_path, reg_image_list)
output_path = os.path.join(results_path, output_path)


#-----------------------------------------------------------------------------
# Import system and nipype Python libraries
#-----------------------------------------------------------------------------
from nipype.pipeline.engine import Workflow, Node
from nipype.interfaces.utility import Function as Fn
from nipype.interfaces.io import DataSink
#from nipype.interfaces.utility import IdentityInterface
#-----------------------------------------------------------------------------
# Import Mindboggle Python libraries
#-----------------------------------------------------------------------------
from mindboggle.measure.measure_functions import pairwise_vector_distances
from mindboggle.evaluate.compare_images import compute_image_histograms, \
    compute_image_similarities, compute_image_overlaps, \
    register_images_to_ref_images, apply_transforms, threshold_images
#-------------------------------------------------------------------------------
# Initialize workflows
#-------------------------------------------------------------------------------
Flow = Workflow(name='Image_comparison_workflow')
Flow.base_dir = temp_path
if not os.path.isdir(temp_path):
    os.makedirs(temp_path)
#-----------------------------------------------------------------------------
# Inputs and Outputs
#-----------------------------------------------------------------------------
fid = open(image_list)
file_list = fid.read()
file_list = file_list.splitlines()
file_list = [x.strip() for x in file_list if len(x)]
fid_reg = open(reg_image_list)
reg_file_list = fid_reg.read()
reg_file_list = reg_file_list.splitlines()
reg_file_list = [x.strip() for x in reg_file_list if len(x)]
#Info = Node(name = 'Inputs',
#            interface = IdentityInterface(fields=['files']))
#Info.iterables = ([('files', file_list)])
Sink = Node(DataSink(), name = 'Results')
Sink.inputs.base_directory = output_path
Sink.inputs.container = 'Results'
if not os.path.isdir(output_path):  os.makedirs(output_path)

#=============================================================================
#   Preprocessing
#=============================================================================
def run_bet_on_files(infiles, indirectory, f_value=0.5):
    """
    Run FSL's bet (brain extraction tool).

    Parameters
    ----------
    infiles : names of input files
    indirectory : name of input directory
    f_value : float

    """
    import os
    from nipype.interfaces.base import CommandLine

    outfiles = []
    for infile in infiles:
        outfile = os.path.join(os.getcwd(), 'brain_' + os.path.basename(infile))
        outfiles.append(outfile)

        cli = CommandLine(command = 'bet2')
        file = os.path.join(indirectory, infile)
        cli.inputs.args = ' '.join([file, outfile, 'f_value' + str(f_value)])
        cli.cmdline
        cli.run()

    return outfiles

if run_bet:

    bet = Node(name = 'Extract_brains',
               interface = Fn(function = run_bet_on_files,
                              iterfield=['infiles'],
                              input_names = ['infiles',
                                             'indirectory',
                                             'f_value'],
                              output_names = ['outfiles']))
    Flow.add_nodes([bet])
    bet.inputs.infiles = file_list
    bet.inputs.indirectory = images_path
    bet.inputs.f_value = 0.25
    Flow.connect([(bet, Sink, [('outfiles', 'brains')])])

#=============================================================================
#   Comparisons
#=============================================================================
#-------------------------------------------------------------------------------
# Compare image histograms
# The images from which the histograms were derived do not need to be coregistered.
#-------------------------------------------------------------------------------
if do_compare_image_histograms:
    compute_histograms = Node(name = 'Compute_histograms',
                              interface = Fn(function = compute_image_histograms,
                                             iterfield=['infiles'],
                                             input_names = ['infiles',
                                                            'indirectory',
                                                            'nbins',
                                                            'threshold'],
                                             output_names = ['histogram_values']))
    Flow.add_nodes([compute_histograms])
    if run_bet:
        Flow.connect([(bet, compute_histograms, [('outfiles','infiles')])])
    else:
        compute_histograms.inputs.infiles = file_list
    compute_histograms.inputs.indirectory = images_path
    compute_histograms.inputs.nbins = 100
    compute_histograms.inputs.threshold = 0

    compare_histograms = Node(name = 'Compare_histograms',
                              interface = Fn(function = pairwise_vector_distances,
                                             input_names = ['vectors',
                                                            'save_file',
                                                            'normalize'],
                                             output_names = ['vector_distances',
                                                             'outfile']))
    Flow.add_nodes([compare_histograms])
    Flow.connect([(compute_histograms, compare_histograms,
                   [('histogram_values','vectors')])])
    compare_histograms.inputs.save_file = True
    compare_histograms.inputs.normalize = True

    Flow.connect([(compare_histograms, Sink, [('outfile', 'histograms')])])

#-------------------------------------------------------------------------------
# Register each image to a reference image
#-------------------------------------------------------------------------------
if do_register_images_to_ref_images:
    register = Node(name = 'Register',
                    interface = Fn(function = register_images_to_ref_images,
                                   input_names = ['files',
                                                  'ref_file_index',
                                                  'max_angle',
                                                  'directory'],
                                   output_names = ['outfiles']))
    Flow.add_nodes([register])
    if run_bet:
        Flow.connect([(bet, register, [('outfiles','files')])])
    else:
        register.inputs.files = reg_file_list
    register.inputs.ref_file_index = 1
    register.inputs.max_angle = max_angle
    register.inputs.directory = images_path
    #Flow.connect([(register, Sink, [('outfiles', 'registrations.@transforms')])])

    transform = Node(name = 'Transform',
                    interface = Fn(function = apply_transforms,
                                   input_names = ['files',
                                                  'ref_file_index',
                                                  'transform_files',
                                                  'directory',
                                                  'interp'],
                                   output_names = ['outfiles']))
    Flow.add_nodes([transform])
    transform.inputs.files = file_list
    transform.inputs.ref_file_index = 1
    Flow.connect([(register, transform, [('outfiles', 'transform_files')])])
    transform.inputs.directory = images_path
    transform.inputs.interp = interp
    Flow.connect([(transform, Sink, [('outfiles', 'registrations.@images')])])

#-------------------------------------------------------------------------------
# Threshold images
#-------------------------------------------------------------------------------
if do_threshold_images:
    threshold = Node(name = 'Threshold',
                     interface = Fn(function = threshold_images,
                                    input_names = ['files',
                                                   'threshold_value',
                                                   'save_files'],
                                    output_names = ['outfiles']))
    Flow.add_nodes([threshold])
    Flow.connect([(transform, threshold, [('outfiles', 'files')])])
    threshold.inputs.threshold_value = threshold_value
    threshold.inputs.save_files = True
    Flow.connect([(threshold, Sink, [('outfiles', 'thresholds')])])

#-------------------------------------------------------------------------------
# Compute image similarities
#-------------------------------------------------------------------------------
if do_compute_image_similarities:
    similarity = Node(name = 'Similarity',
                      interface = Fn(function = compute_image_similarities,
                                    input_names = ['files',
                                                   'masks',
                                                   'metric',
                                                   'save_file'],
                                    output_names = ['pairwise_similarities',
                                                    'outfile']))
    Flow.add_nodes([similarity])
    Flow.connect([(transform, similarity, [('outfiles', 'files')])])
    Flow.connect([(threshold, similarity, [('outfiles', 'masks')])])
    similarity.inputs.metric = 'cc'
    similarity.inputs.save_file = True
    Flow.connect([(similarity, Sink, [('outfile', 'similarities')])])

#-------------------------------------------------------------------------------
# Compute image overlaps
#-------------------------------------------------------------------------------
if do_compute_image_overlaps:
    overlaps = Node(name = 'Overlaps',
                     interface = Fn(function = compute_image_overlaps,
                                    input_names = ['files',
                                                   'list_of_labels',
                                                   'save_file'],
                                    output_names = ['pairwise_overlaps',
                                                    'outfile']))
    Flow.add_nodes([overlaps])
    Flow.connect([(threshold, overlaps, [('outfiles', 'files')])])
    overlaps.inputs.list_of_labels = [1]
    overlaps.inputs.save_file = True
    Flow.connect([(overlaps, Sink, [('outfile', 'overlaps')])])

##############################################################################
if __name__== '__main__':
    #Flow.write_graph(graph2use='flat')
    #Flow.write_graph(graph2use='hierarchical')
    Flow.run()
