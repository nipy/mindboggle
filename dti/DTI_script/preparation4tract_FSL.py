"""

DTI tool - preparation4tract_FSL.py

It performs FSL bedpostx tool.
INPUT: 
- Subject directory (e.g. /mind_cart/denis/50241/). It must already contain:
  * sub_dir/fsl_data/data.nii :the DTI data
  * sub_dir/fsl_data/nodif_brain_mask.nii :binary mask of the brain (NOTE: it must be binary)
  * sub_dir/fsl_data/bvals :text file containing b values in FSL format
  * sub_dir/fsl_data/bvecs :gradient directions in FSL format (NOTE: directions must be normalized)

OUTPUT
The output directory (subject_directory/fsl_data.bedpostX) will contain:
  * merged_th<i>samples - 4D volume - Samples from the distribution on theta
  * merged_ph<i>samples - Samples from the distribution on phi o theta and phi together represent the principal diffusion direction in spherical polar co-ordinates
  * merged_f<i>samples - 4D volume - Samples from the distribution on anisotropic volume fraction (see technical report).
  * mean_th<i>samples - 3D Volume - Mean of distribution on theta
  * mean_ph<i>samples - 3D Volume - Mean of distribution on phi
  * mean_f<i>samples - 3D Volume - Mean of distribution on f anisotropy o Note that in each voxel, fibres are ordered according to a decreasing mean f-value
  * dyads<i> - Mean of PDD distribution in vector form. Note that this file can be loaded into fslview for easy viewing of diffusion directions
  * nodif_brain - brain extracted version of nodif - copied from input directory
  * nodif_brain_mask - binary mask created from nodif_brain - copied from input directory
  
Example:
preparation4tract_FSL /mind_cart/denis/50219
"""
# IMPORT
import os, sys

# SETTINGS
bedpostx_options = '-n 2' #i.e. set the maximum number of accepted fibre's compartment per voxel.
input_dir     = 'fsl_data/'
bvals_file    = 'bvals'
bvecs_file    = 'bvecs'
data_file     = 'data.nii.gz'
mask_file     = 'nodif_brain_mask.nii.gz'

# STEP 1: reading input argument
parent_dir = sys.argv[1]
if parent_dir[len(parent_dir)-1] != '/':
	parent_dir = parent_dir + '/'

# STEP 2: checking for the presence of the data
print 'Checking data...'
error_flag = 0
if not(os.path.isfile(parent_dir + input_dir + bvals_file)):
	print 'WARNING: file %s is not present in %s' %( bvals_file, parent_dir + input_dir)
	error_flag = 1
	
if not(os.path.isfile(parent_dir + input_dir + bvecs_file)):
	print 'WARNING: file %s is not present in %s' %( bvecs_file, parent_dir + input_dir)
	error_flag = 1

if not(os.path.isfile(parent_dir + input_dir + data_file)):
	print 'WARNING: file %s is not present in %s' %( data_file, parent_dir + input_dir)
	error_flag = 1
	
if not(os.path.isfile(parent_dir + input_dir + mask_file)):
	print 'WARNING: file %s is not present in %s' %( mask_file, parent_dir + input_dir)
	error_flag = 1
	
if error_flag:
	print 'Cannot perform the script. Please check the input data.'
	raise

# STEP 3: running bedpostx tool
print 'Running bedpostx...'
cmd = 'bedpostx ' + parent_dir + input_dir + ' ' + bedpostx_options
os.system(cmd)
cmd = 'chmod 777 ' + parent_dir + input_dir[0:len(input_dir)-1] + '.bedpostX'
os.system(cmd)
cmd = 'chmod 666 ' + parent_dir + input_dir[0:len(input_dir)-1] + '.bedpostX/*.*'
os.system(cmd)
