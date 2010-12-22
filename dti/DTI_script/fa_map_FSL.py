"""
DTI tool - fa_map_FSL.py

It fits the tensor to DTI data and compute FA, eigenvalues and eigenvectors maps using FSL's tool DTIFIT.

INPUT: 
- Subject directory (e.g. /mind_cart/denis/50241/). It must already contain:
  * sub_dir/fsl_data/data.nii              :the DTI data
  * sub_dir/fsl_data/nodif_brain_mask.nii  :binary mask of the brain (NOTE: it must be binary)
  * sub_dir/fsl_data/bvals                 :text file containing b values in FSL format
  * sub_dir/fsl_data/bvecs                 :gradient directions in FSL format (NOTE: directions must be normalized)
- Output directory (optional. Default is sub_dir/fsl_data.fa_map/)
- Basenare (optional. Default is last directory name in sub_dir, that is usually the subject ID)

OUTPUT
The output directory will contain:
  * <basename>_V1 - 1st eigenvector
  * <basename>_V2 - 2nd eigenvector
  * <basename>_V3 - 3rd eigenvector
  * <basename>_L1 - 1st eigenvalue
  * <basename>_L2 - 2nd eigenvalue
  * <basename>_L3 - 3rd eigenvalue
  * <basename>_MD - mean diffusivity
  * <basename>_FA - fractional anisotropy
  * <basename>_MO - mode of the anisotropy (oblate ~ -1; isotropic ~ 0; prolate ~ 1)
  * <basename>_S0 - raw T2 signal with no diffusion weighting

Example:
fa_map_fsl /mind_cart/denis/50219 [/mind_cart/denis/50219_results result]

"""

# IMPORT
import os, sys

# SETTINGS
# Default settings
file_sep      = '/'
output_subdir = 'fsl_data.fa_map/'
input_dir     = 'fsl_data/'
bvals_file    = 'bvals'
bvecs_file    = 'bvecs'
data_file     = 'data.nii.gz'
mask_file     = 'nodif_brain_mask.nii.gz'



# STEP 1: reading input data
if len(sys.argv)<2: # NOTE: the first sys.argv element is the script name
	print 'WARNING: not enough input argument. You must specify at least the subject directory.'
	raise 'Input data error'
else:
	sub_dir = sys.argv[1];
	
	# Making sure that sub_dir finishes with /
	if sub_dir[len(sub_dir)-1]==file_sep:
		sub_dir = sub_dir[:len(sub_dir)-1]
	
	# Retriving the last directory name (i.e. subject ID)
	path_words = sub_dir.split(file_sep)
	basename = path_words[len(path_words)-1]
	sub_dir = sub_dir + file_sep
	
if len(sys.argv)>=3: # A new output directory is given
	output_dir = sys.argv[2]
else:
	output_dir = sub_dir + output_subdir

if len(sys.argv)>=3: # A new basename is given
	basename = sys.argv[3]
	# NOTE: if a basename is not given the last directory name in sub_dir is used
	
# Option recap
print 'FA MAP TOOL'
print ' subject path: %s' %(sub_dir)
print ' output path:  %s' %(output_dir)
print ' subject ID:   %s' %(basename)

# STEP 2: checking for the presence of the data
print 'Checking data...'
error_flag = 0
if not(os.path.isfile(sub_dir + input_dir + bvals_file)):
	print 'WARNING: file %s is not present in %s' %( bvals_file, sub_dir + input_dir)
	error_flag = 1
	
if not(os.path.isfile(sub_dir + input_dir + bvecs_file)):
	print 'WARNING: file %s is not present in %s' %( bvecs_file, sub_dir + input_dir)
	error_flag = 1

if not(os.path.isfile(sub_dir + input_dir + data_file)):
	print 'WARNING: file %s is not present in %s' %( data_file, sub_dir + input_dir)
	error_flag = 1
	
if not(os.path.isfile(sub_dir + input_dir + mask_file)):
	print 'WARNING: file %s is not present in %s' %( mask_file, sub_dir + input_dir)
	error_flag = 1
	
if error_flag:
	print 'Cannot perform the script. Please check the input data.'
	raise 'Input data error'

# STEP 3: calling DTIFIT tool
cmd = 'mkdir ' + output_dir
os.system(cmd)
cmd = 'chmod 777 ' + output_dir
os.system(cmd)

cmd = 'dtifit -k ' + sub_dir + input_dir + data_file + ' -m ' + sub_dir + input_dir + mask_file + ' -o ' + output_dir + basename + ' -r ' + sub_dir + input_dir + bvecs_file +  ' -b ' + sub_dir + input_dir + bvals_file
os.system(cmd)

cmd = 'chmod 666 ' + output_dir + basename + '*.*'
os.system(cmd)