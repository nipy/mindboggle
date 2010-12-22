"""
DTI tool - format_data_FSL.py

The script prepares the data in the FSL format for bedpostX and dtifit fsl's tools.

INPUT: 
- Subject directory (e.g. /mind_cart/denis/50241/). It must already contain:
  * sub_dir/original_data/DTIdata.nii   :the DTI data
  * sub_dir/original_data/DTI_mask.nii  :binary mask of the brain (NOTE: it must be binary)
- Output directory (optional. Default is sub_dir/fsl_data/)
- bvals_file (optional. If not given the script import the bvals_fsl file already stored for Columbia DTI data):the file containing the b values (NOTE: it must already be in fsl format)
- bvecs_file (optional. If not given the script import the bvecs_fsl file already stored for Columbia DTI data):the file containing the gradient directions (NOTE: it must already be in fsl format)

OUTPUT
The output directory will contain:
  * data.nii              :the DTI data corrected for the eddy current effect.
  * nodif_brain_mask.nii  :binary mask of the brain   
  * bvals                 :text file containing b values in FSL format
  * bvecs                 :text file containing the gradient directions in FSL format (NOTE: directions must be normalized)
  
NOTE: the final mask could be a little different from the original one because FSL's tools need a strictly binary mask.

Example:
format_data_FSL /mind_cart/denis/50219 [/mind_cart/denis/50219/fsl_data/ /mind_cart/denis/bvals /mind_cart/denis/bvecs]

"""

# IMPORT
import os, sys

# DEFINING BUILT IN FUNCTIONS
def get_script_path(temp_path,file_sep): # Return the script path
	words = temp_path.split(file_sep)
	if len(words)==1:
		return ('')
	else:
		path = ''
		for cont in range(len(words)-1):
			path = path + words[cont] + file_sep
		return (path)
			
# SETTINGS
# Default settings
file_sep      = '/'
script_path   = get_script_path(sys.argv[0],file_sep) # the path where the script is stored
output_subdir = 'fsl_data/'
input_dir     = 'original_data/'
bvals_file_in = 'fsl_bvals'
bvecs_file_in = 'fsl_bvecs'
data_file     = 'DTIdata.nii.gz'
mask_file     = 'DTI_mask.nii.gz'

bvals_file_out = 'bvals'
bvecs_file_out = 'bvecs'
data_file_out  = 'data.nii.gz'
mask_file_out  = 'nodif_brain_mask.nii.gz'

eddy_correction_options = '0' # Volume to use as reference for the eddy correction

# STEP 1: reading input data
if len(sys.argv)<2: # NOTE: the first sys.argv element is the script name
	print 'WARNING: not enough input argument. You must specify at least the subject directory.'
	raise 'Input data error'
else:
	sub_dir = sys.argv[1];
	
	# Making sure that sub_dir finishes with /
	if sub_dir[len(sub_dir)-1]!=file_sep:
		sub_dir = sub_dir + file_sep
	
if len(sys.argv)>=3: # A new output directory is given
	output_dir = sys.argv[2]
else:
	output_dir = sub_dir + output_subdir
if output_dir[len(sub_dir)-1]!=file_sep:
	output_dir = output_dir + file_sep

if len(sys.argv)>=3: # bvals file
	bvals_file = sys.argv[3]
else:
	bvals_file = script_path + bvals_file_in
	
if len(sys.argv)>=4: # bvals file
	bvecs_file = sys.argv[4]
else:
	bvecs_file = script_path + bvecs_file_in

	
# Option recap
print 'FORMATTING DATA FOR FSL TOOLS'
print ' subject path: %s' %(sub_dir)
print ' output path:  %s' %(output_dir)
print ' bvals file:   %s' %(bvals_file)
print ' bvecs file:   %s' %(bvecs_file)
print ' '

# STEP 2: checking for the presence of the data
print 'Checking data...'
error_flag = 0
if not(os.path.isfile(bvals_file)):
	print 'WARNING: file %s is not present' %( bvals_file)
	error_flag = 1
	
if not(os.path.isfile(bvecs_file)):
	print 'WARNING: file %s is not present' %( bvecs_file)
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

# STEP 3: preparing output directory and files
print 'Copying files...'
cmd = 'mkdir ' + output_dir
os.system(cmd)
cmd = 'chmod 777 ' + output_dir
os.system(cmd)

cmd = 'fslmaths ' + sub_dir + input_dir + mask_file + ' -mul -1 -add 1 -thr 0.5 -add -1 -mul -1 -thr 0.5 ' + output_dir + mask_file_out # This command make sure the mask is really binary!
os.system(cmd)
cmd = 'chmod 666 ' + output_dir + mask_file_out
os.system(cmd)

cmd = 'cp ' + bvals_file + ' ' + output_dir + bvals_file_out
os.system(cmd)
cmd = 'chmod 666 ' + output_dir + bvals_file_out
os.system(cmd)

cmd = 'cp ' + bvecs_file + ' ' + output_dir + bvecs_file_out
os.system(cmd)
cmd = 'chmod 666 ' + output_dir + bvecs_file_out
os.system(cmd)

# STEP 4: Eddy correction of the data
print 'Performing eddy correction...'
cmd = 'eddy_correct ' + sub_dir + input_dir + data_file + ' ' + output_dir + data_file_out + ' ' + eddy_correction_options
os.system(cmd)
cmd = 'chmod 666 ' + output_dir + data_file_out 
os.system(cmd)
cmd = 'chmod 666 ' + output_dir + '*.*' 
os.system(cmd)