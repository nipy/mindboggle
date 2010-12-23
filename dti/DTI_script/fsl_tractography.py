"""
DTI tool - fsl_tractography.py
(Denis Peruzzo - 9/29/2010)

The script prepares the seed and perform the tractography on the data provided in input

INPUT: 
 * sub_dir
 * seed_file
 * 
"""

# IMPORT
import os, sys
import numpy as np
import nibabel as nb
from seed_lib import generate_seed_masks
from seed_lib import move_seed2dti_space
from fsl_lib import do_track	
			
# SETTINGS
# Default settings
file_sep           = os.sep
output_subdir      = 'fsl_data.tracts'+file_sep
data_dir           = 'original_data'+file_sep
T12DTItransf       = 'T12DTI_transformation.mat'
allseeds_file      = 'all_seeds_mask.nii.gz'
components_file    = 'labelled_seeds.nii.gz'
seed_mask_basename = 'seed_'
DTImask_file       = 'nodif_brain_mask.nii.gz' # in bedpostX folder
bedpostX_dir       = 'grant_fsl_data.bedpostX'+file_sep
output_seed_dir    = 'seeds'+file_sep
comp_mask          = 'components_mask.nii.gz'
output_tract_dir   = 'tract_'
tract_map_file     = 'fdt_paths.nii.gz'

# STEP 1: reading input data
if len(sys.argv)<3: # NOTE: the first sys.argv element is the script name
	print 'WARNING: not enough input argument. You must specify at least the subject directory and the seed mask'
	raise 'Input data error'
else:
	sub_dir = sys.argv[1]
	
	
	# Making sure that sub_dir finishes with /
	if sub_dir[len(sub_dir)-1]!=file_sep:
		sub_dir = sub_dir + file_sep
	
	seed_file=sys.argv[2]
	
	
if len(sys.argv)>=4: # A new output directory is given
	output_dir = sys.argv[3]
else:
	output_dir = sub_dir + output_subdir
if output_dir[len(sub_dir)-1]!=file_sep:
	output_dir = output_dir + file_sep
	
# Option recap
print 'DTI DATA PRELIMINARY STEPS'
print ' subject path: %s' %(sub_dir)
print ' output path:  %s' %(output_dir)
print ' seed mask:    %s' %(seed_file)

# STEP 2: checking for the presence of the data
print 'Checking data...'
error_flag = 0
T1data = 0
	
if not(os.path.isfile(sub_dir + data_dir + T12DTItransf)):
	print 'WARNING: file %s is not present in %s' %(T12DTItransf, sub_dir + data_dir)
	error_flag = 1
	
if not(os.path.isfile(seed_file)):
	print 'WARNING: file %s not found' %(seed_file)
	error_flag = 1	
	
if not(os.path.exists(sub_dir + bedpostX_dir )):
	print 'WARNING: path %s does not exist in %s' %(bedpostX_dir, sub_dir)
	error_flag = 1		

# STEP 3: preparing the seed masks
print 'preparing the seed masks...'
if not(os.path.exists(output_dir+output_seed_dir)):
	move_seed2dti_space(seed_file,sub_dir+bedpostX_dir+DTImask_file,sub_dir+data_dir+T12DTItransf,output_dir)
	generate_seed_masks(output_dir+allseeds_file,output_dir+output_seed_dir)
else:
	print 'Seed masks already found... the step is skipped'

# STEP 4: preforming tractography
print 'Tractography step'

# Reading the component mask to compute the number of components
fileObj = nb.load(output_dir + output_seed_dir+comp_mask)
data = fileObj.get_data()
nComp = np.max(data)
del data
del fileObj

print '%s seed masks found' %(str(nComp))

for cont in range(nComp):	
	comp=cont+1
	print 'Performing tractography on seed %s' %(str(comp))
	comp_seed=output_dir+output_seed_dir+'seed_'+str(comp)+'_mask.nii.gz'
	comp_output_dir=output_dir+output_tract_dir+str(comp)
	print 'Looking for output file %s ' %(comp_output_dir+file_sep+tract_map_file)
	if os.path.isfile(comp_output_dir+file_sep+tract_map_file):
		print '   output file %s already found in %s: analysis skipped' %(tract_map_file,output_dir+output_tract_dir+str(comp))
	else:
		print '   file not found'
		do_track(comp_seed,sub_dir+bedpostX_dir,comp_output_dir)