"""
DTI tool - prepare_data_old.py
(Denis Peruzzo - 9/24/2010)

The script prepares the T1 and DTI data for future elaboration.

INPUT: 
- Subject directory (e.g. /mind_cart/denis/50241/). It must already contain:
  * sub_dir/original_data/DTIdata.nii       :the DTI data
- Output directory (optional. Default is sub_dir/original_data/)

OUTPUT
The output directory will contain:
  * DTI_mask.nii                 :the DTI data corrected for the eddy current effect.
  * DTI_firstVolume.nii.gz       :brain extracted from the first volume of the DTI data
  
NOTE:
- Fsl's flirt tool performs better if the two images have already been betted.
- Flirt provides a better result if the lower resolution image is coregistered to the higher resolution image.
- Fsl's tool convert_xfm can be used to invert a space transformation provided by flirt.
  
"""

# IMPORT
import os, sys
			
# SETTINGS
# Default settings
file_sep      = os.sep
output_subdir = 'original_data' + file_sep
input_dir     = 'original_data' + file_sep
DTIdata_file  = 'DTIdata.nii.gz'
T1data_file   = 'T1data_masked.nii.gz'
T1mask_file   = 'T1_mask.nii.gz'

DTI2T1transf    = 'DTI2T1_transformation.mat'
T12DTItransf    = 'T12DTI_transformation.mat'
T1data2DTI      = 'T1data_2DTIspace.nii.gz'
DTIdata2T1      = 'DTIdata_2T1space.nii.gz'
DTImask_file    = 'DTI_mask.nii.gz'
DTIvol_file     = 'DTI_firstVolume' # NOTE: the file name must be without the extension
ANTS_T_files    = 'DTI2T1_ants_T_'

bet_options = '-f 0.2 -m' # The threshold is reduced to have larger brain mask (-f option) and the mask is also saved (-m option)

ANTS_T_quantification = 'ANTS 3 -m PR[%s,%s,1,2] -r Gauss[2,0] -o %s -i 100x100x100 -t SyN[0,5] --number-of-affine-iterations 1000x1000x1000x1000x1000'
ANTS_invT_application = 'WarpImageMultiTransform 3 %s %s -R %s -i %sAffine.txt %sInverseWarp.nii.gz'

'ANTS 3 -m PR[DTI_firstVolume.nii.gz,T1data_masked.nii.gz,1,2] -r Gauss[2,0,0.1] -o ants6_ -i 100x100x100 -t SyN[1,2,0.05] --number-of-affine-iterations 1000x1000x1000x1000x1000'


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
	
# Option recap
print 'DTI DATA PRELIMINARY STEPS'
print ' subject path: %s' %(sub_dir)
print ' output path:  %s' %(output_dir)
print ' '

# STEP 2: checking for the presence of the data
print 'Checking data...'
error_flag = 0
T1data = 0
	
if not(os.path.isfile(sub_dir + input_dir + DTIdata_file)):
	print 'WARNING: file %s is not present in %s' %(DTIdata_file, sub_dir + input_dir)
	error_flag = 1
	
if not(os.path.isfile(sub_dir + input_dir + T1data_file)):
	print 'T1 data %s found in %s.' %(T1data_file, sub_dir + input_dir)
	error_flag = 1
	
if not(os.path.isfile(sub_dir + input_dir + T1mask_file)):
	print 'T1 data %s found in %s.' %(T1mask_file, sub_dir + input_dir)
	error_flag = 1
	
if error_flag:
	print 'Cannot perform the script. Please check the input data.'
	raise 'Input data error'
	
# STEP 3: assessment of the transformation from T1 to DTI space
print 'Extracting and cropping the first DTI volume...'

# DTI first volume and brain extraction
tmp_dir = sub_dir + input_dir + 'tmp'
cont = 1
while os.path.isfile(tmp_dir + str(cont)):
	cont = cont +1
tmp_dir = tmp_dir + str(cont) + file_sep
os.system('mkdir ' + tmp_dir)
os.system('fslsplit ' + sub_dir + input_dir + DTIdata_file + ' ' + tmp_dir + 'vol -t')
os.system('bet ' + tmp_dir + 'vol0000.nii.gz ' + sub_dir + input_dir + DTIvol_file + '.nii.gz ' + bet_options)
os.system('chmod 666 ' + sub_dir + input_dir + DTIvol_file + '.nii.gz')
os.system('mv ' + sub_dir + input_dir + DTIvol_file + '_mask.nii.gz ' + sub_dir + input_dir + DTImask_file)
os.system('chmod 666 ' + sub_dir + input_dir + DTImask_file)
os.system('rm -R ' + tmp_dir)

# STEP 4: assessment of the T1 to DTI transformation
# NOTE
# ANTS is used to assess the transformation
# - the transformation is computed on betted images to exclude scalp and skull interference in the coregistration process.
# - the lower resolution image is moved to the higher resolution image to get better transformation estimate.
# - T1 image is moved to DTI space using the inverse transformation.

print 'Computing the T1 -> DTI transformation...'

# Computing the DTI->T1 transformation
cmd = ANTS_T_quantification %(sub_dir+input_dir+T1data_file,sub_dir+input_dir+DTIvol_file+'.nii.gz',sub_dir+input_dir+ANTS_T_files)
os.system(cmd)
cmd = 'chmod 666 ' + sub_dir + input_dir + ANTS_T_files + '*.*'
os.system(cmd)

# Moving the T1 data to DTI space
cmd = ANTS_invT_application %(sub_dir+input_dir+T1data_file, sub_dir+input_dir+T1data2DTI, sub_dir+input_dir+DTIvol_file+'.nii.gz', sub_dir+input_dir+ANTS_T_files,sub_dir+input_dir+ANTS_T_files) 
os.system(cmd)
cmd = 'chmod 666 ' + sub_dir + input_dir + T1data2DTI
os.system(cmd)