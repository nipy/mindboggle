"""
DTI tool - inport_data.py
(Denis Peruzzo - 09/24/2010)

The script imports the T1 and DTI data for future elaboration.

INPUT: 
  * Destination folder: the folder where the subject elaboration will be stored. 
			  Usually it ends with the subject ID. (e.g. /mind_cart/denis/50201/)
  * T1 data:  the file containing the T1 weighted images.
		      NOTE: it is expected to be a cropped brain. (e.g. /data/BI/mindburst/data_pieces/mri_atropos_crop/1/64/20/)
  * T1 mask:  the file containing the T1 mask. It is usually made by using Atropos.
		      (e.g. /data/BI/mindburst/data_pieces/mri_mask/1/94/73/)
  * DTI data: The raw DTI data file. (e.g. /data/BI/human/MRIBANK/50246/raw/s006a001_TENSOR_ASSET_25_Dir.nii.gz)

OUTPUT
In the specified destination a subfolder called 'original_data' is made. 
If such folder already exist a warning message is given and the operator is asked to decide what to do.
The 'original_data' subfolder will contain: 
  * T1data_masked.nii.gz : the T1 data.
  * T1_mask.nii.gz       : the T1 mask.
  * DTIdata.nii.gz       : the DTI data.
  
Example:
import_data /mind_cart/denis/50201/ /data/BI/mindburst/data_pieces/mri_atropos_crop/1/64/20/ /data/BI/mindburst/data_pieces/mri_mask/1/94/73/ /data/BI/human/MRIBANK/50246/raw/s006a001_TENSOR_ASSET_25_Dir.nii.gz

"""

# IMPORT
import os, sys
			
# SETTINGS
# Default settings
existing_data_warning = 1 # If the destination subfolder already exist 

file_sep      = '/'
output_subdir = 'original_data/'
DTIdata_file  = 'DTIdata.nii.gz'
T1data_file   = 'T1data_masked.nii.gz'
T1mask_file   = 'T1_mask.nii.gz'
report_file   = 'import_file_report.txt'

# STEP 1: reading input data
if len(sys.argv)<5: # NOTE: the first sys.argv element is the script name
	print 'WARNING: not enough input argument.'
	print 'You must specify at least the following parameters:'
	print ' - destination folder'
	print ' - T1 data file'
	print ' - T1 mask file'
	print ' - DTI data file'
	raise 'Input data error'
else:
	destination_dir = sys.argv[1];
	
	# Making sure that sub_dir finishes with /
	if destination_dir[len(destination_dir)-1]!=file_sep:
		destination_dir = destination_dir + file_sep

	T1data_file_IN  = sys.argv[2];
	T1mask_file_IN  = sys.argv[3];
	DTIdata_file_IN = sys.argv[4];
	
# Option recap
print 'IMPORT DATA'
print ' output path: %s' %(destination_dir)
print ' T1 data:   %s' %(T1data_file_IN)
print ' T1 mask:   %s' %(T1mask_file_IN)
print ' DTI data:  %s' %(DTIdata_file_IN)
print ' '

# STEP 2: checking for the presence of the data
print 'Checking data...'
error_flag = 0

	
if not(os.path.isfile(T1data_file_IN)):
	print 'WARNING: file %s does not exist' %(T1data_file_IN)
	error_flag = 1	
	
if not(os.path.isfile(T1mask_file_IN)):
	print 'WARNING: file %s does not exist' %(T1mask_file_IN)
	error_flag = 1	
	
if not(os.path.isfile(DTIdata_file_IN)):
	print 'WARNING: file %s does not exist' %(DTIdata_file_IN)
	error_flag = 1	

if os.path.exists(destination_dir + output_subdir) and existing_data_warning:
	print 'WARNING: destination folder %s already exist!' %(destination_dir + output_subdir)
	answer = raw_input('Overwrite data (Y/N)? ')
	if (answer != 'Y') and (answer != 'y'):
		error_flag = 1
	
if error_flag:
	print 'Cannot perform the script. Please check the input data.'
	raise 'Input data error'	
	
# STEP 3: importing data
print 'Importing data...'

if not(os.path.exists(destination_dir)):
	cmd = 'mkdir -p -m 777 ' + destination_dir
	os.system(cmd)
	
cmd = 'mkdir -p -m 777 ' + destination_dir + output_subdir
os.system(cmd)

cmd = 'cp -f ' + T1data_file_IN + ' ' + destination_dir + output_subdir + T1data_file
os.system(cmd)

cmd = 'chmod 666 ' + destination_dir + output_subdir + T1data_file
os.system(cmd)

cmd = 'cp ' + T1mask_file_IN + ' ' + destination_dir + output_subdir + T1mask_file
os.system(cmd)

cmd = 'chmod 666 ' + destination_dir + output_subdir + T1mask_file
os.system(cmd)

cmd = 'cp ' + DTIdata_file_IN + ' ' + destination_dir + output_subdir + DTIdata_file
os.system(cmd)

cmd = 'chmod 666 ' + destination_dir + output_subdir + DTIdata_file
os.system(cmd)

# STEP 4: making the report file
fobj = open(destination_dir + output_subdir + report_file,'w')
fobj.write('IMPORT DATA')
fobj.write(' output path: %s' %(destination_dir))
fobj.write(' T1 data:   %s' %(T1data_file_IN))
fobj.write(' T1 mask:   %s' %(T1mask_file_IN))
fobj.write(' DTI data:  %s' %(DTIdata_file_IN))
fobj.close()