"""
DTI tool - full_preparation_FSL.py
(Denis Peruzzo - 9/24/2010)

The script prepares the T1 and DTI data for future elaboration.

INPUT: 
- Subject directory (e.g. /mind_cart/denis/50241/). It must already contain:

CALL AND PERFORM
- prepare_data          :to coregister T1 and DTI images, move T1 image to DTI space and mask DTI images
- format_data_FSL       :prepare the data in FSL format
- fa_map_FSL            :compute the tensor FIT and the FA map
- preparation4tract_FSL :call bedpostX
"""

# IMPORT
import os, sys

file_sep = '/'

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

# STEP 1: reading input data
script_path   = get_script_path(sys.argv[0],file_sep) # the path where the script is stored

if len(sys.argv)<2: # NOTE: the first sys.argv element is the script name
	print 'WARNING: not enough input argument. You must specify at least the subject directory.'
	raise 'Input data error'
else:
	sub_dir = sys.argv[1];
	
# STEP 2: prepare_data
cmd = 'python ' + script_path + 'prepare_data.py ' + sub_dir
os.system(cmd)

# Step 3: format_data_FSL
cmd = 'python ' + script_path + 'format_data_FSL.py ' + sub_dir
os.system(cmd)

# Step 4: fa_map_FSL
cmd = 'python ' + script_path + 'fa_map_FSL.py ' + sub_dir
os.system(cmd)

# Step 5: prepare4tract_FSL
#cmd = 'python ' + script_path + 'preparation4tract_FSL.py ' + sub_dir
#os.system(cmd)

print 'SUBJECT %s COMPLETE!' %(sub_dir)