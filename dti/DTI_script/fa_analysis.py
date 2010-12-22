"""
FA ANALYSIS
(made by Denis Peruzzo, October, 14th, 2010)
The method analyzes the FA maps (or a given map) on a ROI based approach.
It can import ROI maps from a file, convert it to niftii format, move to 
DTI space and superimpose the ROI map to the FA map to characterize the 
ROI. Finally the results are saved on a text file.

The input is a text file containing a subject in each row with the following layout:
subID label_map FA_map
"""

# IMPORT SECTION
import os,sys
import maps_analysis_lib as ma
import string

# SETUP VARIABLES
# Modify this section to modify the analysis
#input_file = '/mind_cart/DTI/FA_spreadsheets/ArnoAnalysis_subjects.txt'
#input_file = '/mind_cart/DTI/FA_spreadsheets/prova_subjects.txt'
input_file = '/mind_cart/DTI/FA_spreadsheets/Elizabeth_subjects.txt'
#input_file = '/mind_cart/DTI/FA_spreadsheets/KevinOschnerAnalysis_subjects_importWM.txt'
import_file_flag = False # True: import the ROI map files
import_folder = '/mind_cart/DTI/%s/ROImaps/WMparc/' # Where import the ROI maps
coregistration_method = 'ants'


# GLOBAL VARIABLES
comment_character = '#'
separator_character = ' '


# MAIN CODE
# Checking input file
if not(os.path.isfile(input_file)):
	print 'ERROR: %s file not found' %(input_file)
	raise

# Reading the imput file
ROI_map_in_list={}
value_map_list={}
mask_list={}
try:
	file_obj = open(input_file,'r')
except:
	print 'ERROR: can not open file %s' %(input_file)
for line in file_obj:
	if line[0] != comment_character:
		# The row contains a subject
		line = line.strip()
		elements = line.split(separator_character)
		if len(elements) != 4:
			# Error in line, the row doesn't contain the correct number of elements
			print 'ERROR: reading file %s' %(input_file)
			print '       can not interpretate row: %s' %(line)
		else:
			ROI_map_in_list[elements[0]]=elements[1]
			value_map_list[elements[0]]=elements[2]
			mask_list[elements[0]]=elements[3]

# Importing the ROI map files
if import_file_flag:
	# Files in ROI_map_list have to be imported
	print 'Importing ROI map files'
	ROI_map_list = {}
	
	for subID in ROI_map_in_list.keys():
		file_2_import = ROI_map_in_list[subID]
		target_folder = import_folder %(subID)
		
		file_location=ma.import_ROI_map(file_2_import,target_folder,subID,coregistration_method)
		ROI_map_list[subID]=file_location
else:
	ROI_map_list = ROI_map_in_list
	
# Analyzing the maps and making the spreadsheets
# MEAN
result_file_suffix = '_mean.txt'
measure = 'mean'
output_file=os.path.dirname(input_file)+os.sep+os.path.splitext(os.path.basename(input_file))[0]+result_file_suffix
ma.ROI_analysis_list(ROI_map_list,value_map_list,mask_list,output_file,measure)

# SD
result_file_suffix = '_sd.txt'
measure = 'sd'
output_file=os.path.dirname(input_file)+os.sep+os.path.splitext(os.path.basename(input_file))[0]+result_file_suffix
ma.ROI_analysis_list(ROI_map_list,value_map_list,mask_list,output_file,measure)

# N# voxel
result_file_suffix = '_nVoxel.txt'
measure = 'NVOXEL'
output_file=os.path.dirname(input_file)+os.sep+os.path.splitext(os.path.basename(input_file))[0]+result_file_suffix
ma.ROI_analysis_list(ROI_map_list,value_map_list,mask_list,output_file,measure)
