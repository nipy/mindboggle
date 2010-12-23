"""
MAPS_ANALYSIS_LIB

list of python functions to analyze maps (made for FA maps and ROI analysis)
  
"""	



# IMPORT
import os, sys
import string
import numpy as np
import nibabel as nb

# GLOBAL VARIABLE
file_sep = os.sep


ants_inv_warp_template = os.sep+'mind_cart'+os.sep+'DTI'+os.sep+'%s'+os.sep+'original_data'+os.sep+'DTI2T1_ants_T_InverseWarp.nii.gz'
ants_affine_template = os.sep+'mind_cart'+os.sep+'DTI'+os.sep+'%s'+os.sep+'original_data'+os.sep+'DTI2T1_ants_T_Affine.txt'
ants_folder_template = os.sep+'mind_cart'+os.sep+'DTI'+os.sep+'%s'+os.sep+'original_data'+os.sep
ants_invT_application = 'WarpImageMultiTransform 3 %s %s -R %s -i %s %s --use-NN' # %(inFile,outFile,refFile,Affine,InverseWarp)

mri_convert_command = 'mri_convert %s %s -rl %s -ns 1 -nc -rt nearest'
T1data_template = os.sep+'mind_cart'+os.sep+'DTI'+os.sep+'%s'+os.sep+'original_data'+os.sep+'T1data_masked.nii.gz'
DTIdata_template = os.sep+'mind_cart'+os.sep+'DTI'+os.sep+'%s'+os.sep+'original_data'+os.sep+'DTI_firstVolume.nii.gz'

default_output_text_file = 'maps_analysis_output.txt'
new_line_separator = os.linesep
field_text_file_reparator = ' '

freesurfer_label_file = 'mri___aparc+aseg.mgz' # Parcellation based on Desikan-Killiany Atlas
#freesurfer_label_file = 'mri___aparc.a2009s+aseg.mgz' # Parcellation based on Destrieux Atlas
#freesurfer_label_file = 'mri___aseg.mgz' # only raw segmentation

label_subfolder = 'labels' + file_sep
label_import_name = 'freesurfer_aparc+aseg'

label_file_mind_burst = 'mri___aparc+aseg.mgz'
label_file_freeS = 'original_data' + file_sep + 'aparc_aseg_labels_freeSurfer.mgz'
label_file_T1 = 'original_data' + file_sep + 'aparc_aseg_labels_T1.nii.gz'
label_file_DTI = 'original_data' + file_sep + 'aparc_aseg_labels_DTI.nii.gz'

fslT12DTI_transformation_matrix = 'original_data' + file_sep + 'T12DTI_transformation.mat'

fa_map_path = file_sep + 'fsl_data.fa_map' + file_sep
fa_map_file = '_FA.nii.gz'
mask_file = file_sep + 'fsl_data' + file_sep + 'nodif_brain_mask.nii.gz'
DTI_first_volume_file = file_sep + 'original_data' + file_sep + 'DTI_firstVolume.nii.gz'
T1_data_file = file_sep + 'original_data' + file_sep + 'T1data_masked.nii.gz'
main_path = file_sep + 'mind_cart' + file_sep + 'DTI' + file_sep

# METHODS

# --------------------------------------------------------------------------
# FILE_LOCATION_FROM_FILE(file_name,target_folder,coregistration_method)
# Import the given file to the target folder, transform it in a niftii file
# and move it to DTI space
def import_ROI_map(file_name,target_folder,sub_id,coregistration_method):
	# Checking the input data
	if not(os.path.isfile(file_name)):
		print 'ERROR: file %s not found' %(file_name)
	if not(os.path.exists(target_folder)):
		os.system('mkdir -m 777 '+target_folder)
	if target_folder[-1] != os.sep:
		target_folder = target_folder + os.sep
	
	# Preparing data
	file_basename = os.path.basename(file_name)
	file_basename_elements = os.path.splitext(file_basename)
	
	target_file_T1_name = target_folder+file_basename_elements[0] + '_T1.nii.gz'
	target_file_DTI_name = target_folder+file_basename_elements[0] + '_DTI.nii.gz'
	T1data = T1data_template %(sub_id)
	DTIdata= DTIdata_template %(sub_id)
	ants_inv_warp = ants_inv_warp_template %(sub_id)
	ants_affine = ants_affine_template %(sub_id)

	
	# Copying the file to the target location
	cmd = 'cp %s %s' %(file_name,target_folder+file_basename)
	os.system(cmd)
	cmd = 'chmod 666 %s' %(target_folder+file_basename)
	os.system(cmd)
	
	# Moving the file to DTI space		
	# - From original format to niftii
	cmd = mri_convert_command %(target_folder+file_basename,target_file_T1_name,T1data)
	os.system(cmd)
	cmd = 'chmod 666 %s' %(target_file_T1_name)
	os.system(cmd)

	# - From T1 to DTI space
	if coregistration_method == 'ants':
		cmd = ants_invT_application %(target_file_T1_name,target_file_DTI_name,DTIdata,ants_affine,ants_inv_warp)
		# %(inFile,outFile,refFile,Affine,InverseWarp)
		os.system(cmd)
		cmd = 'chmod 666 %s' %(target_file_DTI_name)
		os.system(cmd)
	else:
		print 'ERROR: coregistration method unknown (%s)' %(coregistration_method)
	

	# Preparing the output result
	result_location = target_file_DTI_name
	return result_location
	
# --------------------------------------------------------------------------
# READ_SUBJECTS_FROM_FILE(file_name)
# Read a text file containing the ID of the subjects and the location of the
# freesurfer label files.
# Return a list of lists.
def read_subjects_from_file(file_name):
	final_list=[]
	if os.path.isfile(file_name):
		try:
			fileOBJ = open(file_name,'r')
		except:
			print '\n\nERROR: can not read file %s\n\n' %(file_name)
			raise
		else:			
			for line in fileOBJ:
				line = line.rstrip()
				fields = line.split(' ')
				final_list.append([fields[0],fields[1]])
			fileOBJ.close()
	else:
		print '\n\nERROR: file %s not found\n\n' %(file_name)
		raise
	return final_list

# --------------------------------------------------------------------------
# GET_MAP_FROM_FILE(map_file)
# Read the fa/connectivity/probability/mask map from a file.
# INPUT
#   * map_file :map file path
# OUTPUT
#   * map
def get_map_from_file(map_file):
	map = 0
	if os.path.isfile(map_file):
		# Reading the file

		try:
			map_file_obj = nb.load(map_file)
			map = map_file_obj.get_data()
		except:
			print 'Unexpected error while reading from %s' %(map_file)
	else:
		print '\n\nERROR: file %s not found\n\n' %(map_file)
	return map

# --------------------------------------------------------------------------
# IMPORT_FREESURFER_LABEL_FILE(import_file,destination_file)
# Import the file given by freesurfer and move it to DTI space
def import_freesurfer_label_file(import_file,destination_file):
	if os.path.isfile(import_file):
		# Coping the file
		cmd = 'cp ' + import_file + ' ' + destination_file
		os.system(cmd)
		cmd = 'chmod 666 ' + destination_file
		os.system(cmd)
	else:
		print '\n\nERROR: file %s not found.\n\n' %(import_file)
		raise

# --------------------------------------------------------------------------
# MOVE_LABELS_FREES2DTI(freeS_labels,T1_labels,DTI_labels,T1_data,DTI_data,T12DTImat)
# Move the freeSurfer labels to T1 space by using mri_convert tool and then
# to DTI space using the fsl flirt tool.
def move_labels_frees2dti(freeS_labels,T1_labels,DTI_labels,T1_data,DTI_data,T12DTImat):
	if not(os.path.isfile(freeS_labels)):
		print '\n\nERROR: file %s not found.\n\n' %(freeS_labels)
		raise
	if not(os.path.isfile(T1_data)):
		print '\n\nERROR: file %s not found.\n\n' %(T1_data)
		raise
	if not(os.path.isfile(DTI_data)):
		print '\n\nERROR: file %s not found.\n\n' %(DTI_data)
		raise
	if not(os.path.isfile(T12DTImat)):
		print '\n\nERROR: file %s not found.\n\n' %(T12DTImat)
		raise

	# From freeSurfer to T1
	cmd = 'mri_convert ' +  freeS_labels + ' ' + T1_labels + ' -rl ' + T1_data + ' -ns 1 -nc -rt nearest'
	os.system(cmd)
	cmd = 'chmod 666 ' + T1_labels
	os.system(cmd)
	# From T1 to DTI (FSL)
	#cmd = 'flirt -in ' + T1_labels + ' -ref ' + DTI_data + ' -out ' + DTI_labels + ' -init ' + T12DTImat + ' -applyxfm -interp nearestneighbour'
	#os.system(cmd)
	#cmd = 'chmod 666 ' + DTI_labels
	#os.system(cmd)
	

# --------------------------------------------------------------------------
# IMPORT_FREESURFER_FILES(file_list)
# Import the freeSurfer label files from a list and move them to DTI space
def import_freesurfer_file(file_list):
	for element in file_list:
		sub_ID       = element[0]
		file_name    = element[1]+ label_file_mind_burst
		freeS_labels = main_path + str(sub_ID) + label_file_freeS
		T1_labels    = main_path + str(sub_ID) + label_file_T1
		DTI_labels   = main_path + str(sub_ID) + label_file_DTI
		T1_data      = main_path + str(sub_ID) + T1_data_file
		DTI_data     = main_path + str(sub_ID) + DTI_first_volume_file
		T12DTImat    = main_path + str(sub_ID) + T12DTI_transformation_matrix
		try:
			import_freesurfer_label_file(file_name,freeS_labels)
			move_labels_frees2dti(freeS_labels,T1_labels,DTI_labels,T1_data,DTI_data,T12DTImat)
		except:
			print '\n\nERROR: subject %s not imported correctly.\n\n' %(str(sub_ID))

# --------------------------------------------------------------------------
# SINGLE_ROI_ANALYSIS(map,roi)
# The method analyze the map values in the given single ROI.
# It computes:
# - ROI dimension (as number of voxels)
# - mean value in the ROI
# - standard deviation value in the ROI
# - max value in the ROI
# - min value in the ROI
#
# OUTPUT: 
# [voxel_number, mean_values, std_value, max_value, min_value]

def single_roi_analysis(map,roi,measure_method):
	if map.shape == roi.shape:
		# Dimensions are ok
		roi_indices = roi.nonzero()
		if len(roi_indices[0])>0:
			vett_values = map[roi_indices]
		
			voxel_number = len(vett_values)
			mean_value = vett_values.mean()
			max_value = vett_values.max()
			min_value = vett_values.min()
			std_value = vett_values.std()
		else:
			voxel_number = 0
			mean_value = 0
			max_value = 0
			min_value = 0
			std_value = 0
	else:
		print 'ERROR: map and roi do not have the same dimension'
		raise
		voxel_number = 0
		mean_value = 0
		max_value = 0
		min_value = 0
		std_value = 0
	if measure_method.upper() == 'MEAN':
		return [mean_value]
	elif measure_method.upper() == 'MAX':
		return [max_value]
	elif measure_method.upper() == 'MIN':
		return [min_value]
	elif (measure_method.upper() == 'STD')or(measure_method.upper() == 'SD'):
		return [std_value]
	elif measure_method.upper() == 'NVOXEL':
		return [voxel_number]
	elif measure_method.upper() == 'ALL':
		return [voxel_number, mean_value, std_value, max_value, min_value]
	else:
		print 'ERROR: measurement option unknown (%s)' %(measure_method)
		raise
	
# --------------------------------------------------------------------------
# THREE_D_ROI_ANALYSIS(map,roi)
# The method analyze the map values in the given single ROI.
# It computes:
# - ROI dimension (as number of voxels)
# - mean value in the ROI
# - standard deviation value in the ROI
# - max value in the ROI
# - min value in the ROI
#
# OUTPUT: a list. Each element of the list represents a ROI and it is a list
# containing: [voxel_number, mean_values, std_value, max_value, min_value]

def three_d_roi_analysis(map,roi,measure_method):
	if map.shape == roi.shape:
		# Dimensions are ok
		n_roi = roi.max()
		result_list = []
		for cont in range(n_roi):
			roi_number = cont+1
			roi_mask = roi == roi_number
			roi_results = single_roi_analysis(map,roi_mask,measure_method)
			result_list= result_list+roi_results
	else:
		# map and roi dimensions don't agree
		print 'ERROR: map and roi do not have the same dimension'
		raise
	return result_list



# --------------------------------------------------------------------------
# MAKE_TEXT_FILE(data_list,file_name)
# Create a text file with the given name and fill it with the content of
# data_list.
# data_list is a list of list. Each insider list contains the data which 
# should be saved in the same row of the text file.

def make_text_file(data_list,file_name):
	
	text_list = []
	for row_array in data_list:
		# Each list in data_list contains the data for a subject
		row_string = ''
		for field in row_array:
			row_string = row_string + str(field)
			row_string = row_string + field_text_file_reparator
		text_list.append(row_string)
	
	fileObj = open(file_name,'w')
	fileObj.writelines(['%s%s' %(l,new_line_separator) for l in text_list])
	fileObj.close()

# --------------------------------------------------------------------------
# SUBJECT_ANALYSIS(sub_ID)
def subject_analysis(map_file,roi_file,mask_file,measure_method):
	fa_map = get_map_from_file(map_file)
	roi_map = get_map_from_file(roi_file)
	mask_map = get_map_from_file(mask_file)
	
	if mask_map.shape == roi_map.shape:
		mask_indices = mask_map.nonzero()
		final_roi = np.zeros(roi_map.shape,dtype=roi_map.dtype)
		final_roi[mask_indices] = roi_map[mask_indices]
		sub_results = three_d_roi_analysis(fa_map,final_roi,measure_method)
	else:
		print '\n\nERROR: mask and roi label dimensions do not agree\n\n'
		sub_results = 0
	return sub_results

# --------------------------------------------------------------------------
# SUBJECT_ANALYSIS(sub_ID)
def	numberlist2strlist(number_list):
	str_list=''
	for element in number_list:
		str_list = str_list + str(element)
		str_list = str_list + field_text_file_reparator	
	return str_list
	
# --------------------------------------------------------------------------
# FA_ROI_ANALYSIS_SUBJECTS(list,output_file)
# Apply the ROI based analysis to FA map for subjects in the list.
# Results are save in a text file.
def fa_roi_analysis_subjects(list,output_file):
	results = []
	print 'Subject:'
	for sub_ID in list:
		print ' - %s' %(str(sub_ID))
		#try:
		sub_results = subject_analysis(sub_ID)
		#print '%s : %s' %(str(sub_ID),str(sub_results))
		results.append([sub_ID + ' ' + numberlist2strlist(sub_results)])
		#except:
		#	print '\n\nERRORE: subject %s not correctly analyzed' %(str(sub_ID))
	print 'Saving the results...'
	make_text_file(results,output_file)
	
# --------------------------------------------------------------------------
# INSERT_ELEMENTS(parent_matrix,new_elements,row_ind,col_ind)
# Insert new_elements in parent_matrix starting from the cell [row_ind,col_ind]
def insert_elements(parent_matrix,new_elements,row_ind,col_ind):
	n_elem = len(new_elements)
	if (parent_matrix.shape[1]<n_elem+col_ind)or(row_ind+1>parent_matrix.shape[0]):
		# parent_matrix must be larger
		n_row=np.max([parent_matrix.shape[0],row_ind+1])
		n_col=np.max([parent_matrix.shape[1],col_ind+len(new_elements)])
		new_matrix = np.zeros([n_row,n_col])
		for r in range(parent_matrix.shape[0]):
			for c in range(parent_matrix.shape[1]):
				new_matrix[r][c] = parent_matrix[r][c]
	else:
		new_matrix = parent_matrix
	for e in range(n_elem):
		new_matrix[row_ind][col_ind+e]=new_elements[e]
	return new_matrix

# --------------------------------------------------------------------------
# INSERT_COLUMN_RIGHT
def insert_column_right(table,column):
	n_row,n_col= table.shape
	new_table = np.zeros([n_row,n_col+1])
	for r in range(n_row):
		new_table[r,n_col]=column[r]
		for c in range(n_col):
			new_table[r][c] = table[r][c]
	return new_table

# --------------------------------------------------------------------------
# PREPARE_RESULTS_4_TEXT_FILE(table,sub_list)
def prepare_results_4_text_file(table,sub_list):
	# Reduction of table dimensions

	reduced_table=np.zeros([table.shape[0],1])
	for c in range(table.shape[1]):
		vett = table[1:][:,c]
		if vett.max()>0:
			vett = table[:][:,c]
			reduced_table=insert_column_right(reduced_table,vett)
	
	final_table = []
	for r in range(reduced_table.shape[0]):
		sub_ID = 'Table'
		if r==0:
			sud_ID = 'Table'
		else:
			sub_ID = sub_list[r-1]
		row = reduced_table[r][:]
		final_table.append([sub_ID + field_text_file_reparator + numberlist2strlist(row)])
		
	return final_table
	
	
# --------------------------------------------------------------------------
# FA_ROI_ANALYSIS_SUBJECTS(list,output_file)
# Apply the ROI based analysis to FA map for subjects in the list.
# Results are save in a text file.
def ROI_analysis_list(ROI_list,MAP_list,mask_list,output_file,measure_method):
	results = []
	sub_list = ROI_list.keys()
	result_table=np.zeros([len(sub_list)+1,2])
	cont=1
	print 'ROI analysis'
	for sub_ID in sub_list:
		print ' - Subject %s' %(str(sub_ID))
		print '   ROI %s' %(ROI_list[sub_ID])
		print '   map %s' %(MAP_list[sub_ID])
		print '   Mask %s' %(mask_list[sub_ID])
		map_file = MAP_list[sub_ID]
		ROI_file = ROI_list[sub_ID]
		mask_file = mask_list[sub_ID]
		#try:
		sub_results = subject_analysis(map_file,ROI_file,mask_file,measure_method)
		#print '%s : %s' %(str(sub_ID),str(sub_results))
		#results.append([sub_ID + field_text_file_reparator + numberlist2strlist(sub_results)])
		result_table = insert_elements(result_table,sub_results,cont,1)
		cont = cont+1
		#except:
		#	print '\n\nERRORE: subject %s not correctly analyzed' %(str(sub_ID))
		
	# Reducing the result matix
	result_table[0][:]=range(result_table.shape[1])
	results = prepare_results_4_text_file(result_table,sub_list)
	print 'Saving the results...'
	make_text_file(results,output_file)		