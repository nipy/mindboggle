"""
generate_seed_masks(binary_mask,output_folder):

Take the binary mask for seed ROIs and divide it in one mask for each seed.
Seeds are separate on the basis of connection by using convert3D

INPUT:
  * binary mask
  * output folder
  
---------------------------------------------------------------------------
move_seed2dti_space(seed_file,DTI_file,trans_file,output_folder):

Move the seed mask to the DTI space and save it in the specified folder
  
"""	



# IMPORT
import os, sys
import numpy as np
import nibabel as nb
from pylab import *

# Default Parameters
all_seeds_mask = 'all_seeds_mask.nii.gz'
comp_mask      = 'components_mask.nii.gz'
seed_masks     = 'seed_'
file_sep       = os.sep

freesurfer_lh_wm = 2  # Freesurfer label for the left hemisphere white matter.
freesurfer_rh_wm = 41 # Freesurfer label for the left hemisphere white matter.

# --------------------------------------------------------------------------------------------
# 				GET NEIGHBORS
# 
def get_neighbors(volume,r,c,s):
	volume_shape = volume.shape
	
	edges = 0;
	if ((r <= 0)or(r>=volume_shape[0])):
		edges = edges+1
	if ((c <= 0)or(c>=volume_shape[1])):
		edges = edges+1
	if ((s <= 0)or(s>=volume_shape[2])):
		edges = edges+1
		
	if (edges == 0):
		neighbors = np.zeros((27,1))
	elif (edges == 1):
		neighbors = np.zeros((18,1))
	elif (edges == 2):
		neighbors = np.zeros((12,1))
	elif (edges == 3):
		neighbors = np.zeros((8,1))
		
	cont = 0
	for off_r in [-1,0,1]:
		if ((r+off_r >= 0)and(r+off_r <= volume_shape[0])):
			for off_c in [-1,0,1]:
				if ((c+off_c >= 0)and(c+off_c <= volume_shape[1])):
					for off_s in [-1,0,1]:
						if ((s+off_s >= 0)and(s+off_s <= volume_shape[2])):
							neighbors[cont] = volume[r+off_r][c+off_c][s+off_s]
							cont = cont+1
				
	return neighbors


		
# --------------------------------------------------------------------------------------------
# 				MOVE CORTICAL REGION TO WHITE MATTER
# 

def move_labels2wm(volume):
	volume_shape = volume.shape
	results = np.zeros(volume_shape)
	for s in range(volume_shape[2]):
		print 'slice %s' %(str(s))
		for c in range(volume_shape[1]):
			for r in range(volume_shape[0]):
				if (volume[r][c][s] in [freesurfer_lh_wm,freesurfer_rh_wm]):
					# The voxel is a wm voxel in the left hemisphere
					
					# Get the neighbors
					neighbors = get_neighbors(volume,r,c,s)
					if find(neighbors==volume[r,c,s]).size == len(neighbors.ravel()):
						results[r][c][s] = volume[r][c][s]
					else:
						# Look at the largest non-wm population
						elem_set = set(neighbors.ravel().tolist())
						elem_list = np.zeros(len(elem_set))
						pop_list = np.zeros(len(elem_set))
						cont = 0
						for value in elem_set:
							elem_list[cont] = value
							if (value in [freesurfer_lh_wm,freesurfer_rh_wm,0]):
								pop_list[cont] = 0
							else:
								pop_list[cont] =find(neighbors == value).size
							cont = cont +1
						if (pop_list.max()==0):
							results[r][c][s] = volume[r][c][s]
						else:
							results[r][c][s]=elem_list[find(pop_list==pop_list.max())][0]
				else:
					results[r][c][s] = volume[r][c][s]
	
	return results



# --------------------------------------------------------------------------------------------
# 				MOVE CORTICAL REGION TO WHITE MATTER
# 

def move_labels2wm_2(volume):
	volume_shape = volume.shape
	results = np.zeros(volume_shape)
	for s in range(volume_shape[2]):
		print 'slice %s' %(str(s))
		for c in range(volume_shape[1]):
			for r in range(volume_shape[0]):
				if (volume[r][c][s] in [freesurfer_lh_wm,freesurfer_rh_wm]):
					# The voxel is a wm voxel in the left hemisphere
					
					# Get the neighbors
					neighbors = get_neighbors(volume,r,c,s)
					if find(neighbors==volume[r,c,s]).size == len(neighbors.ravel()):
						results[r][c][s] = volume[r][c][s]
					else:
						# Look at the largest non-wm population
						elem_set = set(neighbors.ravel().tolist())
						elem_list = np.zeros(len(elem_set))
						pop_list = np.zeros(len(elem_set))
						cont = 0
						for value in elem_set:
							elem_list[cont] = value
							if (value in [freesurfer_lh_wm,freesurfer_rh_wm,0]):
								pop_list[cont] = 0
							else:
								pop_list[cont] =find(neighbors == value).size
							cont = cont +1
						if (pop_list.max()==0):
							results[r][c][s] = volume[r][c][s]
						else:
							results[r][c][s]=elem_list[find(pop_list==pop_list.max())][0]
					
	return results



# --------------------------------------------------------------------------------------------
# 				MOVE CORTICAL REGION TO WHITE MATTER
# 

def move_labelFile2wm(fileIn,fileOut):
	dataIN = nb.load(fileIn)
	volume = dataIN.get_data()
	new_volume = move_labels2wm(volume)
	new_data = nb.nifti1.Nifti1Image(new_volume,dataIN.get_affine(),dataIN.get_header(),extra = None,file_map = dataIN.file_map)
	nb.save(new_data,fileOut)


	
# --------------------------------------------------------------------------------------------
# 				EXTRACT SEED FROM A FILE
# 
def extract_seed_from_file(fileIn,fileOut_tmp,seed_list):
	# NOTE: seed_list can be a list
	cmd_tmp = '/data/BI/Toolbox/software/c3d_mac_maci_0.8.0/bin/c3d %s -threshold %s %s 1 0 -o %s' #order: fileIn, threshold threshold fileOut
	for cont in range(len(seed_list)):
		seedNo = seed_list[cont]
		fileOut = fileOut_tmp + '_' + str(seedNo)+'.nii.gz'
		cmd = cmd_tmp %(fileIn,str(seedNo),str(seedNo),fileOut)
		os.system(cmd)
		


# --------------------------------------------------------------------------------------------
# 				GENERATE SEED MASKS
# 

def generate_seed_masks(file_name,output_folder):
	# checking input parameters
	print 'Checking input data ...'
	if not(os.path.isfile(file_name)):
		print 'ERROR: file %s not found' %(file_name)
		raise
	
	if output_folder[len(output_folder)-1] != file_sep:
		output_folder = output_folder + file_sep

	# Create the destination folder, if necessary
	if not(os.path.exists(output_folder)):
		print 'Generating the destination folder...'
		cmd = 'mkdir -p -m 777 ' + output_folder
		os.system(cmd)
	else:
		print 'Destination folder already found'

	# Elaborate the binary image
	# Look for connected components. Each connected component is labelled with an integer
	cmd = '/data/BI/Toolbox/software/c3d_mac_maci_0.8.0/bin/c3d ' + file_name + ' -connected-components -o ' + output_folder + comp_mask
	os.system(cmd)

	# Reading the component mask to compute the number of components
	fileObj = nb.load(output_folder + comp_mask)
	data = fileObj.get_data()
	nComp = np.max(data)
	del data
	del fileObj

	# Generating the mask for each seed
	for comp in range(nComp):
		compNumber = comp +1 # comp start from 0, components start from 1
		cmd = '/data/BI/Toolbox/software/c3d_mac_maci_0.8.0/bin/c3d ' + output_folder + comp_mask + ' -threshold ' + str(compNumber) + ' ' + str(compNumber) + ' 1 0 -o ' + output_folder + seed_masks + str(compNumber) + '_mask.nii.gz'
		os.system(cmd)
	
	cmd = 'chmod 666 ' + output_folder + '*.nii.gz'
	os.system(cmd)
	


# --------------------------------------------------------------------------------------------
# 				GENERATE SEED MASKS (uses fsl)
# 

def move_seed2dti_space(seed_file,DTI_file,trans_file,output_folder):
	# checking input parameters
	print 'Checking input data ...'
	if not(os.path.isfile(seed_file)):
		print 'ERROR: file %s not found' %(seed_file)
		raise
	
	if not(os.path.isfile(DTI_file)):
		print 'ERROR: file %s not found' %(seed_file)
		raise
			
	if not(os.path.isfile(trans_file)):
		print 'ERROR: file %s not found' %(trans_file)
		raise
		
	if output_folder[len(output_folder)-1] != file_sep:
		output_folder = output_folder + file_sep

	# Create the destination folder, if necessary
	if not(os.path.exists(output_folder)):
		print 'Generating the destination folder...'
		cmd = 'mkdir -p -m 777 ' + output_folder
		os.system(cmd)
	else:
		print 'Destination folder already found'
		
	# fsl:FLIRT call to move the mask to DTI space.
	print 'Moving the seed mask to DTI space...'
	cmd = 'flirt -in ' + seed_file + ' -ref ' + DTI_file + ' -out ' + output_folder + all_seeds_mask + ' -init ' + trans_file + ' -applyxfm'
	os.system(cmd)