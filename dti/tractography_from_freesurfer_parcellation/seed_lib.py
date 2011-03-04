"""
SEED LIBRARY

Contains some tools useful to perform tractography.

 --------------------------------------------------------------------------------------------
 				GET NEIGHBORS
 Given a 3D volume and the coordinates (r,c,s) of one voxel it returns the values of the 
 neighbors.
		
 --------------------------------------------------------------------------------------------
 				MOVE LABELS 2 WM
 Given a freesurfer parcellated volume this method grows each ROI in the withe matter by one
 voxel.

 --------------------------------------------------------------------------------------------
 				MOVE LABELS 2 WM 2
 Given a freesurfer parcellated volume, a list of white matter labels and a list of growin 
 labels this method grow the growin labels by one voxel in the white matter ROIs

 --------------------------------------------------------------------------------------------
 				MOVE LABEL FILE 2 WM
 It's like move labels 2 wm 2, but starting from a file in place of a volume. Results are 
 saved in another file.

 --------------------------------------------------------------------------------------------
 				MASKING MAP
 Remove some ROIs from a ROI map. fileIN conaints le name of the file with the ROI map. 
 labels2mask is the list of the label to mask.

 --------------------------------------------------------------------------------------------
 				EXTRACT SEED FROM A FILE
 From a ROI map this method extracts the labels in a given list and save each ROI in a file.

 --------------------------------------------------------------------------------------------
 				EXTRACT ALL ROIs
 From a ROI map this method extracts all non zero ROIs and save each one in a file.

 --------------------------------------------------------------------------------------------
 				GENERATE SEED MASKS
 OLD METHOD. From a ROI map this method extracts all non zero ROIs using C3D and save each 
 one in a file.

 --------------------------------------------------------------------------------------------
 				GENERATE SEED MASKS (uses fsl)
 OLD METHOD. From a ROI map this method extracts all non zero ROIs using FSL and save each 
 one in a file.
"""	



# IMPORT
import os, sys
import numpy as np
#import nibabel as nb
from nifti import NiftiImage
#from pylab import *

# Default Parameters
all_seeds_mask = 'all_seeds_mask.nii.gz'
comp_mask      = 'components_mask.nii.gz'
seed_masks     = 'seed_'
file_sep       = os.sep




# --------------------------------------------------------------------------------------------
# 				GET NEIGHBORS
#
# Given a 3D volume and the coordinates (r,c,s) of one voxel it returns the values of the 
# neighbors.
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
# 				MOVE LABELS 2 WM
# 
# Given a freesurfer parcellated volume this method grows each ROI in the withe matter by one
# voxel.
def move_labels2wm(volume):
	freesurfer_lh_wm = 2  # Freesurfer label for the left hemisphere white matter.
	freesurfer_rh_wm = 41 # Freesurfer label for the right hemisphere white matter.
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
						# All neighbours are of the same type of the [r,c,s] voxel
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
# 				MOVE LABELS 2 WM 2
# 
# Given a freesurfer parcellated volume, a list of white matter labels and a list of growin 
# labels this method grow the growin labels by one voxel in the white matter ROIs

def move_labels2wm_2(volumeIN,whiteMatterLabels,nonGrowingLabels):
# Move all labels in a volume to white matter.
# Note: it supposes to have a Desikan-Kiliany parcellation 
# Note: whiteMatterLabels contains all labels which have to be considered white matter
	# Mask all white matter labels (i.e. all white matter voxel label is set to -1)
	volume = volumeIN.copy()
	for cont in range(len(whiteMatterLabels)):
		ind = (volumeIN==whiteMatterLabels[cont]).nonzero()
		volume[ind]=-1;
	# Mask non growing labels (i.e. non growing voxel label is set to -1)
	for cont in range(len(nonGrowingLabels)):
		ind = (volumeIN==nonGrowingLabels[cont]).nonzero()
		volume[ind]=-2;
	# Prepare the variables for the results
	volume_shape = volume.shape
	results = np.zeros(volume_shape)
	for s in range(volume_shape[2]):
		print 'slice %s' %(str(s))
		for c in range(volume_shape[1]):
			for r in range(volume_shape[0]):
				if (volume[r][c][s] in [-1]):
					# The voxel is a voxel which has to be considered white matter					
					# Get the neighbors
					neighbors = get_neighbors(volume,r,c,s)
					if len((neighbors==volume[r,c,s]).nonzero()[0]) == len(neighbors.ravel()):
						# The neighbors are all of the same type of the main voxels.
						results[r][c][s] = volume[r][c][s]
					else:
						# Look at the largest non-wm population
						elem_set = set(neighbors.ravel().tolist())
						elem_list = np.zeros(len(elem_set))
						pop_list = np.zeros(len(elem_set))
						cont = 0
						for value in elem_set:
							elem_list[cont] = value
							if (value in [-1]):
								pop_list[cont] = 0
							else:
								pop_list[cont] =len((neighbors == value).nonzero()[0])
							cont = cont +1
						if (pop_list.max()==0):
							results[r][c][s] = volume[r][c][s]
						else:
							results[r][c][s]=elem_list[find(pop_list==pop_list.max())][0]
				else:
					results[r][c][s] = volume[r][c][s]
	# Restoring white matter labels
	ind = (results==-1).nonzero()
	results[ind]=volumeIN[ind]
	# Restoring non growing region labels
	ind = (results==-2).nonzero()
	results[ind]=volumeIN[ind]
	return results



# --------------------------------------------------------------------------------------------
# 				MOVE LABEL FILE 2 WM
#
# It's like move labels 2 wm 2, but starting from a file in place of a volume. Results are 
# saved in another file.

def move_labelFile2wm(fileIn,fileOut,whiteMatterLabels,nonGrowingLabels):
# Read a file, move all labels provided by freesurfer to white matter and save the result in 
# an other file.
# - whiteMatterLabels contains the labels which have to be considered white matter, i.e. where the 

# Note: it supposes to have a Desikan-Kiliany parcellation
# Note: whiteMatterLabels
	#dataIN = nb.load(fileIn)
	#volume = dataIN.get_data()
	dataIN = NiftiImage(fileIn)
	volume = dataIN.getDataArray()
	new_volume = move_labels2wm(volume,whiteMatterLabels,nonGrowingLabels)
	#new_data = nb.nifti1.Nifti1Image(new_volume,dataIN.get_affine(),dataIN.get_header(),extra = None,file_map = dataIN.file_map)
	#nb.save(new_data,fileOut)
	new_data = NiftiImage(new_volume,dataIN.header)
	new_data.setFilename(fileOut)
	new_data.save()



# --------------------------------------------------------------------------------------------
# 				MASKING MAP
# 
# Remove some ROIs from a ROI map. fileIN conaints le name of the file with the ROI map. 
# labels2mask is the list of the label to mask.
def masking_map(fileIN,fileOUT,labels2mask):
# Read a file (supposed to be a parcellation file) and eliminate all ROI in labels2mask
	#dataIN = nb.load(fileIN)
	#volumeIN = dataIN.get_data() # With nibabel
	dataIN = NiftiImage(fileIN)
	volumeIN = dataIN.getDataArray()
	volumeOUT = volumeIN.copy()
	for cont in range(len(labels2mask)):
		ind = (volumeIN == labels2mask[cont]).nonzero()
		volumeOUT[ind]=0
	#dataOUT = nb.nifti1.Nifti1Image(volumeOUT,dataIN.get_affine(),dataIN.get_header(),extra = None,file_map = dataIN.file_map)
	#nb.save(dataOUT,fileOUT) # With nibabel
	dataOUT = NiftiImage(volumeOUT,dataIN.header)
	dataOUT.setFilename(fileOUT)
	dataOUT.save()


	
# --------------------------------------------------------------------------------------------
# 				EXTRACT SEED FROM A FILE
#
# From a ROI map this method extracts the labels in a given list and save each ROI in a file.
def extract_seed_from_file(fileIn,fileOut_tmp,seed_list):
	# NOTE: seed_list can be a list
	cmd_tmp = '/data/BI/Toolbox/software/c3d_mac_maci_0.8.0/bin/c3d %s -threshold %s %s 1 0 -o %s' #order: fileIn, threshold threshold fileOut
	for cont in range(len(seed_list)):
		seedNo = seed_list[cont]
		fileOut = fileOut_tmp + '_' + str(seedNo)+'.nii.gz'
		cmd = cmd_tmp %(fileIn,str(seedNo),str(seedNo),fileOut)
		os.system(cmd)



# --------------------------------------------------------------------------------------------
# 				EXTRACT ALL ROIs
#
# From a ROI map this method extracts all non zero ROIs and save each one in a file.
def extract_all_ROIs(fileIN,fileOUT_tmp):
# From a given file (supposed to be a parcellation) extract all ROIs
	#dataIN = nb.load(fileIN) # With nibabel
	#volumeIN = dataIN.get_data() # with nibabel
	dataIN = NiftiImage(fileIN)
	volumeIN = dataIN.getDataArray()
	element_list = np.unique(volumeIN)
	for cont in range(len(element_list)):
		if element_list[cont]>0:
			ind = (volumeIN == element_list[cont]).nonzero()
			volumeOUT = np.zeros(volumeIN.shape)
			volumeOUT[ind]=1
			fileOUT = fileOUT_tmp %(str(np.int(element_list[cont])))
			#dataOUT = nb.nifti1.Nifti1Image(volumeOUT,dataIN.get_affine(),dataIN.get_header(),extra = None,file_map = dataIN.file_map)
			#nb.save(dataOUT,fileOUT) # With nibabel
			dataOUT = NiftiImage(volumeOUT,dataIN.header)
			dataOUT.setFilename(fileOUT)
			dataOUT.save()




# --------------------------------------------------------------------------------------------
# 				GENERATE SEED MASKS
# 
# OLD METHOD. From a ROI map this method extracts all non zero ROIs using C3D and save each 
# one in a file.
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
	#fileObj = nb.load(output_folder + comp_mask)
	#data = fileObj.get_data()
	fileObj = NiftiImage(output_folder + comp_mask)
	data = fileObj.getDataArray()
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
# OLD METHOD. From a ROI map this method extracts all non zero ROIs using FSL and save each 
# one in a file.
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