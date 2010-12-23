"""
SULCAL PITS LIBRARY
Include instruments to read e move the sulcal pits between freesurfer, T1 and DTI spaces.

-----------------------------------------------------------------------------------------
IMPORT_SULCAL_PITS(input_folder,output_folder):
Import the sulcal pits from the freesurfer space to T1 space. 
The following operation are performed:
- sulcal pits are highlighted to 250 (freesurfer space)
- highlighted sulcal pits are superimposed to the brain image (freesurfer space)
- the new image (brain+sulcalPits) is resampled in the T1 space (Nerarest Neighbor)
- threshold is applied to extract the sulcal pits in the T1 space (T1 space)
- sulcal pits are dilated by 3 voxels (T1 space)

-----------------------------------------------------------------------------------------
IMPORT_FROM_TEXT_FILE(fileName):
Read the folder list from a text file e for each folder call import_sulcal_pits

-----------------------------------------------------------------------------------------
MOVE_SULCALPITS_FROM_T1_TO_DTI(T1_sulcalPits,output_folder,ants_inverseWarpFile,ants_affine_file)
Move the sulcal pits from T1 space to DTI one. 
1) Dilate the sulcalPits by 3. (T1 space)
2) Moving the sulcalPits to DTI space (T1 -> DTI)
3) Connected components analysis:
   3.1) Find the center for the connected component
   3.2) Dilate by 3 the connected component center
   3.3) Save the seed file
"""

# ------------
# IMPORT PHASE
# ------------
import os, sys, string


# -------
# METHODS
# -------

# ---------------------------------------------------------------------------------------

def move_sulcalPits_from_T1_to_DTI(T1_sulcalPits,output_folder,DTI_ref_volume,ants_inverseWarpFile,ants_affine_file,seed_dilatation_radius):
	# IMPORT
	import nibabel as nb
	import numpy as np
	import settings
	
	# SETTING PATHS
	paths = settings.settings()
	fsl_path  = paths[0]
	ants_path = paths[1]
	afni_path = paths[2]
	c3d_path  = paths[3]
	if c3d_path=='':
		print "WARNING: convert3D (c3d) not installed for this system"
		quit()
	
	#COMMANDS
	cmd_comp_tmp = c3d_path + 'c3d %s -connected-components -o %s' # Order (input_file,output_file)
	cmd_seed_mask_tmp = c3d_path + 'c3d %s -threshold %s %s 1 0 -o %s' # Order (input_file, min_threshold, max_threshold,output_file)

	cmd_dilate_tmp = ants_path + 'ImageMath 3 %s MD %s %s' #(output_file, input_file, dilatation_radius)
	cmd_moving_tmp = ants_path + 'WarpImageMultiTransform 3 %s %s -R %s -i  %s %s --use-NN'  # Order (in_file, out_file, DTI_first_volume, 
		
	if (output_folder[-1] != os.sep):
		output_folder=output_folder + os.sep
	
	seed_file_tmp = output_folder + 'seed_%s_mask.nii.gz' 
	T1dilated = output_folder + 'T1_sulcalPits_dilated.nii.gz'
	T1dilated2DTI = output_folder + 'T12DTI_sulcalPits_dilated.nii.gz'
	T1dilated2DTIcomponents = output_folder + 'T12DTI_sulcalPits_dilated_components.nii.gz'
	DTIsulcal = output_folder + 'DTI_sulcalPits.nii.gz'
	DTIdilated = output_folder + 'DTI_sulcalPits_dilated.nii.gz'
	dilatation_radius = '2' # Radius for the dilatation steps
	
	if not(os.path.exists(output_folder)):
		os.mkdir(output_folder,0777)
	
	# 1) Dilate the sulcalPits by 3. (T1 space)
 	dilate_cmd = cmd_dilate_tmp %(T1dilated,T1_sulcalPits,dilatation_radius) #(output_file, input_file, dilatation_radius)
	os.system(dilate_cmd)

	# 2) Moving the sulcalPits to DTI space (T1 -> DTI)
	moving_cmd = cmd_moving_tmp %(T1dilated,T1dilated2DTI,DTI_ref_volume,ants_affine_file,ants_inverseWarpFile) # Order (in_file, out_file, DTI_first_volume, affineFile, inverseWarpFile)
	os.system(moving_cmd)
		
	# 3) Generating connected component map
	comp_cmd = cmd_comp_tmp %(T1dilated2DTI,T1dilated2DTIcomponents) # Order (input_file,output_file)
	os.system(comp_cmd)
	
	# 4) Connected components analysis:	
	# Reading the component mask to compute the number of components
	fileObj = nb.load(T1dilated2DTIcomponents)
	components_map = fileObj.get_data()
	nComp = np.max(components_map)
	DTI_sulcalPits_map = np.zeros(components_map.shape)
	
	# Generating the mask for each seed
	for comp in range(nComp):
		# single component analysis
		compNumber = comp +1 # comp start from 0, components start from 1
		# 4.1) Find the center for the connected component
		seed_map = np.zeros(components_map.shape)
		seed_mask = components_map == compNumber
		seed_indices = seed_mask.nonzero()
		mass_center = []
		for cont in range(len(seed_indices)):
			mass_center = mass_center + [np.round(np.mean(seed_indices[cont]))]
		seed_mass_center = np.round(mass_center)
		seed_map[mass_center[0]][mass_center[1]][mass_center[2]]=1
		DTI_sulcalPits_map[mass_center[0]][mass_center[1]][mass_center[2]]=1	
		seed_file = seed_file_tmp %(str(compNumber))
		seed_obj = nb.nifti1.Nifti1Image(seed_map,fileObj.get_affine(),fileObj.get_header(),extra = None,file_map = fileObj.file_map)
		# 4.2) Save the seed file
		nb.save(seed_obj,seed_file)
		# 4.3) Dilate by 3 the connected component center
		dilate_cmd = cmd_dilate_tmp %(seed_file,seed_file,seed_dilatation_radius) #(output_file, input_file, dilatation_radius)
		os.system(dilate_cmd)
		# 4.4) clear temporary variables	
		del(seed_obj)
		del(seed_file)
		del(seed_indices)
		del(seed_mask)
		del(seed_map)
		del(mass_center)
		del(seed_mass_center)
	
	# Saving a summary map	
	map_obj = nb.nifti1.Nifti1Image(DTI_sulcalPits_map,fileObj.get_affine(),fileObj.get_header(),extra = None,file_map = fileObj.file_map)
	nb.save(map_obj,DTIsulcal)




# ---------------------------------------------------------------------------------------
def import_sulcal_pits(input_folder,output_folder):
	# SETTING PATHS
	paths = settings.settings()
	fsl_path  = paths[0]
	ants_path = paths[1]
	afni_path = paths[2]
	c3d_path  = paths[3]
	
	# Summary
	print 'Subject: %s' %(sub)
	print 'Subject: %s' %(input_folder)
	print 'Subject: %s' %(output_folder)
	
	# Making the output folder
	os.system('mkdir -p -m 777 ' +output_folder)

	# Commands execution
	cmd=ants_path+'ImageMath 3 ' + output_folder + 'pits_times_250.nii.gz m ' + input_folder + 'brain_SulcalPits.nii.gz 250'
	os.system(cmd) # Highlight the sulcal pits
	
	cmd=ants_path+'ImageMath 3 ' + output_folder + 'pits_plus_brains.nii.gz + ' + input_folder + 'brain.nii.gz ' + output_folder + 'pits_times_250.nii.gz'
	os.system(cmd) # sum the sulcal pits to the anatomical image
	
	cmd='mri_convert ' + output_folder + 'pits_plus_brains.nii.gz ' + output_folder + 'temp.nii.gz -rl ' + input_folder + 'brain_n3.nii.gz -rt nearest -ns 1 -nc' 
	os.system(cmd) # Reslice the image from freeSurfer space to T1 space
	
	cmd=c3d_path+'c3d ' + output_folder + 'temp.nii.gz -threshold -inf 249 0 1 -o ' + output_folder + 'brain_n3_sulcalPits.nii.gz'
	os.system(cmd) # Select the sulcal pits using a threshold
	
	cmd=ants_path+'ImageMath 3 ' + output_folder + 'pits_dilate3.nii.gz MD ' + output_folder + 'brain_n3_sulcalPits.nii.gz 3'
	os.system(cmd) # Dilate the sulcal pits for tractography
	
	# Permissions
	cmd='chmod 666 '+output_folder+'*.*'
	os.system(cmd)




# ---------------------------------------------------------------------------------------
def import_from_text_file(fileName):

	fileObj = open(fileName,'r')
	for line in fileObj:
		line = line.strip()
		arguments = line.split(' ')
		sub = arguments[0]
		input_folder = arguments[1]
		#folder=line[:len(line)-1] # Remove the newline character
		#sub=line[0:5] # Remove the newline character
		#input_folder=line[5:len(line)-1] # Remove the newline character
		output_folder = '/mind_cart/DTI/' + sub + '/SulcalPits/'
		import_sulcal_pits(input_folder,output_folder)
	fileObj.close()
