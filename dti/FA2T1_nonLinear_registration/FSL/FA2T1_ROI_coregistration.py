""" This is a method to move a set of ROIs from T1 to DTI space using the trilinear
interpolation method in place of the nearest-neighbor.

- this method assumes that the coregistration files have already been computed using 
FSL.
- the method assumes a file name configuration and requires as input only the folder 
containing the files. See Defining local constant section.

NOTE: the Loading setup paths section call a file to load the configuration for the 
software paths.
"""

#                                                                 Loading setup paths
# -----------------------------------------------------------------------------------
run /mind_cart/DTI/DTI_script/settings.py


#                                                                      Import section
# -----------------------------------------------------------------------------------
from nifti import NiftiImage # To manage niftii files
import numpy as np
import os


#                                                            Reading Input Parameters
# -----------------------------------------------------------------------------------
main_folder = sys.argv[1]
if (main_folder[-1] != os.sep):
	main_folder = main_folder + os.sep


#                                                             Defining local constant
# -----------------------------------------------------------------------------------
# Input files
DTI_ref_volume   = main_folder + "DTI_firstVolume.nii.gz"
fsl_warp_FA_2_T1 = main_folder + "FA_2_T1_warpcoef.nii.gz"
fsl_warp_T1_2_FA = main_folder + "T1_2_FA_warpcoef.nii.gz"
T1_data          = main_folder + "T1_cropped.nii.gz"
ROI_on_T1        = main_folder + "roi_T1_cropped.nii.gz"

# Output files
temporary_folder   = main_folder + "trilinear_coregistration" + os.sep
ROI_tmp            = temporary_folder + 'ROI_original_%s.nii.gz'
ROI_warped_tmp     = temporary_folder + 'ROI_warped_%s.nii.gz'
final_map_filename = temporary_folder + 'ROI_warped_2_FA.nii.gz'

mask_value        = 100
minimum_threshold = 50 #Minimum value in an interpolated voxel to be considered (maximum is 100)

#                                                             Commands template
# -----------------------------------------------------------------------------------
fsl_apply_warp = fsl_path + 'applywarp --ref=%s --in=%s --warp=%s --out=%s --interp=%s' #(reference volume, input volume, warp file, output_file, interpolation method {nn,trilinear,sinc})


#                                                                        Main section
# -----------------------------------------------------------------------------------
if not(os.path.exists(temporary_folder)):
	os.mkdir(temporary_folder)
	os.chmod(temporary_folder,0777)

# Separating and moving each ROI
# ------------------------------
ROI_map_obj = NiftiImage(ROI_on_T1)
ROI_map     = ROI_map_obj.getDataArray()
file_list    = []
for cont in range(ROI_map.max()):
	# Single ROI estraction
	ROI_index = cont +1
	indices = np.nonzero(ROI_map == ROI_index)
	if (len(indices[0])>0):
		single_ROI_map = np.zeros(ROI_map.shape)
		single_ROI_map [indices] = mask_value
		# Saving the file
		single_ROI_obj = NiftiImage(single_ROI_map,header=ROI_map_obj.header)
		single_ROI_filename = ROI_tmp %(str(ROI_index))
		single_ROI_obj.save(single_ROI_filename)
		# Moving the ROI
		moved_ROI_filename = ROI_warped_tmp %(str(ROI_index))
		mov_cmd = fsl_apply_warp %(DTI_ref_volume,single_ROI_filename,fsl_warp_T1_2_FA,moved_ROI_filename,'trilinear')
		result = os.system(mov_cmd)
		if (result == 0):
			file_list.append([moved_ROI_filename,ROI_index])
		else:
			print "WARNING: can't coregister file: "+single_ROI_filename
del(ROI_map_obj)
del(ROI_map)
del(single_ROI_obj)
del(single_ROI_map)
	
# Reading all coregistered ROIs
# -----------------------------
map_list = []
ROI_index_list = []
for cont in range(len(file_list)):
	ROI_map_obj = NiftiImage(file_list[cont][0])
	map_list.append(ROI_map_obj.getDataArray())
	ROI_index_list.append(file_list[cont][1])
del(ROI_map_obj)

# Computing a consensus ROI map
# -----------------------------
template_obj = NiftiImage(file_list[0][0])
template_map = template_obj.getDataArray()
dimension = template_map.shape
nMap = len(map_list)
final_map = np.zeros(dimension)
for c1 in range(dimension[0]):
	print(c1)
	for c2 in range(dimension[1]):
		for c3 in range(dimension[2]):
			# Generating the vector with the map values
			map_value_vector = np.zeros(nMap)
			for cm in range(nMap):
				map_value_vector[cm]=map_list[cm][c1][c2][c3]
			# Chosing the ROI
			if not(map_value_vector.max()>minimum_threshold):
				indices = np.nonzero(map_value_vector==(map_value_vector.max()))
				# Writing the final_map
				if (len(indices)==1):
					final_map[c1][c2][c3]=ROI_index_list[indices[0]]
				else:
					print "WARNING: multiple candidate in a voxel"
					final_map[c1][c2][c3]=ROI_index_list[indices[0]]
					
# Saving the consensus ROI map
# ----------------------------
final_map_obj = NiftiImage(final_map,header=template_obj.header)
final_map_obj.save(final_map_filename)
