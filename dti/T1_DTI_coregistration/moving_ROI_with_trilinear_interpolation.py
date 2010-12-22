#                                                                      Import section
# -----------------------------------------------------------------------------------
from nifti import NiftiImage # To manage niftii files
import numpy as np
import os,sys, platform


#                                                                 Loading setup paths
# -----------------------------------------------------------------------------------
sw_base = "/data/BI/Toolbox/software/"

if platform.system() == "Darwin":
	# Mac settings
	fsl_path  = sw_base + "fsl_maci_4.1.5/bin/"
	ANTS_path = sw_base + "ants_maci_svn460/bin/"
	c3d_path  = sw_base + "c3d_mac_maci_0.8.0/bin/"
	
	if platform.architecture()[0] == "64bit":
		afni_path = sw_base + "afni_maci_os10.5_64bit_v1710/"
	elif platform.architecture()[0] == "32bit":
		afni_path = sw_base + "afni_maci_os10.5_32bit_v1710/"
	else:
		print("Error: running on unsupported architecture " + platform.system())
		quit()
elif platform.system() == "Linux":
	# Linux settings
	fsl_path  = sw_base + "fsl_lnx_4.1.4/bin/"
	ANTS_path = sw_base + "ants_lnx64_svn460/bin/"
	c3d_path  = "No c3d folder for linux"
	afni_path = sw_base + "afni_lnx_2009_12_31_1431/"
else:
	print("Error: running on unsupported architecture " + platform.system())
	quit()


#                                                            Reading Input Parameters
# -----------------------------------------------------------------------------------
main_folder = sys.argv[1]
if (main_folder[-1] != os.sep):
	main_folder = main_folder + os.sep



#                                                             Defining local constant
# -----------------------------------------------------------------------------------
# Input files
DTI_ref_volume    = main_folder + "DTI_firstVolume.nii.gz"
ants_warp_FA_2_T1 = main_folder + "FA_2_T1_Warp.nii.gz"
ants_warp_T1_2_FA = main_folder + "FA_2_T1_InverseWarp.nii.gz"
ants_affine       = main_folder + "FA_2_T1_Affine.txt"
T1_data           = main_folder + "T1_cropped.nii.gz"
ROI_on_T1         = main_folder + "roi_T1_cropped.nii.gz"

# Output files
temporary_folder   = main_folder + "trilinear_coregistration_ants" + os.sep
ROI_tmp            = temporary_folder + 'ROI_original_%s.nii.gz'
ROI_warped_tmp     = temporary_folder + 'ROI_warped_%s.nii.gz'
final_map_filename = temporary_folder + 'ROI_warped_2_FA_ants_TRI.nii.gz'

mask_value        = 100
minimum_threshold = 50 #Minimum value in an interpolated voxel to be considered (maximum is 100)


#                                                                  Commands template
# -----------------------------------------------------------------------------------
ants_apply_warp = ANTS_path + 'WarpImageMultiTransform 3 %s %s -R %s -i %s %s' #(moving image, output image, reference image, affine file, warp file )


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
		mov_cmd = ants_apply_warp %(single_ROI_filename,moved_ROI_filename,DTI_ref_volume,ants_affine,ants_warp_T1_2_FA) #(moving image, output image, reference image, affine file, warp file )
		#mov_cmd = fsl_apply_warp %(DTI_ref_volume,single_ROI_filename,fsl_warp_T1_2_FA,moved_ROI_filename,'trilinear')
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
			if (map_value_vector.max()>minimum_threshold):
				indices = np.nonzero(map_value_vector==(map_value_vector.max()))
				# Writing the final_map
				if (len(indices)==1):
					final_map[c1][c2][c3]=ROI_index_list[indices[0][0]]
				else:
					print "WARNING: multiple candidate in a voxel"
					final_map[c1][c2][c3]=ROI_index_list[indices[0][0]]
final_map_obj = NiftiImage(final_map,header=template_obj.header)
final_map_obj.save(final_map_filename)
