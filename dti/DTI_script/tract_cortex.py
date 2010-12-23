import nibabel as nb
import numpy as np
import os,sys
import seed_lib as sl

# -----------------------------------------------------------------------------------------
# IMPORT PHASE
# Subject folder, containing all necessary folders and files
main_folder = sys.argv[1]
if (main_folder[-1] != os.sep):
	main_folder=main_folder + os.sep
	
index = sys.argv[2]

# -----------------------------------------------------------------------------------------
# CONSTANT DEFINITION
#lh_cortex = [1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035]
lh_cortex = [1008, 1022]

#rh_cortex = [2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029, 2030, 2031, 2032, 2033, 2034, 2035]
rh_cortex = [2008, 2022]

#other_structures = [10, 11, 12, 13, 16, 17, 18, 26, 28, 49, 50, 51, 52, 53, 54, 58, 60]
other_structures = [16, 28]

seed_list= lh_cortex + rh_cortex + other_structures

target_list = seed_list 
# Seeds and targets are the same

cortex_threshold = 999 
# All ROIs below this threshold are considered in the tractography domain (i.e. streamlines allowded voxels

csf_threshold = 0.8
# threshold to binarize the csf probability map. csf_mask is used as exlusion mask in tractography



# -----------------------------------------------------------------------------------------
# COMMANDS TEMPLATE

# Making streamlines domain
cmd_tract_domain_tmp = 'c3d %s -threshold 1 %s 1 0 -o %s' #order: fileIn, threshold threshold fileOut

# Performing tractography
cmd_tract_tmp = 'probtrackx --mode=seedmask -x %s  -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --avoid=%s --forcedir --opd -s %s/merged -m %s  --dir=%s --targetmasks=%s --os2t' #Order (seed_file, exclusion_mask, bedpostX, domain, output_folder, target_file)

# Moving image from T1 to DTI space
cmd_move_csf1_tmp = 'WarpImageMultiTransform 3 %s %s -R %s -i  %s %s --use-NN'  # Order (T1_csf_file, out_file, DTI_first_volume, affineFile, inverseWarpFile)

# Extracting ROI from parcellation file
cmd_move_csf2_tmp = 'fslmaths %s -thr %s -bin %s' # Order (DTI_csf_prob, csf_threshold, DTI_csf_mask)



# -----------------------------------------------------------------------------------------
# FOLDERS AND FILES ORGANIZATION
# INPUT
# Tractography
bedpostX_folder = '%sfsl_data.bedpostX' 

# Coregistration
Ants_affine = '%soriginal_data/DTI2T1_ants_T_Affine.txt' %(main_folder)
Ants_invWarp = '%soriginal_data/DTI2T1_ants_T_InverseWarp.nii.gz' %(main_folder)
Ants_Warp = '%soriginal_data/DTI2T1_ants_T_Warp.nii.gz' %(main_folder) 

# Seeds/Targets
labels_DTI_file = '%sROImaps/mri___aparc+aseg_DTI.nii.gz' %(main_folder) 

# Other
T1_csf_map_file = '%sfsl_segmentation/fsl_seg_prob_0.nii.gz' %(main_folder)
DTI_first_volume = '%soriginal_data/DTI_firstVolume.nii.gz' %(main_folder)

# OUTPUT
# Tractography
tract_domain = 'tract_domain.nii.gz' # Where streamlines can propagate
target_file = 'targets.txt' # Text files containing the target file list for tractography
tract_out_folder_tmp = '%stract_cortex/tract_from_%s'  #Template for tract results

# Seeds/Targets
whole_ROI_file = 'cortical_parcellation.nii.gz'
moved_ROI_file = 'dilated_parcellation.nii.gz'
ROI_name_tmp = 'ROI_%s.nii.gz' # Template for seed ROI
ROI_folder = '%stract_cortex/ROIs/' # Seed/target folder

# Coregistration
DTI_csf_mask = '%sDTI_csf_mask80.nii.gz' %(ROI_folder)
DTI_csf_prob = '%sDTI_csf_prob.nii.gz' %(ROI_folder)

# Other
temp_folder_IN = '/tmp/tract%s/IN/' %(index)# Temporary folder to use while performing tractography
temp_folder_OUT = '/tmp/tract%s/OUT/' %(index) # Temporary folder to use while performing tractography
temp_result_folder = '/tmp/tract%s/OUT/tract_cortex' %(index)
final_result_folder = '%stract_cortex/' %(main_folder)

label_WM_file = 'corticalROI2WM.nii.gz'








# -----------------------------------------------------------------------------------------
# PERFORMING TRACTOGRAPHY

# 1) Copying data in the temporary folder
# ---------------------------------------
print 'Copying data in the temporary folder...'
# DebpostX data
cmd = 'mkdir -m 777 -p ' + bedpostX_folder %(temp_folder_IN)
os.system(cmd)
cmd = 'cp ' + bedpostX_folder %(main_folder) + '/*.* ' + bedpostX_folder %(temp_folder_IN)
os.system(cmd)

cmd = 'mkdir -m 777 -p ' + ROI_folder %(temp_folder_OUT)
os.system(cmd)
cmd = 'cp ' +  labels_DTI_file + ' ' + ROI_folder %(temp_folder_OUT) + whole_ROI_file
os.system(cmd)


# 2) Generating the cortex ROI
print 'Preparing seed and target ROIs...'
sl.move_labelFile2wm(ROI_folder %(temp_folder_OUT) + whole_ROI_file,ROI_folder %(temp_folder_OUT) + moved_ROI_file) # Compute the projection of the cortical ROI to the white matter

# Extracting the single ROI mask to use as seed and target
sl.extract_seed_from_file(ROI_folder %(temp_folder_OUT) + moved_ROI_file,ROI_folder %(temp_folder_OUT)+'ROI',lh_cortex) # Dilated ROI
sl.extract_seed_from_file(ROI_folder %(temp_folder_OUT) + moved_ROI_file,ROI_folder %(temp_folder_OUT)+'ROI',rh_cortex) # Dilated ROI
sl.extract_seed_from_file(ROI_folder %(temp_folder_OUT) + whole_ROI_file,ROI_folder %(temp_folder_OUT)+'ROI',other_structures) # NOTE: internal structures aren't dilated


# 3) Making the tractography domain
print 'Generating the tractography domain'
# NOTE: seed ROI must be added to tract domain in order to run tractography.
cmd = cmd_tract_domain_tmp %(ROI_folder %(temp_folder_OUT) + whole_ROI_file,cortex_threshold,ROI_folder %(temp_folder_OUT) + tract_domain)
os.system(cmd)


# 4) Generating CSF mask
print 'Generating the CSF mask'
# Move the csf probability map to DTI space and binarize it.
cmd_move_csf1 = cmd_move_csf1_tmp %(T1_csf_map_file, DTI_csf_prob %(temp_folder_OUT), DTI_first_volume, Ants_affine, Ants_invWarp)
os.system(cmd_move_csf1)
cmd_move_csf2 = cmd_move_csf2_tmp %(DTI_csf_prob %(temp_folder_OUT), str(csf_threshold), DTI_csf_mask %(temp_folder_OUT))
os.system(cmd_move_csf2)


# 5) Generating target text file
print 'Generating target file list'
target_file_lines = []
for target_ind in range(len(target_list)):
	seedNo = target_list[target_ind]
	target_file_name = ROI_folder %(temp_folder_OUT) + ROI_name_tmp %(str(seedNo))
	target_file_lines.append(target_file_name)
fobj = open(ROI_folder %(temp_folder_OUT) + target_file,'w')
fobj.writelines(['%s\n' %(line) for line in target_file_lines])
fobj.close()


# 6) Performing Tractography
print 'Performing tractography...'
exclusion_mask = DTI_csf_mask %(temp_folder_OUT)
bedpostX_temp_folder = bedpostX_folder %(temp_folder_IN)
target_text_file = ROI_folder %(temp_folder_OUT) + target_file
domain_file = ROI_folder %(temp_folder_OUT) + 'temp_domain.nii.gz'
for seed_ind in range(len(seed_list)):
	seedNo = seed_list[seed_ind]
	print 'Tractography from seed ' + str(seedNo)
	seed_file = ROI_folder %(temp_folder_OUT) + ROI_name_tmp %(str(seedNo))
	output_folder = tract_out_folder_tmp %(temp_folder_OUT,str(seedNo))
	
	cmd = 'fslmaths ' + ROI_folder %(temp_folder_OUT) + tract_domain + ' -add ' + seed_file + ' -bin ' + ROI_folder %(temp_folder_OUT) + 'temp_domain.nii.gz'
	os.system(cmd) 
	
	cmd_tract = cmd_tract_tmp %(seed_file, exclusion_mask, bedpostX_temp_folder, domain_file, output_folder, target_text_file)
	os.system(cmd_tract)
	
	
# 7) Copying the results
cmd = 'find ' + temp_folder_OUT +'/'+' -type d -exec chmod 777 "{}" \;'
os.system(cmd)
cmd = 'find ' + temp_folder_OUT +'/'+' -type f -exec chmod 666 "{}" \;'
os.system(cmd)

shutil.copytree(temp_result_folder,final_result_folder)
shutil.rmtree(temp_folder)




#cmd = 'cp -R -f ' + temp_folder_OUT + '/* ' + main_folder
#os.system(cmd)
#cmd = 'rm -rf '+ temp_folder_IN
#os.system(cmd)
#cmd = 'rm -rf '+ temp_folder_OUT
#os.system(cmd)