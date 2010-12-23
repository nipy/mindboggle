import nibabel as nb
import numpy as np
import os,sys
import seed_lib as sl

#main_folder = '/Users/denis/Desktop/DTI_tractography/'
main_folder = sys.argv[1]
seed_list = [1014, 2014, 16, 26, 58]
target_list =[]
csf_threshold = 0.8
if (main_folder[-1] != os.sep):
	main_folder=main_folder + os.sep

# Main variables
# Input folders and files
labels_DTI_file = '%sROImaps/mri___aparc+aseg_DTI.nii.gz' %(main_folder)
CBFmap_file = '%sfsl_segmentation/fsl_seg_prob_0.nii.gz' %(main_folder)
Ants_affine = '%soriginal_data/DTI2T1_ants_T_Affine.txt' %(main_folder) 
Ants_invWarp = '%soriginal_data/DTI2T1_ants_T_InverseWarp.nii.gz' %(main_folder)
Ants_Warp = '%soriginal_data/DTI2T1_ants_T_Warp.nii.gz' %(main_folder)
bedpostX_folder = '%sfsl_data.bedpostX' %(main_folder)
T1_csf_file = '%sfsl_segmentation/fsl_seg_prob_0.nii.gz' %(main_folder)
DTI_first_volume = '%soriginal_data/DTI_firstVolume.nii.gz' %(main_folder)

# Output folders and files
label_WM_file = 'corticalROI2WM.nii.gz'
seed_folder = '%stract_doreen/ROIs/' %(main_folder)
ROI_name_tmp = 'ROI_%s.nii.gz'
tract_out_folder_tmp = '%stract_doreen/tract_from_%s' 
DTI_csf_prob = '%sfsl_segmentation/DTI_csf_prob.nii.gz' %(main_folder)
DTI_csf_mask = '%sfsl_segmentation/DTI_csf_mask80.nii.gz' %(main_folder)
csf_mask = '%scsf_80.nii.gz' %(seed_folder)
target_file = '%stargets.txt' %(seed_folder)
temp_folder = '/tmp/tract'

# Command templates
cmd_tract_tmp = 'probtrackx --mode=seedmask -x %s  -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --stop=%s --forcedir --opd -s %s/merged -m %s/nodif_brain_mask  --dir=%s --targetmasks=%s --os2t' #Order (seed_file, stop_mask, bedpostX_folder, bedpostX_folder, output_folder, targets_file)
cmd_move_csf1_tmp = 'WarpImageMultiTransform 3 %s %s -R %s -i  %s %s --use-NN'  # Order (T1_csf_file, out_file, DTI_first_volume, affineFile, inverseWarpFile)
cmd_move_csf2_tmp = 'fslmaths %s -thr %s -bin %s' # Order (DTI_csf_prob, csf_threshold, DTI_csf_mask)


# Making the folders
print 'Making folders...'
cmd = 'mkdir -p -m 777 ' + seed_folder
os.system(cmd)

# Moving the labels from the cortex to the WM
print 'Preparing seed and target ROIs...'
sl.move_labelFile2wm(labels_DTI_file,seed_folder+label_WM_file)

# Extracting the ROI
sl.extract_seed_from_file(seed_folder+label_WM_file,seed_folder+'ROI',seed_list+target_list)

# Moving CSF mask
cmd_move_csf1 = cmd_move_csf1_tmp %(T1_csf_file, DTI_csf_prob, DTI_first_volume, Ants_affine, Ants_invWarp)
os.system(cmd_move_csf1)
cmd_move_csf2 = cmd_move_csf2_tmp %(DTI_csf_prob, str(csf_threshold), DTI_csf_mask)
os.system(cmd_move_csf2)
cmd = 'cp '+DTI_csf_mask+' '+csf_mask
os.system(cmd)

# Making targets file
target_file_lines = []
for seed_ind in range(len(seed_list)):
	seedNo = seed_list[seed_ind]
	seed_file = seed_folder + ROI_name_tmp %(str(seedNo))
	target_file_lines.append(seed_file)
fobj = open(target_file,'w')
fobj.writelines(['%s\n' %(line) for line in target_file_lines])
fobj.close()

# Performing tractography
print 'Performing tractography...'
for seed_ind in range(len(seed_list)):
	seedNo = seed_list[seed_ind]
	print 'Tractography from seed ' + str(seedNo)
	seed_file = seed_folder + ROI_name_tmp %(str(seedNo))
	stop_mask = csf_mask
	output_folder = tract_out_folder_tmp %(main_folder,str(seedNo))
	
	# Preparing temporary output folder
	cmd = 'mkdir -p -m 777 '+temp_folder
	os.system(cmd)
	cmd_tract = cmd_tract_tmp %(seed_file, stop_mask, bedpostX_folder, bedpostX_folder, temp_folder, target_file)
	os.system(cmd_tract)
	
	# Moving the results in the right folder
	cmd = 'mkdir -p -m 777 '+output_folder
	os.system(cmd)
	cmd = 'cp '+temp_folder+'/*.* '+output_folder
	os.system(cmd)
	cmd = 'rm -rf '+temp_folder
	os.system(cmd)