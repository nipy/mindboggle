""" This method performs the tractography between cortical ROIs from freesurfer parcellation.
It assumes that the 

INPUT
- main folder: containing all subfolders and files
(See Folder and File Organization section to find the required configuration for the input files)

OUTPUT
- all tractography maps from each sulcal pit and the number of tracts connecting each couple of sulcal pits.
(See the FSL Classification targets option for probtractx to see the output files)

NOTE: in the setting path section the path for the software is set.
"""


import os,sys
import seed_lib
import shutil
import settings
import string
from subprocess import call

# -----------------------------------------------------------------------------------------
#                                                                              SETTING PATH
paths = settings.settings()
fsl_path  = paths[0]
ants_path = paths[1]
afni_path = paths[2]
c3d_path  = paths[3]


# -----------------------------------------------------------------------------------------
#                                                                   READING INPUT VARIABLES
# Subject folder, containing all necessary folders and files
main_folder = sys.argv[1]
if (main_folder[-1] != os.sep):
	main_folder=main_folder + os.sep

tree_words = string.split(main_folder,os.sep)
cont = -1
sub_ID = tree_words[cont]
while sub_ID=='':
	cont = cont-1
	sub_ID = tree_words[cont]

temp_folder_template = '/tmp/tract_freesurfer_rois_%s_' %(str(sub_ID))
cont =1
temp_folder_candidate = temp_folder_template + str(cont)
while os.path.exists(temp_folder_candidate):
	cont = cont +1
	temp_folder_candidate = temp_folder_template + str(cont)
temp_folder = temp_folder_candidate + os.sep

# -----------------------------------------------------------------------------------------
#                                                                       CONSTANT DEFINITION
# List of label of ROIs that are not seed, but white matter parcellation.
labels2mask = [2,4,5,7,14,15,24,30,31,41,43,44,46,62,63,72,77,80,85,251,252,253,254,255]
# Freesurfer labels not used as seed:
# 2 - left cerebral white matter
# 4 - left lateral ventricle
# 5 - left inf-lat ventricle
# 7 - left cerebellum white matter
# 14 - 3rd ventricle
# 15 - 4th ventricle
# 24 - CSF
# 30 - left vessel
# 31 - left choroid plexus
# 41 - right cerebral white matter
# 43 - right lateral ventricle
# 44 - right inf-lat-ventricle
# 46 - right cerebellum white matter
# 62 - right vessel
# 63 - right coroid plexus
# 72 - 5th ventricle
# 77 - WM hypointensities
# 80 - non WM-hypointensities
# 85 - Optic chiasm
# 251 - CC posterior
# 252 - CC mid posterior
# 253 - CC central
# 254 - CC mid anterior
# 255 - CC Anterior


# -----------------------------------------------------------------------------------------
#                                                                         COMMANDS TEMPLATE

# Performing tractography
cmd_tract_tmp = fsl_path+'probtrackx --mode=seedmask -x %s  -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd -s %s/merged -m %s  --dir=%s --targetmasks=%s --os2t' #Order (seed_file,  bedpostX, domain/mask, output_folder, target_file)

# Fsl environment
fsl_env = {
"FSLDIR": fsl_path,
"PATH": "/bin:/usr/bin:" + fsl_path + "/bin",
"FSLOUTPUTTYPE":"NIFTI_GZ"
       }

# Performing erosion
cmd_masking_FA_tmp = fsl_path+'fslmaths %s -ero -bin %s' # Order (input FA file, output mask file)

# Masking data
cmd_masking_data_tmp = fsl_path+'fslmaths %s -mul %s %s' # Order (input file, mask file, output file)

# -----------------------------------------------------------------------------------------
#                                                            FOLDERS AND FILES ORGANIZATION
# INPUT
# ROI folder
parcellation_file = '%sparcellation/ROI_aparc+aseg_warped_2_FA.nii.gz' %(main_folder)
# bedpostX
bedpostX_IN = '%sfsl_data.bedpostX' %(main_folder) 
# FA map
FA_input_file = '%sfsl_data.fa_map/%s_FA.nii.gz' %(main_folder,str(sub_ID))

# OUTPUT 
seed_folder = '%sOUT/tract_freesurfer_parcellation/ROIs/' %(temp_folder)
seed_file_tmp = temp_folder + 'OUT/tract_freesurfer_parcellation/ROIs/seed_%s_mask.nii.gz'
seed_map = temp_folder + 'OUT/tract_freesurfer_parcellation/ROIs/seed_map.nii.gz' 
mask_map = temp_folder + 'OUT/tract_freesurfer_parcellation/ROIs/tract_mask.nii.gz' 
bedpostX_OUT = '%sIN/bedpostX' %(temp_folder)
tract_out_folder_tmp = temp_folder + 'OUT/tract_freesurfer_parcellation/tract_from_seed_%s/'
temp_result_folder = '%sOUT/tract_freesurfer_parcellation/' %(temp_folder)
target_file = '%sOUT/tract_freesurfer_parcellation/ROIs/targets.txt' %(temp_folder)
final_result_template = '%stract_freesurfer_parcellation_' %(main_folder)
parcellation_map = '%sOUT/tract_freesurfer_ROIs/ROIs/parcellation_map.nii.gz' %(temp_folder)
dilated_parcellation_map = '%sOUT/tract_freesurfer_ROIs/ROIs/parcellation_dilated.nii.gz' %(temp_folder)

parcellation_DTI = '%sOUT/tract_freesurfer_ROIs/ROIs/DTI_sulcal_pits.nii.gz' %(temp_folder)

# -----------------------------------------------------------------------------------------
#                                                                                 MAIN CODE

cmd = 'mkdir -m 777 -p '+temp_folder
os.system(cmd)
cmd = 'mkdir -m 777 -p '+seed_folder	
os.system(cmd)

#                                                               Preparing seeds and targets
# -----------------------------------------------------------------------------------------
# NOTE: 
# 1) cortical region are not dilated because otherwise some cortical ROIs get almost completely 
# closed by the other ROIs.
# 2) Some ROI in the parcellation are not used as seed (e.g. ventricles, CSF, white matter...)
# 2) FA is eroded to compute the mask to use as domain for the tractography. 
# 3) Cortical ROI maps are masked with the eroded FA to discard voxels outside the brain.

# Masking the parcellation map
seed_lib.masking_map(parcellation_file,seed_map,labels2mask)

# Eroding the FA to get the mask
cmd_mask = cmd_masking_FA_tmp %(FA_input_file,mask_map)
#os.system(cmd_mask)
try:
	retcode = call(cmd_mask, shell=True, env=fsl_env)
	if retcode != 0:
		print("Error: command " + cmd_mask + " failed with code " + str(retcode))
		quit()
except OSError, e:
	print("Error: can't execute command " + cmd_mask + ": " + str(e))
	quit()  


# Masking the seed map
cmd_masking_data=cmd_masking_data_tmp %(seed_map,mask_map,seed_map)
#os.system(cmd_masking_data)
try:
	retcode = call(cmd_masking_data, shell=True, env=fsl_env)
	if retcode != 0:
		print("Error: command " + cmd_masking_data + " failed with code " + str(retcode))
		quit()
except OSError, e:
	print("Error: can't execute command " + cmd_masking_data + ": " + str(e))
	quit()  
		
# Extracting the seed
seed_lib.extract_all_ROIs(seed_map,seed_file_tmp)

#                                                                   Performing tractography
# -----------------------------------------------------------------------------------------
# 1) Copying bedpostX data in the temporary folder to perform tractography
cmd = 'mkdir -m 777 -p ' + bedpostX_OUT
os.system(cmd)
cmd = 'cp ' + bedpostX_IN + '/*.* ' + bedpostX_OUT
os.system(cmd)

# 2) Generating target text file and counting the number of the seeds
target_file_lines = []
seedNo_list = []
for counter in range(3000):
	file_name = seed_file_tmp %(str(counter))
	if os.path.isfile(file_name): # file_name exists
		target_file_lines.append(file_name)
		seedNo_list.append(counter)
fobj = open(target_file,'w')
fobj.writelines(['%s\n' %(line) for line in target_file_lines])
fobj.close()

# 3) Performing tractography
for counter in range(len(seedNo_list)):
	seedNo = seedNo_list[counter]
	seed_file = seed_file_tmp %(str(seedNo))
	output_folder = tract_out_folder_tmp %(str(seedNo))
	cmd_tract = cmd_tract_tmp %(seed_file, bedpostX_OUT, mask_map, output_folder, target_file)
	try:
		retcode = call(cmd_tract, shell=True, env=fsl_env)
		if retcode != 0:
			print("Error: command " + cmd_tract + " failed with code " + str(retcode))
			quit()
	except OSError, e:
		print("Error: can't execute command " + cmd_tract + ": " + str(e))
		quit()  
	#os.system(cmd_tract)

#                        Moving the results to the right place and deleting temporary files
# -----------------------------------------------------------------------------------------
# Changing file and folder permissions
cmd = 'find ' + temp_folder+'OUT/'+' -type d -exec chmod 777 "{}" \;'
os.system(cmd)
cmd = 'find ' + temp_folder+'OUT/'+' -type f -exec chmod 666 "{}" \;'
os.system(cmd)


final_result_template
cont =1
final_result_candidate = final_result_template + str(cont)
while os.path.exists(final_result_candidate):
	cont = cont +1
	final_result_candidate = final_result_template + str(cont)
final_result_folder = final_result_candidate + os.sep

shutil.copytree(temp_result_folder,final_result_folder)
shutil.rmtree(temp_folder)