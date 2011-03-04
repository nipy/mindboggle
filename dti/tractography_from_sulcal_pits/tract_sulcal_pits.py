""" This method performs the tractography between the sulcal pits.

INPUT
- main folder: containing all subfolders and files
(See Folder and File Organization section to find the required configuration for the input files)

OUTPUT
- all tractography maps from each sulcal pit and the number of tracts connecting each couple of sulcal pits.
(See the FSL Classification targets option for probtractx to see the output files)

NOTE: in the setting path section the path for the software is set.
"""

import nibabel as nb
import numpy as np
import os,sys
import shutil
import sulcal_pits_lib as spl
import settings
import string

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

temp_folder_template = '/tmp/tract_SP_new_%s_' %(str(sub_ID))
cont =1
temp_folder_candidate = temp_folder_template + str(cont)
while os.path.exists(temp_folder_candidate):
	cont = cont +1
	temp_folder_candidate = temp_folder_template + str(cont)
temp_folder = temp_folder_candidate + os.sep

# -----------------------------------------------------------------------------------------
#                                                                       CONSTANT DEFINITION
dilatation_radius = 2 # Radius for the sulcal pit dilatition


# -----------------------------------------------------------------------------------------
#                                                                         COMMANDS TEMPLATE

# Performing tractography
cmd_tract_tmp = fsl_path+'probtrackx --mode=seedmask -x %s  -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd -s %s/merged -m %s/nodif_brain_mask.nii.gz  --dir=%s --targetmasks=%s --os2t' #Order (seed_file,  bedpostX, domain, output_folder, target_file)

cmd_dilate_tmp = ants_path + 'ImageMath 3 %s MD %s %s' # Order (OutputImage, InputImage, DilatationRadius)

# -----------------------------------------------------------------------------------------
#                                                            FOLDERS AND FILES ORGANIZATION
# INPUT
# Seed
sulcal_pits_T1 = '%sSulcalPits/brain_n3_sulcalPits.nii.gz' %(main_folder)
dti_first_vol = '%soriginal_data/DTI_firstVolume.nii.gz' %(main_folder)
# Coregistration
Ants_affine = '%soriginal_data/DTI2T1_ants_T_Affine.txt' %(main_folder)
Ants_invWarp = '%soriginal_data/DTI2T1_ants_T_InverseWarp.nii.gz' %(main_folder)
Ants_Warp = '%soriginal_data/DTI2T1_ants_T_Warp.nii.gz' %(main_folder) 
# bedpostX
bedpostX_IN = '%sfsl_data.bedpostX' %(main_folder) 


# OUTPUT
sulcal_pits_DTI = '%sOUT/tract_sulcalpits/ROIs/DTI_sulcal_pits.nii.gz' %(temp_folder)
sulcal_pits_DTI_dilated = '%sOUT/tract_sulcalpits/ROIs/DTI_sulcal_pits_dilated.nii.gz' %(temp_folder)
seed_folder = '%sOUT/tract_sulcalpits/ROIs/' %(temp_folder)
seed_file_tmp = temp_folder + 'OUT/tract_sulcalpits/ROIs/seed_%s_mask.nii.gz'
bedpostX_OUT = '%sIN/bedpostX' %(temp_folder)
tract_out_folder_tmp = temp_folder + 'OUT/tract_sulcalpits/tract_from_seed_%s/'
temp_result_folder = '%sOUT/tract_sulcalpits/' %(temp_folder)
target_file = '%sOUT/tract_sulcalpits/ROIs/targets.txt' %(temp_folder)

final_result_template = '%stract_sulcalpits_' %(main_folder)



# -----------------------------------------------------------------------------------------
#                                                                                 MAIN CODE

cmd = 'mkdir -m 777 -p '+temp_folder
os.system(cmd)
cmd = 'mkdir -m 777 -p '+seed_folder	
os.system(cmd)

# Move the sulcal pits from T1 to DTI space, 
spl.move_sulcalPits_from_T1_to_DTI(sulcal_pits_T1,seed_folder,dti_first_vol,Ants_invWarp,Ants_affine,dilatation_radius)

# Copying bedpostX data in the temporary folder to perform tractography
cmd = 'mkdir -m 777 -p ' + bedpostX_OUT
os.system(cmd)
cmd = 'cp ' + bedpostX_IN + '/*.* ' + bedpostX_OUT
os.system(cmd)

# Generating target text file and counting the number of the seeds
target_file_lines = []
exit_flag = False
counter = 1
while not(exit_flag):
	file_name = seed_file_tmp %(str(counter))
	if os.path.isfile(file_name): # file_name exists
		target_file_lines.append(file_name)
		counter = counter +1
	else: # file_name doesn't exist
		exit_flag = True	
fobj = open(target_file,'w')
fobj.writelines(['%s\n' %(line) for line in target_file_lines])
fobj.close()
nSeed = counter -1

# Performing tractography
for counter in range(nSeed):
	seedNo = counter +1
	seed_file = seed_file_tmp %(str(seedNo))
	output_folder = tract_out_folder_tmp %(str(seedNo))
	cmd_tract = cmd_tract_tmp %(seed_file, bedpostX_OUT, bedpostX_OUT, output_folder, target_file)
	os.system(cmd_tract)

# Moving the results to the right place and deleting temporary files
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