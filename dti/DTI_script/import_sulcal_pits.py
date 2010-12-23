"""
Read the folder list from a text file e for each folder perform the following commands

ImageMath 3 pits_times_250.nii.gz m brain_SulcalPits.nii.gz 250
ImageMath 3 pits_plus_brains.nii.gz + brain.nii.gz pits_times_250.nii.gz
mri_convert pits_plus_brains.nii.gz temp.nii.gz -rl brain_n3.nii.gz -rt nearest -ns 1 -nc
c3d temp.nii.gz -threshold -inf 249 0 1 -o brain_n3_sulcalPits.nii.gz
ImageMath 3 pits_dilate3.nii.gz MD brain_n3_sulcalPits.nii.gz 3

"""

import os, sys, string

fileName = sys.argv[1];

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
	os.system('mkdir -p -m 777 ' +output_folder)
	
	print 'Subject: %s' %(sub)
	print 'Subject: %s' %(input_folder)
	print 'Subject: %s' %(output_folder)

	# Commands execution
	cmd='ImageMath 3 ' + output_folder + 'pits_times_250.nii.gz m ' + input_folder + 'brain_SulcalPits.nii.gz 250'
	os.system(cmd) # Highlight the sulcal pits
	
	cmd='ImageMath 3 ' + output_folder + 'pits_plus_brains.nii.gz + ' + input_folder + 'brain.nii.gz ' + output_folder + 'pits_times_250.nii.gz'
	os.system(cmd) # sum the sulcal pits to the anatomical image
	
	cmd='mri_convert ' + output_folder + 'pits_plus_brains.nii.gz ' + output_folder + 'temp.nii.gz -rl ' + input_folder + 'brain_n3.nii.gz -rt nearest -ns 1 -nc' 
	os.system(cmd) # Reslice the image from freeSurfer space to T1 space
	
	cmd='c3d ' + output_folder + 'temp.nii.gz -threshold -inf 249 0 1 -o ' + output_folder + 'brain_n3_sulcalPits.nii.gz'
	os.system(cmd) # Select the sulcal pits using a threshold
	
	cmd='ImageMath 3 ' + output_folder + 'pits_dilate3.nii.gz MD ' + output_folder + 'brain_n3_sulcalPits.nii.gz 3'
	os.system(cmd) # Dilate the sulcal pits for tractography
	
	# Permissions
	cmd='chmod 666 '+output_folder+'*.*'
	os.system(cmd)
	
fileObj.close()

# python script_arno.py fileName.txt
#
# fileName.txt
# 50475 //data/BI/mindburst/data_pieces/mri_freesurfer_output/2/66/85/