"""
FSL_lib

list of python functions to preform tractography using fsl

-------------------------------------------------------------
do_track(seed_file,bedpostx_folder,output_folder)
	Perform tractography using fsl function protrackx.
	INPUT
	 * seed_file       :the file containing the seed mask
	 * bedpostx_folder :the folder containing the bedpostx results
	 * output_folder   :the folder where store the results
	OUTPUT
	The output_folder will contain
	 * fdt_paths.nii.gz :a 3D image file containing the output connectivity distribution to the seed mask
	 * fdt.log          :a log of the setup of the FDT GUI when the analysis was run
	 * probtrackx.log   :a text record of the command that was run
	 * waytotal         :a text file containing a single number corresponding to the total number of generated tracts that have not been rejected by inclusionfile_sepexclusion mask criteria
  
"""	
import os,sys,subprocess
file_sep='/'

def do_track(seed_file,bedpostx_folder,output_folder):

	if not(os.path.isfile(seed_file)):
		print 'ERROR: file %s not found' %(seed_file)
		raise
		
	if bedpostx_folder[len(bedpostx_folder)-1] != file_sep:
		bedpostx_folder = bedpostx_folder + file_sep
			
	if output_folder[len(output_folder)-1] == file_sep:
		output_folder = output_folder[:len(output_folder)-1]

	# Definition of the fsl command to perform tractography
	cmd = ['probtrackx']
	
	# Definition of the arguments
	arg = ['--mode=seedmask'] # Declare that we use an ROI maskfs
	arg = arg + ['-x '+seed_file] # Define the seed mask file
	arg = arg + ['-V 0'] # Verbose level [0..2]
	arg = arg + ['-l'] # Perform loopcheck on paths 
	arg = arg + ['-c 0.2'] # Curvature threshold (default=0.2)
	arg = arg + ['-S 2000'] # Number of step per sample (default=2000)
	arg = arg + ['--steplength=0.5'] # Step length (default=0.5)
	arg = arg + ['-P 5000'] # Number of samples drawn per voxel (default=5000)
	arg = arg + ['--forcedir'] # Use the actual directory name given - i.e. don't add + to make a new directory
	arg = arg + ['--opd'] # output path distribution 
	arg = arg + ['-s '+bedpostx_folder+'merged'] # basename for samples files
	arg = arg + ['-m '+bedpostx_folder+'nodif_brain_mask'] # brain mask
	arg = arg + ['--dir='+output_folder] # Output folder
	
	print cmd
	print arg
	
	comando = cmd[0]
	for opt in arg:
		comando = comando+' '+opt
	
	os.system(comando)
	#subprocess.call(cmd+arg)
	
	# Folder and file permission changement
#	cmd = ['chmod']
#	arg = ['777',output_folder]
#	subprocess.call(cmd+arg)
#	arg = ['666',output_folder+file_sep+"*.*"]
#	subprocess.call(cmd+arg)
	cmd = 'chmod 777 ' + output_folder
	os.system(cmd)
	cmd = 'chmod 666 ' + output_folder + file_sep + '*.*'
	os.system(cmd)
	