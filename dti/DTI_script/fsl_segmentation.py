""" Compute the segmentation by using FSL FAST tool.
1. Compute the transformation between the T1 data and the template to use the priors in the segmentation step
2. Perform the segmentation using the priors
"""

import os,sys

T1_data_TMP = '/mind_cart/DTI/%s/original_data/T1data_masked.nii.gz'
result_folder_TMP = '/mind_cart/DTI/%s/fsl_segmentation/'
T12templ  = 'T12templ_transf.mat'
templ2T1  = 'templ2T1_transf.mat'
T1inTempl = 'T1data2templ.nii.gz'
segm_out  = 'fsl_seg'


move2template='flirt -in %s -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -out %s -omat %s -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear' #[T1data_masked, movedImage, transfFile]

get_inverseT='convert_xfm -omat %s -inverse %s' # [outTransFile inTransFile]

perform_seg='fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -a %s -o %s -p %s' # [transformation, outPath, inData]

subID = sys.argv[1]
T1_data = T1_data_TMP %(subID)
result_folder = result_folder_TMP %(subID)

print 'Subject %s segmentation' %(subID)
print 'Input data: %s' %(T1_data)
print 'Output folder: %s' %(result_folder)

if not(os.path.exists(result_folder)):
	print ' Making the result folder...'
	os.system('mkdir -p -m 777 ' +result_folder)

print ' Computing the template coregistration...'
cmd = move2template %(T1_data, result_folder+T1inTempl, result_folder+T12templ)
os.system(cmd)

cmd = get_inverseT %(result_folder+templ2T1, result_folder+T12templ)
os.system(cmd)

print ' Performing the segmentation...'
cmd = perform_seg %(result_folder+templ2T1, result_folder+segm_out, T1_data)
os.system(cmd)

os.system('chmod 666 '+result_folder+'*.*')
print ' Complete!\n\n'