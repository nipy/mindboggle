# *Reslice a conformed volume in native ("transformed") space:

import os

subject_path = 'ManualSurfandVolLabels/subjects/'
output_path = 'native_label_volumes/'
finalized_list = "finalized.txt"
f=open(finalized_list)
subject_list = f.readlines()

for subject in subject_list:
  subject = subject.strip('\n')
  conformed_volume_mgz = subject_path + subject + '/mri/aparcNMMjt+aseg.mgz'
  native_volume_mgz = subject_path + subject + '/mri/orig/001.mgz'
  transformed_volume_nii = output_path + subject + '_aparcNMMjt+aseg_native.nii.gz'
  #args = ['mri_vol2vol --mov', conformed_volume_mgz, '--targ', native_volume_mgz, '--regheader --o', transformed_volume_nii]
  args = ['mri_convert -rl',native_volume_mgz,'-rt nearest',conformed_volume_mgz,transformed_volume_nii]
  print(" ".join(args)); os.system(" ".join(args));  # p = Popen(args); p.close()

