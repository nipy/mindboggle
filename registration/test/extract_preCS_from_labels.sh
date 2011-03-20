datapath=/hd2/data/Archive/landmark_registration_evaluation_2011_data_output/data/


c3d ${datapath}S20_to_S05_labels.nii.gz -threshold 28 28 1 0 -o ${datapath}S20_to_S05_rightPreCG.nii.gz 
c3d ${datapath}S20_to_S05_rightPreCG.nii.gz -dilate 1 1x1x1mm -o ${datapath}S20_to_S05_rightPreCG_dilate1x1x1mm.nii.gz 

c3d ${datapath}S20_to_S05_labels.nii.gz -threshold 22 26 1 0 -o ${datapath}S20_to_S05_rightFG.nii.gz 
c3d ${datapath}S20_to_S05_rightFG.nii.gz -replace 23 0 25 0 -o ${datapath}S20_to_S05_rightFG.nii.gz 
c3d ${datapath}S20_to_S05_rightFG.nii.gz -dilate 1 1x1x1mm -o ${datapath}S20_to_S05_rightFG_dilate1x1x1mm.nii.gz 

c3d ${datapath}S20_to_S05_rightPreCG_dilate1x1x1mm.nii.gz ${datapath}S20_to_S05_rightFG_dilate1x1x1mm.nii.gz -multiply -o ${datapath}S20_to_S05_rightPreCS.nii.gz


c3d ${datapath}S05_labels.nii.gz -threshold 28 28 1 0 -o ${datapath}S05_rightPreCG.nii.gz 
c3d ${datapath}S05_rightPreCG.nii.gz -dilate 1 1x1x1mm -o ${datapath}S05_rightPreCG_dilate1x1x1mm.nii.gz 

c3d ${datapath}S05_labels.nii.gz -threshold 22 26 1 0 -o ${datapath}S05_rightFG.nii.gz 
c3d ${datapath}S05_rightFG.nii.gz -replace 23 0 25 0 -o ${datapath}S05_rightFG.nii.gz 
c3d ${datapath}S05_rightFG.nii.gz -dilate 1 1x1x1mm -o ${datapath}S05_rightFG_dilate1x1x1mm.nii.gz 

c3d ${datapath}S05_rightPreCG_dilate1x1x1mm.nii.gz ${datapath}S05_rightFG_dilate1x1x1mm.nii.gz -multiply -o ${datapath}S05_rightPreCS.nii.gz


fslview ${datapath}S20_to_S05_rightPreCS.nii.gz ${datapath}S05_rightPreCS.nii.gz ${datapath}S20_to_S05_rightPreCS_in_axial196.nii.gz ${datapath}S05_rightPreCS_in_axial196.nii.gz &
