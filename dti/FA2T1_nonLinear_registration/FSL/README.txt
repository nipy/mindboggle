COREGISTRATION USING FSL

See the common preprocessing pipeline for registering DTI and T1 images.

FSL commands:

1) Preliminar affine transformation
flirt -ref <T1_cropped> -in <FA_eroded> -omat <FA_2_T1_Affine> -out <flirt_checking>

2) Non linear registration
fnirt --config=FA_2_same_subject_T1.cnf --aff=<FA_2_T1_Affine> --ref=<T1_cropped> --in=<FA_eroded> --cout=<FA_2_T1_warpcoef> --iout=<FA_warped_2_T1>

3) Computing the inverse warp
invwarp --ref=<FA_eroded> --warp=<FA_2_T1_warpcoef> --out=<T1_2_FA_warpcoef>

4) Moving the T1 image to DTI space
applywarp --ref=<FA_eroded> --in=<T1_cropped> --warp=<T1_2_FA_warpcoef> --out=<T1_warped_2_DTI>

NOTE: To move T1 ROI maps to DTI space see the python script FA2T1_ROI_coregistration.py. Otherwise it is possible to use the nearest neighbor interpolation:
applywarp --ref=<FA_eroded> --in=<roi_T1_cropped> --warp=<T1_2_FA_warpcoef> --out=<roi_warped_2_FA> --interp=nn