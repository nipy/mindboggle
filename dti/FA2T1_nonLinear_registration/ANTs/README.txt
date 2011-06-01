COREGISTRATION USING ANTS

See the common preprocessing pipeline for registering DTI and T1 images.

ANTS commands:

1) Computing the transformation. It automatically compute the linear transformation before the non linear one and both direct and inverse transformation warp files.
ANTS 3 -m CC[<T1_cropped>,<FA_eroded>,1,2] -r Gauss[3,0] -o <FA2DTI_transf_> -i 100x100x100 -t SyN[0.5] --number-of-affine-iterations 10000x10000x10000x10000x10000 

2) Moving the T1 image to DTI space
WarpImageMultiTransform 3 <T1_cropped> <T1_warped_2_FA> -R <FA_eroded> -i <FA2DTI_transf_Affine.txt> <FA2DTI_transf_InverseWarp>

NOTE: To move T1 ROI maps to DTI space see the python script FA2T1_ROI_coregistration_ANTS
.py. Otherwise it is possible to use the nearest neighbor interpolation:

WarpImageMultiTransform 3 <roi_T1_cropped> <roi_warped_2_FA> -R <FA_eroded> -i <FA2DTI_transf_Affine.txt> <FA2DTI_transf_InverseWarp> --use-NN
