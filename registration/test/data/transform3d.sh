
ANTS 3 -m CC[S05.nii.gz, S20.nii.gz, 1, 2] -o S20_to_S05.nii.gz -r Gauss[2,0] -t SyN[0.5] -i 30x100x10 --use-Histogram-Matching --number-of-affine-iterations 10000x10000x10000x10000x10000

WarpImageMultiTransform 3 S20.nii.gz S20_to_S05.nii.gz -R S05.nii.gz S20_to_S05Warp.nii.gz S20_to_S05Affine.txt 

WarpImageMultiTransform 3 S20_labels.nii.gz S20_to_S05_labels.nii.gz -R S05.nii.gz S20_to_S05Warp.nii.gz S20_to_S05Affine.txt --use-NN


