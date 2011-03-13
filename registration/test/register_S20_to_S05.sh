
ANTS 3 -m CC[data/S05.nii.gz, data/S20.nii.gz, 1, 2] -o data/S20_to_S05.nii.gz -r Gauss[2,0] -t SyN[0.5] -i 30x100x10 --use-Histogram-Matching --number-of-affine-iterations 10000x10000x10000x10000x10000

WarpImageMultiTransform 3 data/S20.nii.gz data/S20_to_S05.nii.gz -R data/S05.nii.gz data/S20_to_S05Warp.nii.gz data/S20_to_S05Affine.txt 

WarpImageMultiTransform 3 data/S20_labels.nii.gz data/S20_to_S05_labels.nii.gz -R data/S05.nii.gz data/S20_to_S05Warp.nii.gz data/S20_to_S05Affine.txt --use-NN


