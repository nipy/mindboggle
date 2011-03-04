COREGISTRATION BETWEEN DTI AND T1 IMAGES

NOTE: this is a non-linear registration project between images from the same subject, but not from the same acquisition modality. The correct approach is a linear rigid-body registration and the field map correction should be performed to correct artifacts in EPI images. A non linear registration approach is not the right way to deal the coregistration problem, it is a way to save the data when you don't have the field map.

PRE PROCESSING PIPELINE
These tricks should speed up the computation and improve the results:
1) Compute the coregistration between FA and T1 images
2) Move the lower resolution image (FA) to the higher resolution one (T1) and then compute the inverse transformation
3) Erode the FA map to remove the bright voxels on the edge of the brain (they may mislead the coregistration)
4) Crop the T1 image to reduce the computational time.

COREGISTRATION
See readme file for ants and fsl