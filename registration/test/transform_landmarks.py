
import os
from subprocess import call
from numpy import random

ANTSPATH = '' #os.environ.get("ANTSPATH")
dim = 3
data_dir = "/hd2/data/Archive/landmark_registration_evaluation_2011_data_output/data/"
labels = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,41,42,43,44,45,46,
          47,48,49,50,61,62,63,64,65,66,67,68,81,82,83,84,85,86,87,88,
          89,90,91,92,101,102,121,122,161,162,163,164,165,166,181,182]
source_landmarks_path = data_dir+"S20_labelshells/lpba_b_mask"
output_landmarks_path = data_dir+"S20_to_S05_labelshells/lpba_b_mask"
target = data_dir+"S05.nii.gz"
ext = ".nii.gz"
output = data_dir+"S20_to_S05"

for label in labels:
    source_landmarks = source_landmarks_path + '_' + str(label) + ext
    output_landmarks = output_landmarks_path + '_' + str(label) + ext
    apply_warp = ANTSPATH + "WarpImageMultiTransform " + str(dim)
    args = " ".join([apply_warp, source_landmarks, output_landmarks, '-R ' + target, output+'Warp'+ext, output+'Affine.txt', '--use-NN'])
    print(args); print('')
    p = call(args, shell="True")

source_landmarks = source_landmarks_path + ext
output_landmarks = output_landmarks_path + ext
apply_warp = ANTSPATH + "WarpImageMultiTransform " + str(dim)
args = " ".join([apply_warp, source_landmarks, output_landmarks, '-R ' + target, output+'Warp'+ext, output+'Affine.txt', '--use-NN'])
print(args); print('')
p = call(args, shell="True")
