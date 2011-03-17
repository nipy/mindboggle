
from subprocess import call

#labels = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,41,42,43,44,45,46,47,48,49,50,61,62,63,64,65,66,67,68,81,82,83,84,85,86,87,88,89,90,91,92,101,102,121,122,161,162,163,164,165,166,181,182]
labels = [21,22,23,24,27,28,41,42,43,44,45,46,47,48,49,50]
#source = "output/test2D/Gauss_MSQ_MSQ/test1_0.6_1_0.5_0_0.5_labels.nii.gz"
source = "data/S20_to_S05_labels_axial196.nii.gz"
#"data/S20_to_S05_labels_axial196.nii.gz"
target = "data/S05_labels_axial196.nii.gz"

average_dice = 0
average_jacc = 0
for label in labels:
    args = " ".join(['c3d', source, target, '-overlap', str(label), '>temp_overlap.txt'])
    p = call(args, shell="True")
    f = open('temp_overlap.txt','r')
    temp = f.read()
    print(temp)
    average_dice += float(temp.split()[-2].split(',')[0])
    average_jacc += float(temp.split()[-1].split(',')[0])
average_dice = average_dice/len(labels)
average_jacc = average_jacc/len(labels)
print('Average Dice: ' + str(average_dice))
print('Average Jacc: ' + str(average_jacc))

