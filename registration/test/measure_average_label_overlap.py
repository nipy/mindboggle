
labels = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,41,42,43,44,45,46,47,48,49,50,61,62,63,64,65,66,67,68,81,82,83,84,85,86,87,88,89,90,91,92,101,102,121,122,161,162,163,164,165,166,181,182]
source = "data/S20_to_S05_labels_axial196.nii.gz"
target = "data/S05_labels_axial196.nii.gz"

average_dice = 0
average_jaccard = 0
for label in labels:
    args = " ".join(['c3d', source, target, '-overlap', str(label), '>temp_overlap.txt'])
    p = call(args, shell="True")
    f = open('output/temp_overlap.txt','r')
    temp = f.read()
    average_overlaps = 1
    if average_overlaps:
        average_dice += float(temp.split()[-2].split(',')[0])
        average_jaccard = float(temp.split()[-1].split(',')[0])
    else:
        print_out = ', '.join(['test '+str(count), ' '.join(temp.split()[-6:]), '"'+output_file+'"', '"'+args1+'"\n'])
        print(print_out)
        f_eval.write(print_out)
if average_overlaps:
    average_dice = average_dice/len(labels)
    average_jaccard = average_jaccard/len(labels)
    print(average_dice)
    print(average_jaccard)

