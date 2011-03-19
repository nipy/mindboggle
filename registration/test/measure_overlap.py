"""
Measure overlap between individual label regions in a source and target image.

Input:
arg1, arg2:  source, target images, consisting of index-labeled pixels/voxels
arg3:  list of label indices

(c) Arno Klein . arno@mindboggle.info . 2011 . MIT license
"""
def measure_overlap(source, target, labels):

    from os.path import exists
    if exists(source) and exists(target):

        from subprocess import call
        from numpy import float, isnan

        if len(labels) == 0:
            labels = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,41,42,43,44,45,46,
                      47,48,49,50,61,62,63,64,65,66,67,68,81,82,83,84,85,86,87,88,
                      89,90,91,92,101,102,121,122,161,162,163,164,165,166,181,182]
            #labels = [21,22,23,24,27,28,41,42,43,44,45,46,47,48,49,50]
        if len(source) == 0:
            source = "data/S20_to_S05_labels_axial196.nii.gz"
            #source = "output/test2D/Gauss_MSQ_MSQ/test1_0.6_1_0.5_0_0.5_labels.nii.gz"
        if len(target) == 0:
            target = "data/S05_labels_axial196.nii.gz"
            #"data/S20_to_S05_labels_axial196.nii.gz"

        average_dice = 0
        average_jacc = 0
        for label in labels:
            args = " ".join(['c3d', source, target, '-overlap', str(label), '>temp_overlap.txt'])
            p = call(args, shell="True")
            f = open('temp_overlap.txt','r')
            temp = f.read()
            dice = float(temp.split()[-2].split(',')[0])
            jacc = float(temp.split()[-1].split(',')[0])
            print(' '.join(['Label:', str(label), 'Dice:', str(dice), 'Jaccard:', str(jacc)]))
            if isnan(dice):
                dice = 0
            if isnan(jacc):
                jacc = 0
            average_dice += dice
            average_jacc += jacc
        average_dice = average_dice/len(labels)
        average_jacc = average_jacc/len(labels)
        print('Average Dice: ' + str(average_dice))
        print('Average Jacc: ' + str(average_jacc))

        return average_dice, average_jacc

    else:
        raise NameError('Check input file names.')


"""
# MeasureImageSimilarity ImageDimension whichmetric image1.ext image2.ext 
#                        {logfile} {outimage.ext}  {target-value}   {epsilon-tolerance}
# target-value and epsilon-tolerance set goals for the metric value 
# -- if the metric value is within epsilon-tolerance of the target-value, then the test succeeds 
# Metric 0 - MeanSquareDifference, 1 - Cross-Correlation, 2-Mutual Information, 3-SMI

args = " ".join([ANTSPATH+'MeasureImageSimilarity', str(dim), '2', output+'_landmarks'+ext, target, 'output/temp_intensity_similarity.txt'])
if verbose: print(args)
p = call(args, shell="True")
f = open('output/temp_intensity_similarity.txt','r')
temp = f.read()

intensity_similarity = temp.split()[-1]
args = " ".join([ANTSPATH+'ImageMath', str(dim), 'output/temp_dtransform_warped.nii.gz', 'D', output+'_landmarks'+ext])
if verbose: print(args)
p = call(args, shell="True")

args = " ".join([ANTSPATH+'ImageMath', str(dim), 'output/temp_dtransform_target.nii.gz', 'D', target_landmarks])
if verbose: print(args)
p = call(args, shell="True")

args = " ".join([ANTSPATH+'MeasureImageSimilarity', str(dim), '2', 'output/temp_dtransform_warped.nii.gz', 'output/temp_dtransform_target.nii.gz', 'output/temp_landmark_similarity.txt'])
if verbose: print(args)
p = call(args, shell="True")
f = open('output/temp_landmark_similarity.txt','r')
temp = f.read()
landmark_similarity = temp.split()[-1]
"""
