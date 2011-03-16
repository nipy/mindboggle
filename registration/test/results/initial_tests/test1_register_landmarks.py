"""
Test of landmark-based registration on 2D slices.

NOTE: This is the first test to see what range of parameters
successfully register two images together based on a manual landmark.
The # DMFFD iterations was accidentally set to ['6x12'] instead of ['12x6'].

parameter ranges suggested by Brian Avants (2/28/2011):
"
--- Gaussian:  [3,0], [4,0], [5,0], [6,0], [9,0]
--- DMFFD <see below>
--- gradient step size: [0.1, 0.15, 0.2, 0.25, 0.35, 0.5]
--- MSQ: weight relative to intensity should be studied.
--- PSE is problem-dependent.  
    need to experiment and find out useful parameter range.  
    weight relative to intensity should be studied.
       PSE/point-set-expectation/PointSetExpectation[
       fixedImage,movingImage,
       fixedPoints,movingPoints,weight,
       pointSetPercentage,pointSetSigma,
       boundaryPointsOnly,kNeighborhood,
       PartialMatchingIterations=100000]
    the only value i'd suggest fixing is boundaryPointsOnly=0 and
    PartialMatchingIterations=100000 .... unless you want to try
    asymmetric feature-based matching.
"

Nicholas Tustison, re: DMFFD options (3/10/2011):
"
-r DMFFD[option1,option2,option3]

option1:  mesh size for regularization of the gradient field.  For example,
a mesh size of 10x10x5 would specify a mesh size with that resolution
encompassing the entire domain of the 3-D image where the first number is
the number of mesh elements along the first dimension of the image, etc.
The mesh resolution doubles at every image resolution level starting from the
mesh size specified by the user.  So, yes, it is doing a multi-resolution fit.
Re: resolution:
It just depends on what you're analyzing.  There's negligible difference in 
processing time between a mesh size of 1x1x1 and 100x100x100 since
the underlying B-spline fitting algorithm is O(n) where n is the number of
points to which you're trying to fit the B-spline object and the mesh size
doesn't change that.
As regards brain segmentation, I've always followed Brian's lead
and used the standard gaussian settings that Brian uses since 
they've pretty much always worked.  And I really haven't used 
landmarks in conjunction with intensity-based metrics in brain 
registration so it's hard to say.  However, if I remember correctly,
and you might know better than I do, Rueckert uses something like
15x15x15 control points for a typical brain registration and I would 
use the same settings, i.e. a mesh size of 12x12x12 which would
double at every resolution level.  With that in mind, though, I would 
try to adjust the mesh size such that you have an isotropic physical 
mesh element spacing over your entire image.  So, if your image
size is something like 256 mm x 256 mm x 128 mm, I would start with
something like:

-r DMFFD[12x12x6,0,3]

option2: mesh size for regularization of the total field.  The same definition
as for option 1 applies although I typically set this value to 0 since the
DMFFD option is such a heavy regularizer to begin with.  I really haven't
explored specifying a large mesh resolution since the Gaussian option
works for those cases anyway.

option3:  Is the order of spline (default = 3).  The spline order indicates the
degree of polynomial and the higher the order, the higher the regularization.
Typically everybody uses 3 although FNIRT has the option for using quadratic
splines (spline order = 2).  Although spline orders higher than 3 are allowed,
I would advise against it as ITK's B-spline kernel function is written such that
only orders 0 through 3 are explicitly written out.  I had to write a generalized
B-spline kernel function which works, but because it's generalized, it's much
slower.

Re: when would you recommend gaussian regularization over dmffd?
Anytime I have to calculate strain or calculate properties of a deformable body 
like the lungs or tendon or cardiac, I prefer DMFFD regularization using the 
Elast transformation.  Or anytime I have to match up landmarks and I don't care
about symmetry. Collaborators at UPenn in the mechanical engineering department 
just put out an abstract where they investigated different ANTs configurations 
for disc mechanics and they preferred DMFFD to the Gaussian regularization.
"

(c) Arno Klein (2011)  .  arno@mindboggle.info  .  http://www.mindboggle.info
"""
                                                                               
import os
from subprocess import call
from numpy import float, isnan

#
# Inputs
#
verbose = 1
ANTSPATH = '' #os.environ.get("ANTSPATH")
source = "data/S20_to_S05_axial196.nii.gz"
target = "data/S05_axial196.nii.gz"
#labels = [23,27]
labels = [21,22,23,24,27,28,41,42,43,44,45,46,47,48,49,50]
source_labels = "data/S20_to_S05_labels_axial196.nii.gz"
target_labels = "data/S05_labels_axial196.nii.gz"
source_landmarks = "data/S20_to_S05_axial196_leftPreCS.nii.gz"
target_landmarks = "data/S05_axial196_leftPreCS.nii.gz"
outpath = "output/"
results = 'results/test1_registration_average_overlaps.csv'
ext = ".nii.gz"

#
# Registration parameters
#
dim = 2
gradient_step_sizes = [0.1, 0.25, 0.5]
iterations = "100x100"
options = " --use-Histogram-Matching"
initialize = ""  # --number-of-affine-iterations 10000x10000x10000x10000x10000"

#
# Regularization parameters
#
regularizer = "DMFFD"
if regularizer == "DMFFD":
    regularizer_settings = ['6x12']  # DMFFD iterations
    bspline_order = 3
elif regularizer == "Gauss":
    regularizer_settings = [2,3,4]  # similarity gradient sigma
    deformation_field_sigma = 0

#
# Intensity parameters
#
intensity_weights = [0, 0.25, 0.5, 0.75, 1.0]
intensity_measure = "CC"
if intensity_measure == "CC":
    intensity_settings = [2,3,5]  # radius
elif intensity_measure == "MSQ":
    intensity_settings = [0]

#
# Landmark parameters
#
landmark_weights = [0.25, 0.50, 0.75, 1.0]
landmark_measure = "PSE"
if landmark_measure == "PSE":
    percent = 0.99  # real number: 0.99 = 100%
    boundary = 0  # 0: not only boundaries
    sigmas = [5,20] # big numbers are nearly uniform distributions
    neighbors = [5,10]
    matching_iters = [0,100000] # partial matching iterations


#
# Begin!
#                 
f_eval = open(results, 'w')

count = 0
for gradient_step_size in gradient_step_sizes:
  for regularizer_setting in regularizer_settings:
    for intensity_weight in intensity_weights:
      for intensity_setting in intensity_settings:
        for landmark_weight in landmark_weights:
          for sigma in sigmas:
            for neighbor in neighbors:
              for matching_iter in matching_iters:
                #if count<1:
                count += 1
                #args0 = ['test'+str(count),gradient_step_size,regularizer_setting,intensity_weight,intensity_setting,landmark_weight]
                args0 = ['test'+str(count),gradient_step_size,regularizer_setting,intensity_weight,intensity_setting,landmark_weight,sigma,neighbor,matching_iter]
                output = outpath + '_'.join([str(s) for s in args0])
                output_file = output + ext
                if verbose: print("test " + str(count) + ": " + output_file)

                #
                # Arguments for ANTS
                #

                # Output:
                out = "-o " + output_file

                # Transform:
                warp = ANTSPATH + "ANTS " + str(dim)
                apply_warp = ANTSPATH + "WarpImageMultiTransform " + str(dim)

                transform = "-t SyN[" + str(gradient_step_size) +"] -i " + \
                            str(iterations) + options + initialize

                # Regularize:
                if regularizer == "Gauss":
                    regularize = "-r Gauss[" + str(regularizer_setting) + ", " + \
                                 str(deformation_field_sigma) + "]"
                else:
                    regularize = "-r DMFFD["+ regularizer_setting + ", 0, " + str(bspline_order) + "]"

                # Intensity similarity:
                intensity = [target, source, intensity_weight, intensity_setting]
                intensity = "-m "+intensity_measure+"[" + ", ".join([str(s) for s in intensity]) + "]"

                # Landmark similarity:
                landmarks = ""
                if landmark_measure == "MSQ":
                    args = [target_landmarks, source_landmarks, landmark_weight, 0]
                    landmarks = " ".join([landmarks, "-m MSQ[" + ", ".join([str(s) for s in args]) + "]"])
                elif landmark_measure == "PSE":
                    args = [target, source, target_landmarks, source_landmarks,
                            landmark_weight, percent, sigma,
                            boundary, neighbor, matching_iter]
                    landmarks = ", ".join(["-m PSE[" + ", ".join([str(s) for s in args]) + "]"])

                #
                # Run commands
                #
                if os.path.exists(source) and os.path.exists(target) and \
                   os.path.exists(source_landmarks) and os.path.exists(target_landmarks):

                    args1 = " ".join([warp, out, transform, regularize, intensity, landmarks])
                    if verbose: print(args1)
                    #p = call(args1, shell="True")

                    args2 = " ".join([apply_warp, source, output_file, '-R ' + target, output+'Warp'+ext, output+'Affine.txt'])
                    if verbose: print(args2)
                    #p = call(args2, shell="True")

                    args3 = " ".join([apply_warp, source_landmarks, output+'_landmarks'+ext, '-R ' + target, output+'Warp'+ext, output+'Affine.txt', '--use-NN'])
                    if verbose: print(args3)
                    #p = call(args3, shell="True")

                    args3 = " ".join([apply_warp, source_labels, output+'_labels'+ext, '-R ' + target, output+'Warp'+ext, output+'Affine.txt', '--use-NN'])
                    if verbose: print(args3)
                    #p = call(args3, shell="True")

                    compare_labels = 1
                    average_overlaps = 1
                    if compare_labels:
                        average_dice = 0
                        average_jaccard = 0
                        for label in labels:
                            args4 = " ".join(['c3d', output+'_labels'+ext, target_labels, '-overlap', str(label), '>output/temp1_overlap.txt'])
                            if verbose: print(args4)
                            p = call(args4, shell="True")
                            f = open('output/temp1_overlap.txt','r')
                            temp = f.read()
                            if average_overlaps:
                                dice = float(temp.split()[-2].split(',')[0])
                                jaccard = float(temp.split()[-1].split(',')[0])
                                if isnan(dice):
                                    dice = 0
                                if isnan(jaccard):
                                    jaccard = 0
                                average_dice += dice
                                average_jaccard += jaccard
                            else:
                                print_out = ', '.join(['test '+str(count), ' '.join(temp.split()[-6:]), '"'+output_file+'"', '"'+args1+'"\n'])
                                print(print_out)
                                f_eval.write(print_out)
                        if average_overlaps:
                            average_dice = average_dice/len(labels)
                            average_jaccard = average_jaccard/len(labels)
                            print_out = ', '.join(['Test '+str(count), str(average_dice), str(average_jaccard), '"'+output_file+'"', '"'+args1+'"\n'])
                            print(print_out)
                            f_eval.write(print_out)
                    else:
                        # MeasureImageSimilarity ImageDimension whichmetric image1.ext image2.ext 
                        #                        {logfile} {outimage.ext}  {target-value}   {epsilon-tolerance}
                        # target-value and epsilon-tolerance set goals for the metric value 
                        # -- if the metric value is within epsilon-tolerance of the target-value, then the test succeeds 
                        # Metric 0 - MeanSquareDifference, 1 - Cross-Correlation, 2-Mutual Information, 3-SMI

                        args4 = " ".join([ANTSPATH+'MeasureImageSimilarity', str(dim), '2', output+'_landmarks'+ext, target, 'output/temp_intensity_similarity.txt'])
                        if verbose: print(args4)
                        p = call(args4, shell="True")
                        f = open('output/temp_intensity_similarity.txt','r')
                        temp = f.read()

                        intensity_similarity = temp.split()[-1]
                        args5 = " ".join([ANTSPATH+'ImageMath', str(dim), 'output/temp_dtransform_warped.nii.gz', 'D', output+'_landmarks'+ext])
                        if verbose: print(args5)
                        p = call(args5, shell="True")

                        args6 = " ".join([ANTSPATH+'ImageMath', str(dim), 'output/temp_dtransform_target.nii.gz', 'D', target_landmarks])
                        if verbose: print(args6)
                        p = call(args6, shell="True")

                        args7 = " ".join([ANTSPATH+'MeasureImageSimilarity', str(dim), '2', 'output/temp_dtransform_warped.nii.gz', 'output/temp_dtransform_target.nii.gz', 'output/temp_landmark_similarity.txt'])
                        if verbose: print(args7)
                        p = call(args7, shell="True")
                        f = open('output/temp_landmark_similarity.txt','r')
                        temp = f.read()
                        landmark_similarity = temp.split()[-1]

                        # Write test output to file
                        f_eval.write(' '.join([str(landmark_similarity), str(intensity_similarity), '"'+output_file+'"', '"'+args1+'"', '\n']))
