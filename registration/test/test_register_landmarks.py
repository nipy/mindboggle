"""
Test of landmark-based registration

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

#
# Inputs
#
ANTSPATH = '' #os.environ.get("ANTSPATH")
source = "data/s1"
target = "data/s2"
outpath = "output/"
resultpath = "results/"
ext = ".nii.gz"

#
# Registration parameters
#
dim = 2
gradient_step_sizes = [0.1, 0.2, 0.3, 0.4, 0.5]
iterations = '100x100'
options = " --use-Histogram-Matching"
initialize = " --number-of-affine-iterations 10x10x10x10x10" #10000x10000x10000x10000x10000"

#
# Regularization parameters
#
regularizer = "DMFFD"
if regularizer == "DMFFD":
    regularizer_settings = ['12x12']  # DMFFD iterations
    bspline_order = 3
elif regularizer == "Gauss":
    regularizer_settings = [2,3,4]  # similarity gradient sigma
    deformation_field_sigma = 0

#
# Landmark parameters
#
labels = ["marks"]
weights = [0.25,0.50,1.0]
landmark_measure = "MSQ"
if landmark_measure == "PSE":
    percent = 0.99  # real number: 0.99 = 100%
    boundary = 0  # 0: not only boundaries
    sigma = 100 # big numbers are nearly uniform distributions
    neighbors = [5,10,20]
    matching_iters = [0,100,100000] # partial matching iterations


#
# Begin!
#                 
source_file = source + ext
target_file = target + ext
eval_file = resultpath + 'landmark_intensity_similarities.txt'
f_eval = open(eval_file, 'w')

count = 0
for gradient_step_size in gradient_step_sizes:
  for regularizer_setting in regularizer_settings:
    for weight in weights:
         #for neighbor in neighbors:
          #for matching_iter in matching_iters:
              if count<10:
                count += 1
                args0 = ['test'+str(count),gradient_step_size,regularizer_setting,weight]
                #args0 = ['test'+str(count),gradient_step_size,regularizer_setting,weight,neighbor,matching_iter]
                output = outpath + '_'.join([str(s) for s in args0])
                output_file = output + ext
                print("test " + str(count) + ": " + output_file)

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
                intensity = [target_file, source_file, intensity_weight, intensity_setting]
                intensity = "-m "+intensity_measure+"[" + ", ".join([str(s) for s in intensity]) + "]"

                # Landmark similarity:
                landmarks = ""
                for label in labels:
                    if landmark_measure == "MSQ":
                        args = [target+label+ext, source+label+ext, weight, 0]
                        landmarks = " ".join([landmarks, "-m MSQ[" + ", ".join([str(s) for s in args]) + "]"])
                    elif landmark_measure == "PSE":
                        args = [target_file, source_file, target+label+ext, source+label+ext,
                                weight, percent, sigma,
                                boundary, neighbor, matching_iter]
                        landmarks = " ".join([landmarks, "-m PSE[" + ", ".join([str(s) for s in args]) + "]"])

                #
                # Run commands
                #
                if os.path.exists(source_file) and os.path.exists(target_file):
                    args1 = " ".join([warp, out, transform, regularize, intensity, landmarks])
                    print(args1)
                    #p = call(args1, shell="True")
                    
                    args2 = " ".join([apply_warp, source_file, output_file, '-R ' + target_file, output+'Warp'+ext, output+'Affine.txt'])
                    print(args2)
                    #p = call(args2, shell="True")

                    args3 = " ".join([apply_warp, source+label+ext, output+'_'+label+ext, '-R ' + target_file, output+'Warp'+ext, output+'Affine.txt'])
                    print(args3)
                    #p = call(args3, shell="True")

                    args4 = " ".join([ANTSPATH+'MeasureImageSimilarity', str(dim), '2', output+'_'+label+ext, target_file, 'temp/intensity_similarity.txt'])
                    print(args4)
                    #p = call(args4, shell="True")
                    f = open('temp/intensity_similarity.txt','r')
                    temp = f.read()
                    intensity_similarity = temp.split()[-1]

                    # MeasureImageSimilarity ImageDimension whichmetric image1.ext image2.ext 
                    #                        {logfile} {outimage.ext}  {target-value}   {epsilon-tolerance}
                    # target-value and epsilon-tolerance set goals for the metric value 
                    # -- if the metric value is within epsilon-tolerance of the target-value, then the test succeeds 
                    # Metric 0 - MeanSquareDifference, 1 - Cross-Correlation, 2-Mutual Information, 3-SMI
                    
                    args5 = " ".join([ANTSPATH+'ImageMath', str(dim), 'temp/dtransform_warped.nii.gz', 'D', output+'_'+label+ext])
                    print(args5)
                    p = call(args5, shell="True")

                    args6 = " ".join([ANTSPATH+'ImageMath', str(dim), 'temp/dtransform_target.nii.gz', 'D', target+label+ext])
                    print(args6)
                    p = call(args6, shell="True")

                    args7 = " ".join([ANTSPATH+'MeasureImageSimilarity', str(dim), '2', 'temp/dtransform_warped.nii.gz', 'temp/dtransform_target.nii.gz', 'temp/landmark_similarity.txt'])
                    print(args7)
                    p = call(args7, shell="True")
                    f = open('temp/landmark_similarity.txt','r')
                    temp = f.read()
                    landmark_similarity = temp.split()[-1]

                    # Write test output to file
                    f_eval.write(' '.join([str(landmark_similarity), str(intensity_similarity), '"'+output_file+'"', '"'+args1+'"', '\n']))
