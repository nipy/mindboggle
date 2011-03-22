"""
Test of (2-D or 3-D) landmark-based registration in 2D and in 3D.

See correspondence at the bottom regarding parameter settings.

(c) Arno Klein . arno@mindboggle.info . 2011 . MIT license
"""
                                                                               
import os
from subprocess import call
from numpy import random

from measure_overlap import measure_overlap

#
# Inputs
#
verbose = 1
save_results = 0
overwrite = 0
ANTSPATH = '' #os.environ.get("ANTSPATH")
regularizer       = "Gauss"  # Gauss or DMFFD
intensity_measure = "MSQ"    # MSQ or CC
landmark_measure  = "MSQ"    # MSQ or PSE
dim = 3
landmark_dim = 3

#
# Files
#
data_dir = "/hd2/data/Archive/landmark_registration_evaluation_2011_data_output/data/"
output_dir = "/hd2/data/Archive/landmark_registration_evaluation_2011_data_output/output/"
#outpath = output_dir+"test_register_"+str(landmark_dim)+"D_landmark_in_"+str(dim)+"D/"+regularizer+"_"+intensity_measure+"_"+landmark_measure+"/"
#results_dir = "results/test_register_"+str(landmark_dim)+"D_landmark_in_"+str(dim)+"D/"
outpath = output_dir+"test_register_"+str(landmark_dim)+"D_labelshell_in_"+str(dim)+"D/"+regularizer+"_"+intensity_measure+"_"+landmark_measure+"/"
results_dir = "results/test_register_"+str(landmark_dim)+"D_labelshell_in_"+str(dim)+"D/"
temp_dir = output_dir+"temp/"
ext = ".nii.gz"
if dim == 2:
    labels = [21,22,23,24,27,28,41,42,43,44,45,46,47,48,49,50] #labels = [24,28]
    source = data_dir+"S20_to_S05_axial196.nii.gz"
    target = data_dir+"S05_axial196.nii.gz"
    source_labels = data_dir+"S20_to_S05_labels_axial196.nii.gz"
    target_labels = data_dir+"S05_labels_axial196.nii.gz"
    source_landmarks = data_dir+"S20_to_S05_axial196_rightPreCS.nii.gz"
    target_landmarks = data_dir+"S05_axial196_rightPreCS.nii.gz"
elif dim == 3:
    labels = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,41,42,43,44,45,46,
              47,48,49,50,61,62,63,64,65,66,67,68,81,82,83,84,85,86,87,88,
              89,90,91,92,101,102,121,122,161,162,163,164,165,166,181,182]
    #source = data_dir+"S20_to_S05.nii.gz"
    source = data_dir+"S20_to_S05_Gauss_MSQ_0.25_200x200x200_3_0.5.nii.gz"
    target = data_dir+"S05.nii.gz"
    #source_labels = data_dir+"S20_to_S05_labels.nii.gz"
    #source_labels = output_dir+"test_register_3D_labelshell_in_3D/Gauss_MSQ_MSQ/test_0.1_200x200x200_5_0.0_0_0.5_labels.nii.gz"
    source_labels = data_dir+"S20_to_S05_Gauss_MSQ_0.25_200x200x200_3_0.5_labels.nii.gz"
    target_labels = data_dir+"S05_labels.nii.gz"
    if landmark_dim == 2:
        source_landmarks = data_dir+"S20_to_S05_rightPreCS_in_axial196.nii.gz"
        target_landmarks = data_dir+"S05_rightPreCS_in_axial196.nii.gz"
    elif landmark_dim == 3:
        #source_landmarks = data_dir+"S20_to_S05_rightPreCS.nii.gz"
        #target_landmarks = data_dir+"S05_rightPreCS.nii.gz"
        #source_landmarks = data_dir+"S20_to_S05_labelshells/lpba_b_mask_binary.nii.gz"
        #target_landmarks = data_dir+"S05_labelshells/lpba_b_mask_binary.nii.gz"
        #source_landmarks = output_dir+"test_register_3D_labelshell_in_3D/Gauss_MSQ_MSQ/test_0.1_200x200x200_5_0.0_0_0.5_landmarks.nii.gz"
        #target_landmarks = data_dir+"S05_labelshells/lpba_b_mask_binary.nii.gz"
        source_landmarks = data_dir+"S20_to_S05_Gauss_MSQ_0.25_200x200x200_3_0.5_labelshell.nii.gz"
        target_landmarks = data_dir+"S05_labelshells/lpba_b_mask_binary.nii.gz"

#
# Registration parameters
#
gradient_step_sizes = [0.25] #[0.1, 0.15, 0.25, 0.35, 0.5]
if dim == 2:
    iterations = "200x200x200x200"
elif dim == 3:
    iterations = "200x200x201"
options = " --use-Histogram-Matching"
initialize = " --number-of-affine-iterations 0 --continue-affine 0" #10000x10000x10000x10000x10000"

intensity_weights = [0.5] #[0.0, 0.25, 0.5]
landmark_weights = [0.5] #[0.5, 0.75, 1.0]

#
# Initialize settings below
#
intensity_settings = [0]
sigmas = [0]
neighbors = [0]
matching_iters = [0]
#
# Regularization parameters
#
if regularizer == "Gauss":
    regularizer_settings = [3] #[3, 4, 5, 6, 9]  # similarity gradient sigma
    deformation_field_sigma = 0
elif regularizer == "DMFFD":  # DMFFD iterations
    if dim == 2:
        regularizer_settings = ['32x32']
    elif dim == 3:
        regularizer_settings = ['32x32x32']
    bspline_order = 3
#
# Intensity parameters
#
if intensity_measure == "CC":  # radius
    intensity_settings = [5] #[2, 5]
#
# Landmark parameters
#
if landmark_measure == "PSE":
    percent = 1.0 # 0.99 # real number: 1.0 = 100%
    boundary = 0  # 0: not only boundaries
    sigmas = [5] #[5, 100] #[2,100] # big numbers are nearly uniform distributions
    neighbors = [10] #[1, 10] # 1=ICP
    matching_iter = 100000 # partial matching iterations

#
# Begin!
#
print(source, target, source_landmarks, target_landmarks)
if os.path.exists(source) and os.path.exists(target) and \
   os.path.exists(source_landmarks) and os.path.exists(target_landmarks):
    if save_results and dim == 2:
        results_file = results_dir+regularizer+"_"+intensity_measure+"_"+landmark_measure+".csv"
        f_eval = open(results_file, 'w')
    print('Number of tests: ' + str(len(gradient_step_sizes)*len(regularizer_settings)*len(intensity_weights)*len(intensity_settings)*len(landmark_weights)*len(sigmas)*len(neighbors)*len(matching_iters)))
    count = 0
    for gradient_step_size in gradient_step_sizes:
      for regularizer_setting in regularizer_settings:
        for intensity_weight in intensity_weights:
          for intensity_setting in intensity_settings:
            for landmark_weight in landmark_weights:
              for sigma in sigmas:
                for neighbor in neighbors:
                  count += 1
                  if landmark_measure == 'PSE':
                      args0 = [gradient_step_size,iterations,regularizer_setting,intensity_weight,intensity_setting,landmark_weight,sigma,neighbor]
                  elif landmark_measure == 'MSQ':
                      args0 = [gradient_step_size,iterations,regularizer_setting,intensity_weight,intensity_setting,landmark_weight]
                  output = outpath + 'test_' + '_'.join([str(s) for s in args0])
                  output_file = output + ext
                  if save_results and os.path.exists(output_file) and overwrite == 0:
                      raise NameError('File already exists: '+output_file)
                  else:
                    if verbose: print("test " + str(count) + ": " + output_file)
                    if save_results and dim == 3:
                        results_file = results_dir +regularizer+"_"+intensity_measure+"_"+landmark_measure+"_" + '_'.join([str(s) for s in args0])+".csv"
                        f_eval = open(results_file, 'w')

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
                    args1 = " ".join([warp, out, transform, regularize, intensity, landmarks])
                    if verbose: print(args1); print('')
                    p = call(args1, shell="True")

                    args2 = " ".join([apply_warp, source, output_file, '-R ' + target, output+'Warp'+ext, output+'Affine.txt'])
                    if verbose: print(args2); print('')
                    p = call(args2, shell="True")

                    args3 = " ".join([apply_warp, source_landmarks, output+'_landmarks'+ext, '-R ' + target, output+'Warp'+ext, output+'Affine.txt', '--use-NN'])
                    if verbose: print(args3); print('')
                    p = call(args3, shell="True")

                    args4 = " ".join([apply_warp, source_labels, output+'_labels'+ext, '-R ' + target, output+'Warp'+ext, output+'Affine.txt', '--use-NN'])
                    if verbose: print(args4); print('')
                    p = call(args4, shell="True")

                    average_dice, average_jacc = measure_overlap(output+'_labels'+ext, target_labels, labels)
                    print_out = ', '.join(['Test '+str(count), str(average_dice), str(average_jacc), '"'+output_file+'"', '"'+args1+'"\n'])
                    print(print_out)
                    if save_results:
                        f_eval.write(print_out)
                        if dim == 3:
                            f_eval.close()


    if save_results and dim == 2:
        f_eval.close()
else:
    raise NameError('Check input files.')

"""
Parameter ranges suggested by Brian Avants (2/28/2011):
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
    Re PSE sigma: 0 is like a delta function, 100 is like a flat function.
    Re PSE #neighbors: large means you average all neighbor displacements.
                       (#neighbors = 1 is just like icp.)
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
"""