#!/usr/bin/python

"""
Functions for comparing images.


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

"""
def rotate_volume(infile, axis, nrotations, outfile):
    ""
    Rotate image volume.

    Parameters
    ----------
    infile : string
        input file name
    axis : integer
        axis of rotation
    nrotations : integer
        number of times to rotate new xy (ccwise)
    outfile : string
        output file name

    ""
    import numpy as np
    import nibabel as nb

    # Load image
    img = nb.load(infile)
    data = img.get_data()

    # Rotate image
    if nrotations1 > 0:
        data = np.rot90(data, nrotations1)
    if nrotations2 > 0:
        data = np.rot90(data, nrotations2)

    # Save output
    img = nb.Nifti1Image(data, img.get_affine())
    img.to_filename(outfile)
"""

def compare_image_histograms(files):
    """
    Measure similarity between histograms of values from nibabel-readable images.

    The images do not need to be coregistered.

    Parameters
    ----------
    files : list of strings
        file names

    Returns
    -------
    all_bins : numpy array of integers
        histogram bin values
    pairwise_hist_distances : numpy array of floats
        distances between each pair of histograms

    """
    import numpy as np
    import nibabel as nb
    from pylab import plot #, hist

    plot_histograms = False

    # Initialize output
    nbins = 100
    remove_nbins = 1
    all_bins = np.zeros((len(files), nbins - remove_nbins))
    pairwise_hist_distances = np.zeros((len(files), len(files)))

    #---------------------------------------------------------------------------
    # Compute histogram for each image
    #---------------------------------------------------------------------------
    # For each image
    for ifile, file in enumerate(files):

        # Load image
        data = nb.load(file).get_data().ravel()

        # Compute histogram
        bins, bin_edges = np.histogram(data, bins=nbins)
        if remove_nbins > 0:
            bins = bins[remove_nbins::]
        all_bins[ifile, :] = bins

        # Plot histogram:
        if plot_histograms:
            #a,b,c = hist(data, bins=nbins)
            plot(range(len(bins)), bins, '-')

    #---------------------------------------------------------------------------
    # Compute distance between each pair of histograms
    #---------------------------------------------------------------------------
    # Loop through every pair of images
    for ifile1 in range(len(files)):
        for ifile2 in range(len(files)):
            if ifile2 >= ifile1:

                # Store pairwise distances between histogram values
                pairwise_hist_distances[ifile1, ifile2] = np.sqrt(
                    sum((all_bins[ifile1] - all_bins[ifile2])**2)) / len(data)

    return all_bins, pairwise_hist_distances

def register_images(files):
    """
    Register all images to the first image.

    Parameters
    ----------
    files : list of strings
        file names

    """
    import numpy as np
    from nipy.algorithms.registration import HistogramRegistration, resample
    from nipy.utils import example_data
    from nipy import load_image, save_image

    similarity = 'crl1'
    interp = 'pv'
    optimizer = 'powell'

    # Get data
    target_file = files[0]
    J = load_image(target_file)

    # Copy first image as reference
    os.system('mkdir output')
    outimg = os.path.join('output',
             os.path.basename(files[0]) + '0_to_' + os.path.basename(files[0]))
    os.system('cp ' + target_file + ' ' + outimg)

    for isource, source_file in enumerate(files[1::]):
        I = load_image(source_file)

        # Perform affine registration
        # The output is an array-like object such that
        # np.asarray(T) is a customary 4x4 matrix
        R = HistogramRegistration(I, J, similarity=similarity, interp=interp)
        T = R.optimize('affine', optimizer=optimizer)

        # Resample source image
        It = resample(I, T.inv(), reference=J)

        # Save resampled source
        outimg = os.path.join('output',
                 os.path.basename(source_file) + \
                 str(isource + 1) + '_to_' + \
                 os.path.basename(target_file))
        print ('Saving resampled source in: %s' % outimg)
        save_image(It, outimg)

        # Save transformation matrix
        outparams = outimg + '.npy'
        np.save(outparams, np.asarray(T))

def threshold_images(files, threshold_value=0.2):
    """
    Mask images.

    Parameters
    ----------
    files : list of strings
        file names of coregistered files

    """
    import nibabel as nb
    #import nipype.interfaces.mrtrix as mrt

    # Loop through every pair of images
    ref_file = files[0]
    coreg_dir = "output"
    for isource, source_file in enumerate(files):
        infile = os.path.join(coreg_dir,
                 os.path.basename(source_file) + \
                 str(isource) + '_to_' + \
                 os.path.basename(ref_file))
        outfile = os.path.join(coreg_dir,
                 'mask_' + os.path.basename(source_file) + \
                 str(isource) + '_to_' + \
                 os.path.basename(ref_file))

        #thresh = mrt.Threshold()
        #thresh.inputs.in_file = infile
        #thresh.inputs.out_filename = outfile
        #thresh.inputs.absolute_threshold_value = threshold_value
        #thresh.run()
        img = nb.load(infile)
        data = img.get_data()
        data[data < threshold_value] = 0
        data[data >= threshold_value] = 1
        img = nb.Nifti1Image(data, img.get_affine())
        img.to_filename(outfile)

def compute_image_similarities(files, metric='cc'):
    """
    Measure similarity between coregistered nibabel-readable images.

    Parameters
    ----------
    files : list of strings
        file names of coregistered files
    metric: integer
        Cost-function for assessing image similarity. If a string,
        one of 'cc': correlation coefficient, 'cr': correlation
        ratio, 'crl1': L1-norm based correlation ratio, 'mi': mutual
        information, 'nmi': normalized mutual information, 'slr':
        supervised log-likelihood ratio. If a callable, it should
        take a two-dimensional array representing the image joint
        histogram as an input and return a float.

    """
    import os
    import numpy as np
    from nipype.interfaces.nipy.utils import Similarity

    # Initialize output
    pairwise_similarities = np.zeros((len(files), len(files)))

    # Loop through every pair of images
    ref_file = files[0]
    coreg_dir = "output"
    for isource1, source1_file in enumerate(files):
        img1 = os.path.join(coreg_dir,
                 os.path.basename(source1_file) + \
                 str(isource1) + '_to_' + \
                 os.path.basename(ref_file))
        mask1 = os.path.join(coreg_dir,
                 'mask_' + os.path.basename(source1_file) + \
                 str(isource1) + '_to_' + \
                 os.path.basename(ref_file))

        for isource2, source2_file in enumerate(files[isource1+1::]):
            img2 = os.path.join(coreg_dir,
                     os.path.basename(source2_file) + \
                     str(isource2) + '_to_' + \
                     os.path.basename(ref_file))
            mask2 = os.path.join(coreg_dir,
                     'mask_' + os.path.basename(source2_file) + \
                     str(isource2) + '_to_' + \
                     os.path.basename(ref_file))

            # Store pairwise similarities
            similarity = Similarity()
            similarity.inputs.volume1 = img1
            similarity.inputs.volume2 = img2
            similarity.inputs.mask1 = mask1
            similarity.inputs.mask2 = mask2
            similarity.inputs.metric = metric
            res = similarity.run()
            pairwise_similarities[isource1, isource2] = res.outputs.similarity

    return pairwise_similarities

def compute_image_overlaps(files, list_of_labels):
    """
    Measure volume overlaps between coregistered nibabel-readable images.

    Parameters
    ----------
    files : list of strings
        file names of coregistered files
    list_of_labels : list of integers
        label numbers to compute overlap on

    """
    import os
    import numpy as np
    from mindboggle.evaluate.evaluate_labels import measure_volume_overlap

    # Initialize output
    pairwise_overlaps = np.zeros((len(files), len(files),
                                  len(list_of_labels)))

    # Loop through every pair of images
    ref_file = files[0]
    coreg_dir = "output"
    for isource1, source1_file in enumerate(files):
        mask1 = os.path.join(coreg_dir,
                 'mask_' + os.path.basename(source1_file) + \
                 str(isource1) + '_to_' + \
                 os.path.basename(ref_file))

        for isource2, source2_file in enumerate(files[isource1+1::]):
            mask2 = os.path.join(coreg_dir,
                     'mask_' + os.path.basename(source2_file) + \
                     str(isource2) + '_to_' + \
                     os.path.basename(ref_file))

            # Compute and store pairwise overlaps
            overlaps, out_file = measure_volume_overlap(list_of_labels, mask1, mask2)
            pairwise_dice = [x[1] for x in overlaps]
            pairwise_overlaps[isource1, isource2, :] = pairwise_dice

    return pairwise_overlaps


if __name__ == "__main__":

    import os, sys
    import numpy as np
    from mindboggle.evaluate.compare_images import compare_image_histograms, \
        compute_image_similarities, compute_image_overlaps

    do_compare_image_histograms = True
    do_register_images = True
    do_threshold_images = True
    do_compute_image_similarities = True
    do_compute_image_overlaps = True

    dir = '/Users/arno/Dropbox/phantoms/DTI_Phantom'
    files = [os.path.join(dir,'CUDP01CUMR2R1_20120712/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'CUDP01CUMR3R1_20120806/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'CUDP01CUMR4R1_20120915/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'CUDP01CUMR5R1_20121013/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'MGDP01MGMR1R1MGMR1R1_20120801/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'MGDP01MGMR2R1MGMR1R1_20120905/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'MGDP01MGMR3R1MGMR1R1_20121003/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'TXDP01TXMR2R1_20120911/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'TXDP01TXMR3R1_20121002/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'TXDP01TXMR4R1_20121102/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'UMDP01UMMR1R1_20120706/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'UMDP01UMMR2R1_20120803/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'UMDP01UMMR3R1_20120907/DTI_Tensor_FA.nii.gz'),
        os.path.join(dir,'UMDP01UMMR5R1_20121105/DTI_Tensor_FA.nii.gz')]

    dir = '/Users/arno/Dropbox/phantoms/Structural_ADNI_Phantom'
    files = [os.path.join(dir,'CU_SP/CUSP01CUMR2R1_20120711/s003a1001.nii.gz'),
        os.path.join(dir,'CU_SP/CUSP01CUMR3R1_20120806/s003a1001.nii.gz'),
        os.path.join(dir,'CU_SP/CUSP01CUMR4R1_20120915/s003a1001.nii.gz'),
        os.path.join(dir,'CU_SP/CUSP01CUMR5R1_20121013/s003a1001.nii.gz'),
        os.path.join(dir,'MG_SP/MGSP01MGMR1R1_20120627/s021a1001.nii.gz'),
        os.path.join(dir,'MG_SP/MGSP01MGMR1R1MGMR1R1_20120801/s004a1001.nii.gz'),
        os.path.join(dir,'MG_SP/MGSP01MGMR2R1MGMR1R1_20120905/s004a1001.nii.gz'),
        os.path.join(dir,'MG_SP/MGSP01MGMR3R1MGMR1R1_20121003/s004a1001.nii.gz'),
        os.path.join(dir,'TX_SP/TXSP01TXMR2R1_20120907/s401a1004.nii.gz'),
        os.path.join(dir,'TX_SP/TXSP01TXMR3R1_20121002/s601a1006.nii.gz'),
        os.path.join(dir,'TX_SP/TXSP01TXMR4R1_20121108/s501a1005.nii.gz'),
        os.path.join(dir,'UM_SP/UMSP01UMMR1R1_20120706/s1001a1010.nii.gz'),
        os.path.join(dir,'UM_SP/UMSP01UMMR2R1_20120803/s301a1003.nii.gz'),
        os.path.join(dir,'UM_SP/UMSP01UMMR3R1_20120907/s301a1003.nii.gz'),
        os.path.join(dir,'UM_SP/UMSP01UMMR4R1_20121002/s401a1004.nii.gz'),
        os.path.join(dir,'UM_SP/UMSP01UMMR5R1_20121105/s401a1004.nii.gz')]

    if do_compare_image_histograms:
        all_bins, pairwise_hist_distances = compare_image_histograms(files)
        out_file = os.path.join(os.getcwd(), 'pairwise_hist_distances.txt')
        np.savetxt(out_file, pairwise_hist_distances,
                   fmt=len(files) * '%.4f ', delimiter='\t', newline='\n')

    if do_register_images:
        register_images(files)

    if do_threshold_images:
        threshold_images(files, threshold_value=0.2)

    if do_compute_image_similarities:
        pairwise_similarities = compute_image_similarities(files, metric='cc')
        out_file = os.path.join(os.getcwd(), 'output', 'pairwise_similarities.txt')
        np.savetxt(out_file, pairwise_similarities,
                   fmt=len(files) * '%.4f ', delimiter='\t', newline='\n')

    if do_compute_image_overlaps:
        list_of_labels = [1]
        pairwise_overlaps = compute_image_overlaps(files, list_of_labels)
        pairwise_overlap_averages = np.mean(pairwise_overlaps, axis=2)
        out_file = os.path.join(os.getcwd(), 'output', 'pairwise_overlap_averages.txt')
        np.savetxt(out_file, pairwise_overlap_averages,
                   fmt=len(files) * '%.4f ', delimiter='\t', newline='\n')


"""
# Anticipating that there will be a rapidly decreasing distribution
# of low intensity values with a long tail of higher values,
# smooth the bin values (Gaussian), convolve to compute slopes,
# and find the value for the first bin with slope = 0.
from scipy.ndimage.filters import gaussian_filter1d
bins_smooth = gaussian_filter1d(bins.tolist(), 5)
window = [-1, 0, 1]
bin_slopes = np.convolve(bins_smooth, window, mode='same') / (len(window) - 1)
ibin = np.where(bin_slopes == 0)[0]
if len(ibin):
    threshold = bin_edges[ibin[0]]
else:
    threshold = np.median(data)

# Plot histograms:
plot(range(len(bins)), bins, '.', range(len(bins)), bins_smooth,'-')

# Remove first lobe of low values in the data
indices = [i for i,x in enumerate(bins) if x >= threshold]
"""

