#!/usr/bin/python

"""
Compare images.


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def compare_image_histograms(list_of_files):
    """
    Measure similarity between histograms of values from nibabel-readable images.

    Parameters
    ----------
    list_of_files : list of strings
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
    all_bins = np.zeros((len(list_of_files), nbins - remove_nbins))
    pairwise_hist_distances = np.zeros((len(list_of_files), len(list_of_files)))

    #---------------------------------------------------------------------------
    # Compute histogram for each image
    #---------------------------------------------------------------------------
    # For each image
    for ifile, file in enumerate(list_of_files):

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
    for ifile1 in range(len(list_of_files)):
        for ifile2 in range(len(list_of_files)):
            if ifile2 > ifile1:

                # Store pairwise distances between histogram values
                pairwise_hist_distances[ifile1, ifile2] = np.sqrt(
                    sum((all_bins[ifile1] - all_bins[ifile2])**2)) / len(data)

    return all_bins, pairwise_hist_distances

def compute_image_similarities(list_of_files, dim=3, metric=2):
    """
    Measure similarity between nibabel-readable images.

    From ANTS:
    "MeasureImageSimilarity ImageDimension whichmetric image1.ext image2.ext
        {logfile} {outimage.ext}  {target-value}   {epsilon-tolerance}
    outimage (Not Implemented for MI yet)  and logfile are optional
    target-value and epsilon-tolerance set goals for the metric value --
    if the metric value is within epsilon-tolerance of the target-value,
        then the test succeeds"

    Parameters
    ----------
    list_of_files : list of strings
        file names
    dim : integer
        dimensions of image (e.g., 3 for 3-D)
    metric: integer
        0 - MeanSquareDifference
        1 - Cross-Correlation
        2 - Mutual Information
        3 - SMI

    """
    import os
    #import numpy as np

    # Initialize output
#    pairwise_similarities = np.zeros((len(list_of_files), len(list_of_files)))

    # Loop through every pair of images
    for ifile1 in range(len(list_of_files)):
        for ifile2 in range(len(list_of_files)):
            # Assume symmetry and only calculate one direction per pair
            if ifile2 > ifile1:

                # Store pairwise similarities
                cmd = "MeasureImageSimilarity {0} {1} {2}".format(
                      dim, metric, ' '.join(list_of_files))
                print(cmd); os.system(cmd)
#                pairwise_similarities[ifile, ifile2] = similarity

def compute_image_overlaps(list_of_files, list_of_labels):
    """
    Measure volume overlaps between nibabel-readable images.

    Parameters
    ----------
    list_of_files : list of strings
        file names
    list_of_labels : list of integers
        label numbers to compute overlap on

    """
    import numpy as np
    from mindboggle.evaluate.evaluate_labels import measure_volume_overlap

    # Initialize output
    pairwise_overlaps = np.zeros((len(list_of_files), len(list_of_files),
                                  len(list_of_labels)))

    # Loop through every pair of images
    for ifile1 in range(len(list_of_files)):
        for ifile2 in range(len(list_of_files)):
            if ifile2 > ifile1:

                # Compute and store pairwise overlaps
                overlaps, out_file = measure_volume_overlap(list_of_labels,
                    list_of_files[ifile1], list_of_files[ifile2])
                pairwise_dice = [x[1] for x in overlaps]
                pairwise_overlaps[ifile1, ifile2, :] = pairwise_dice

    return pairwise_overlaps


if __name__ == "__main__":

    import os
    import numpy as np
    from mindboggle.evaluate.compare_images import compare_image_histograms, \
        compute_image_similarities, compute_image_overlaps

    do_compare_image_histograms = True
    do_compute_image_similarities = False
    do_compute_image_overlaps = False

    list_of_files = ['/desk/hln1.nii.gz','/desk/hln2.nii.gz']
#    list_of_files = ['/Users/arno/Dropbox/MB/data/subjects/MMRR-21-1/labels/labels.manual.nii.gz',
#                     '/Users/arno/Dropbox/MB/data/subjects/MMRR-21-1/labels/labels.manual.nii.gz']
#    list_of_labels = [1002, 1003, 1005]

    if do_compare_image_histograms:
        all_bins, pairwise_hist_distances = compare_image_histograms(list_of_files)
        out_file = os.path.join(os.getcwd(), 'pairwise_hist_distances.txt')
        np.savetxt(out_file, pairwise_hist_distances,
                   fmt=len(list_of_files) * '%.4f ', delimiter='\t', newline='\n')

    if do_compute_image_similarities:
        compute_image_similarities(list_of_files, metric=2)

    if do_compute_image_overlaps:
        pairwise_overlaps = compute_image_overlaps(list_of_files, list_of_labels)
        pairwise_overlap_averages = np.mean(pairwise_overlaps, axis=2)
        out_file = os.path.join(os.getcwd(), 'pairwise_overlap_averages.txt')
        np.savetxt(out_file, pairwise_overlap_averages,
                   fmt=len(list_of_files) * '%.4f ', delimiter='\t', newline='\n')


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

