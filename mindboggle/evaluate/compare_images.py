#!/usr/bin/python

"""
Functions for comparing images.


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def rotate_image_volume(infile, axis, nrotations, outfile):
    """
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

    """
    import os
    import numpy as np
    import nibabel as nb

    # Load image
    img = nb.load(infile)
    data = img.get_data()

    # Rotate image
    if nrotations > 0:
        data = np.rot90(data, nrotations)

    # Save output
    img = nb.Nifti1Image(data, img.get_affine())
    img.to_filename(os.path.join(os.getcwd(), outfile))

def compute_image_histogram(infile, nbins=100, remove_first_nelements=0):
    """
    Compute histogram values from nibabel-readable image.

    Parameters
    ----------
    infile : string
        input file name
    nbins : integer
        number of bins

    Returns
    -------
    histogram_values : numpy array
        histogram bin values

    """
    import numpy as np
    import nibabel as nb
    #from pylab import plot #, hist

    #---------------------------------------------------------------------------
    # Compute histogram
    #---------------------------------------------------------------------------
    # Load image
    data = nb.load(infile).get_data().ravel()

    # Compute histogram
    histogram_values, bin_edges = np.histogram(data, bins=nbins)

    # plot(range(len(histogram_values)), histogram_values, '-')
    ##a,b,c = hist(data, bins=nbins)

    return histogram_values[remove_first_nelements::]

def register_images_to_first_image(files):
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
    Threshold images.

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

