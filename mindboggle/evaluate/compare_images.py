#!/usr/bin/python

"""
Functions for comparing images.


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def compute_image_histogram(infile, nbins=100, threshold=0.0):
    """
    Compute histogram values from nibabel-readable image.

    Parameters
    ----------
    infile : string
        input file name
    nbins : integer
        number of bins
    threshold : float
        remove values lower than threshold

    Returns
    -------
    histogram_values : numpy array
        histogram bin values

    Examples
    --------
    >>> import os
    >>> from mindboggle.evaluate.compare_images import compute_image_histogram
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> infile = os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                  'labels', 'labels.manual.nii.gz')
    >>> compute_image_histogram(infile, nbins=100, threshold=0.1)

    """
    import numpy as np
    import nibabel as nb
    #from pylab import plot #, hist

    #---------------------------------------------------------------------------
    # Compute histogram
    #---------------------------------------------------------------------------
    # Load image
    data = nb.load(infile).get_data().ravel()

    # Threshold image
#    if threshold > 0:
#        data = data / max(data)
#        data = data[data >= threshold]

    # Compute histogram
    histogram_values, bin_edges = np.histogram(data, bins=nbins)

    # plot(range(len(histogram_values)), histogram_values, '-')
    ##a,b,c = hist(data, bins=nbins)

    return histogram_values

def compute_image_histograms(infiles, indirectory='', nbins=100, threshold=0.0):
    """
    Compute histogram values from multiple nibabel-readable images.

    Parameters
    ----------
    infiles : list of strings
        input file names
    indirectory : string
        path to input files
    nbins : integer
        number of bins
    threshold : float
        remove values lower than threshold

    Returns
    -------
    histogram_values_list : list of numpy arrays
        histogram bin values for each file

    Examples
    --------
    >>> import os
    >>> from mindboggle.evaluate.compare_images import compute_image_histogram
    >>> data_path = os.environ['MINDBOGGLE_DATA']
    >>> infiles = [os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                    'labels', 'labels.manual.nii.gz'),
    >>>            os.path.join(data_path, 'subjects', 'MMRR-21-1',
    >>>                                    'labels', 'labels.manual.nii.gz')]
    >>> compute_image_histograms(infiles, indirectory='', nbins=100, threshold=0.1)

    """
    import os
    from mindboggle.evaluate.compare_images import compute_image_histogram

    histogram_values_list = []

    #---------------------------------------------------------------------------
    # Compute histograms
    #---------------------------------------------------------------------------
    for infile in infiles:

        infile_path = os.path.join(indirectory, infile)
        histogram_values = compute_image_histogram(infile_path, nbins, threshold)

        histogram_values_list.append(histogram_values)

    return histogram_values_list

def register_images_to_first_image(files, directory=''):
    """
    Compute registration transforms from each image to the first image.

    Parameters
    ----------
    files : list of strings
        input file names
    directory : string
        path to input files

    Returns
    -------
    outfiles : list of strings
        output transform file names

    """
    import os

    run_ants = False
    if run_ants:
        """
        from nipype.interfaces.ants import Registration
        reg = Registration()
        reg.inputs.transforms = ['Affine']
        reg.inputs.transform_parameters = [(2.0,)]
        reg.inputs.number_of_iterations = [[1500, 200]]
        reg.inputs.dimension = 3
        reg.inputs.write_composite_transform = True
        reg.inputs.collapse_output_transforms = False
        reg.inputs.metric = ['Mattes']
        reg.inputs.metric_weight = [1] # Default (value ignored currently by ANTs)
        reg.inputs.radius_or_number_of_bins = [32]
        reg.inputs.sampling_strategy = ['Random']
        reg.inputs.sampling_percentage = [0.05]
        reg.inputs.convergence_threshold = [1.e-8]
        reg.inputs.convergence_window_size = [20]
        reg.inputs.smoothing_sigmas = [[1,0]]
        reg.inputs.shrink_factors = [[2,1]]
        reg.inputs.use_estimate_learning_rate_once = [True]
        reg.inputs.use_histogram_matching = [True] # This is the default
        """
    """
    else:
        from nipype.interfaces.freesurfer import RobustRegister
        reg = RobustRegister()
        reg.inputs.auto_sens = True
        reg.inputs.init_orient = True
    """

    # Get data
    target_file = os.path.join(directory, files[0])

    outfiles = []
    for isource, source_filename in enumerate(files):  #(files[1::]):

        source_file = os.path.join(directory, source_filename)

        # Save transformation matrix
        prefix = 'registered' + str(isource) + '_'
        out_prefix = os.path.join(os.getcwd(), prefix)
        outfile = out_prefix + 'Affine.txt'
        outfiles.append(outfile)
        print('Save registration transform: {0}'.format(outfile))
        print(target_file,source_file,outfile)

        if run_ants:
            """
            reg.inputs.fixed_image = [target_file]
            reg.inputs.moving_image = [source_file]
            #reg.inputs.output_transform_prefix = prefix
            #reg.inputs.initial_moving_transform = outxfm
            reg.inputs.output_warped_image = outfile
            reg.run()
            """
            iterations = '10000x10000x10000x10000x10000'
            cmd = ' '.join(['ANTS 3',
                '-m  MI[' + target_file + ',' + source_file + ',1,32]',
                '-o', out_prefix, '-i 0 --use-Histogram-Matching',
                '--number-of-affine-iterations', iterations])
            print(cmd)
            os.system(cmd)
        else:
            cmd = ' '.join(['flirt -in', source_file,
                            '-ref', target_file,
                            '-dof 7',
                            '-omat', outfile])
            print(cmd)
            os.system(cmd)
            """
            reg.inputs.source_file = source_file
            reg.inputs.target_file = target_file
            #reg.inputs.out_reg_file = outxfm
            reg.inputs.registered_file = outfile
            reg.run()
            """

    return outfiles

def apply_transforms(image_files, transform_files, directory=''):
    """
    Apply transforms to register all images to the first image.

    Parameters
    ----------
    image_files : list of strings
        input image file names
    transform_files : list of strings
        input transform file names
    directory : string
        path to image files

    Returns
    -------
    outfiles : list of strings
        output registered image file names

    """
    import os

    run_ants = False

    # Get data
    target_file = os.path.join(directory, image_files[0])

    outfiles = []
    for isource, source_filename in enumerate(image_files):  #(image_files[1::]):

        source_file = os.path.join(directory, source_filename)
        transform_file = transform_files[isource]

        # Save registered image
        prefix = 'registered' + str(isource) + '_'
        outfile = os.path.join(os.getcwd(),
                  prefix + os.path.basename(image_files[isource]))
        outfiles.append(outfile)
        print('Save registered image: {0}'.format(outfile))

        if run_ants:
            cmd = ' '.join(['WarpImageMultiTransform 3',
                source_file, outfile, '-R', target_file, transform_file])
        else:
            cmd = ' '.join(['flirt -in', source_file,
                            '-ref', target_file,
                            '-applyxfm -init', transform_file,
                            '-out', outfile])
        print(cmd)
        os.system(cmd)

    return outfiles

def threshold_images(files, threshold_value=0.1, save_files=False):
    """
    Threshold images.

    Parameters
    ----------
    files : list of strings
        file names of coregistered files
    threshold_value : float
        threshold value
    save_files : Boolean
        save files?

    Returns
    -------
    outfiles : list of strings
        output file names

    """
    import os
    import nibabel as nb
    #import nipype.interfaces.mrtrix as mrt

    # Loop through every pair of images
    ref_file = files[0]
    coreg_dir = "output"
    outfiles = []
    for ifile, file in enumerate(files):

        #thresh = mrt.Threshold()
        #thresh.inputs.in_file = infile
        #thresh.inputs.out_filename = outfile
        #thresh.inputs.absolute_threshold_value = threshold_value
        #thresh.run()
        img = nb.load(file)
        data = img.get_data()
        if threshold_value > 0:
            data = data / max(data.ravel())
            data[data < threshold_value] = 0
            data[data >= threshold_value] = 1
        img = nb.Nifti1Image(data, img.get_affine())

        if save_files:
            outfile = os.path.join(os.getcwd(),
                      'thresholded_' + os.path.basename(files[ifile]))
            outfiles.append(outfile)
            img.to_filename(outfile)

    return outfiles

def compute_image_similarities(files, masks=[], metric='cc', save_file=False):
    """
    Measure similarity between coregistered nibabel-readable images.

    Parameters
    ----------
    files : list of strings
        file names of coregistered files
    masks : list of strings
        file names of thresholded (mask) files
    metric: integer
        Cost-function for assessing image similarity. If a string,
        one of 'cc': correlation coefficient, 'cr': correlation
        ratio, 'crl1': L1-norm based correlation ratio, 'mi': mutual
        information, 'nmi': normalized mutual information, 'slr':
        supervised log-likelihood ratio. If a callable, it should
        take a two-dimensional array representing the image joint
        histogram as an input and return a float.
    save_file : Boolean
        save file?

    Returns
    -------
    pairwise_similarities : array of floats
        pairwise similarity measures
    outfile : string [optional]
        output filename

    Examples
    --------
    >>> # Missing: include image files and masks
    >>> from mindboggle.evaluate.compare_images import compute_image_similarities
    >>> compute_image_similarities([file1,file2], [mask1,mask2], 'cc', False)

    """
    import os
    import numpy as np
    import nibabel as nb
    #from nipype.interfaces.nipy.utils import Similarity

    # Initialize output
    pairwise_similarities = np.zeros((len(files), len(files)))

    # Loop through every pair of images
    for i1, volume1 in enumerate(files):
        volume1 = nb.load(volume1)
        volume1 = volume1.get_data().ravel()
        if len(masks):
            mask1 = masks[i1]
            mask1 = nb.load(mask1)
            mask1 = mask1.get_data().ravel()

        for i2, volume2 in enumerate(files):
            if i2 == i1:
                pairwise_similarities[i1, i2] = 1
            elif i2 > i1:
                volume2 = nb.load(volume2)
                volume2 = volume2.get_data().ravel()
                if len(masks):
                    mask2 = masks[i2]
                    mask2 = nb.load(mask2)
                    mask2 = mask2.get_data().ravel()

                # Store pairwise similarities
                """
                similarity = Similarity()
                similarity.inputs.volume1 = volume1
                similarity.inputs.volume2 = volume2
                similarity.inputs.mask1 = mask1
                similarity.inputs.mask2 = mask2
                similarity.inputs.metric = metric
                res = similarity.run()
                pairwise_similarities[i1, i2] = res.outputs.similarity
                """

                if len(masks):
                    volume1 = mask1 * volume1
                    volume2 = mask2 * volume2
                similarity = np.corrcoef(volume1,volume2)[0,1]
                pairwise_similarities[i1, i2] = similarity

    if save_file:
        outfile = os.path.join(os.getcwd(), 'pairwise_similarities.txt')
        np.savetxt(outfile, pairwise_similarities,
                   fmt=len(files) * '%.4f ', delimiter='\t', newline='\n')
    else:
        outfile = ''

    return pairwise_similarities, outfile

def compute_image_overlaps(files, list_of_labels, save_file=False):
    """
    Measure volume overlaps between coregistered nibabel-readable images.

    Parameters
    ----------
    files : list of strings
        file names of coregistered files
    list_of_labels : list of integers
        label numbers to compute overlap on
    save_file : Boolean
        save file?

    Returns
    -------
    pairwise_overlaps : array of floats
        pairwise overlap measures
    outfile : string [optional]
        output filename

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
    for ifile1, file1 in enumerate(files):
        for ifile2, file2 in enumerate(files):
            if ifile1 == ifile2:
                pairwise_overlaps[ifile1, ifile2, :] = 1
            if ifile2 > ifile1:

                # Compute and store pairwise overlaps
                overlaps, out_file = measure_volume_overlap(list_of_labels, file1, file2)
                pairwise_dice = [x[1] for x in overlaps]
                pairwise_overlaps[ifile1, ifile2, :] = pairwise_dice

    if save_file:
        #average_pairwise_overlaps = np.mean(pairwise_overlaps, axis=2)
        outfile = os.path.join(os.getcwd(), 'pairwise_overlaps.txt')
        np.savetxt(outfile, pairwise_overlaps,
                   fmt=len(files) * '%.4f ', delimiter='\t', newline='\n')
    else:
        outfile = ''

    return pairwise_overlaps, outfile
