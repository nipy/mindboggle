#!/usr/bin/env python

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
    print(infile)
    data = nb.load(infile).get_data().ravel()

    # Threshold image
    if threshold > 0:
        data = data / max(data)
        data = data[data >= threshold]

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

def register_images_to_ref_images(files, ref_file_index=1, max_angle=90, directory=''):
    """
    Compute registration transforms from each image to its reference image.

    Parameters
    ----------
    files : list of strings
        input image file names
    ref_file_index : integer
        index to reference input image file
    max_angle : integer
        maximum search angle
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

    outfiles = []
    target_file = os.path.join(directory, files[ref_file_index])
    for isource, source_filename in enumerate(files):  #(files[1::]):

        source_file = os.path.join(directory, source_filename)

        # Save transformation matrix
        prefix = 'registered' + str(isource) + '_'
        out_prefix = os.path.join(os.getcwd(), prefix)
        outfile = out_prefix + 'Affine.txt'
        outfiles.append(outfile)
        print('Save registration transform: {0}'.format(outfile))

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
            min_angle = '-' + str(max_angle)
            max_angle = str(max_angle)
            cmd = ' '.join(['flirt -in', source_file,
                            '-ref', target_file,
                            '-dof 7',
                            '-searchrx', min_angle, max_angle,
                            '-searchry', min_angle, max_angle,
                            '-searchrz', min_angle, max_angle,
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

def apply_transforms(files, ref_file_index, transform_files,
                     directory='', interp='trilinear'):
    """
    Apply transforms to register all images to a reference image
    (else the first of the image files).

    Parameters
    ----------
    files : list of strings
        input image file names
    ref_file_index : integer
        index to reference input image file
    transform_files : list of strings
        input transform file names
    directory : string
        path to image files
    interp : string
        interpolation {'trilinear', 'nearestneighbour', 'sinc', 'spline'}

    Returns
    -------
    outfiles : list of strings
        output registered image file names

    """
    import os

    run_ants = False

    outfiles = []
    target_file = os.path.join(directory, files[ref_file_index])
    for isource, source_filename in enumerate(files):  #(files[1::]):

        source_file = os.path.join(directory, source_filename)
        transform_file = transform_files[isource]

        # Save registered image
        prefix = 'registered' + str(isource) + '_'
        outfile = os.path.join(os.getcwd(),
                  prefix + os.path.basename(files[isource]))
        outfiles.append(outfile)
        print('Save registered image: {0}'.format(outfile))

        if run_ants:
            cmd = ' '.join(['WarpImageMultiTransform 3',
                source_file, outfile, '-R', target_file, transform_file])
        else:
            cmd = ' '.join(['flirt -in', source_file,
                            '-ref', target_file,
                            '-applyxfm -init', transform_file,
                            '-interp', interp,
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
            #data = data / max(data.ravel())
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
                overlaps, out_file = measure_volume_overlap(list_of_labels,
                                                            file1, file2)
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

def dot_plot_no_overlaps(table_file, pairs_to_process,
                         category_labels=[], no_ones=False, no_zeros=True,
                         plot_xlabel='', plot_ylabel='', plot_title='',
                         figsize=4, dpi=150):
    """
    Create a dot plot for multiple categories, without overlapping dots.
    Modified after:
    http://stackoverflow.com/questions/8671808/matplotlib-preventing-
                                               overlaying-datapoints

    Parameters
    ----------
    table_file : string
        text file containing space-delimited NxN table of values
    pairs_to_process : list of lists of pairs of integers
        each list if for a category, and each pair contains the indices
        to the (rows/columns of the) table to be compared
    category_labels : list of strings
        names of categories
    no_ones : Boolean
        ignore ones in data?
    no_zeros : Boolean
        ignore zeros in data?

    """
    import numpy as np
    from matplotlib import pyplot as plt
    from itertools import groupby

    # Load file
    data = np.loadtxt(table_file)

    # Collect all within-category data as a separate list
    data_lists = []
    c = 0
    for pairs in pairs_to_process:
        data_list = []
        for pair in pairs:
            data_pair = data[pair[0],pair[1]]
            # Remove ones and zeros from data if specified
            test = True
            if no_ones and data_pair == 1:
                test = False
            if no_zeros and data_pair == 0:
                test = False
            if test:
                data_list.append(data_pair)
        data_lists.append(data_list)

    # Compute means and standard deviations for each list
    means = [np.mean(x) for x in data_lists]
    stds = [np.std(x) for x in data_lists]
    print("Means = {0}".format(means))
    print("SDs = {0}".format(stds))

    # Remove ones and zeros from data if specified
    if no_ones or no_zeros:
        new_data_lists = []
        for data_list in data_lists:
            if no_ones:
                if no_zeros:
                    new_data_list = [x for x in data_list if x !=0 if x != 1]
                else:
                    new_data_list = [x for x in data_list if x != 1]
            elif no_zeros:
                new_data_list = [x for x in data_list if x != 0]
            new_data_lists.append(new_data_list)
        data_lists = new_data_lists


    # Plot dots in vertical lines for each category,
    # offsetting otherwise overlapping dots
    x = []
    y = []
    for indx, klass in enumerate(data_lists):
        klass = groupby(sorted(klass))
        for item, objt in klass:
            objt = list(objt)
            points = len(objt)
            pos = 1 + indx + (1 - points) / 50.
            for item in objt:
                x.append(pos)
                y.append(item)
                pos += 0.04
    plt.figure(figsize=[figsize,figsize], dpi=dpi, facecolor='white', edgecolor=None,
        linewidth=0.0, frameon=True, subplotpars=None, tight_layout=None)

    plt.plot(x, y, 'o')
    plt.xlim((0,len(data_lists) + 1))
    #plt.ylim((0,1))
    if category_labels:
        locs, labels = plt.xticks()
        labels = ['']
        labels.extend(category_labels)
        plt.xticks(locs, labels)
        #plot_file = 'dotplot.png'
        #plt.savefig(plot_file, bbox_inches=0)
    plt.xlabel = plot_xlabel
    plt.ylabel = plot_ylabel
    plt.title = plot_title
    plt.show()

#------------------------------------------------------------------------------
# Example
#------------------------------------------------------------------------------
if __name__ == "__main__":

    import os
    from mindboggle.evaluate.compare_images import dot_plot_no_overlaps

    #-----------------------------------------------------------------------------
    # Data to run
    #-----------------------------------------------------------------------------
    process_phantoms = False  # use phantom (or human) data?
    process_DTI = False  # use diffusion (or structural) data?
    do_similarity = True  # use image similarity (or histogram) table?
    # Choose one:
    pairs_type = 1  # {1,2,3}
                    # 1: different_subject_pairs_per_site
                    # 2: same_subject_pairs_per_site
                    # 3: same_subject_intersite_pairs_per_site
    site_names = ['CU','MG','TX','UM']

    if process_phantoms:
        if process_DTI:
            base_path = '/desk/workspace_DTI_phantoms/Image_comparison_workflow'
            pairs_file = '/Users/arno/Data/Phantom_DTI/Human_retest_pairs.py'
        else:
            base_path = '/desk/workspace_phantoms/Image_comparison_workflow'
            pairs_file = '/Users/arno/Data/Phantom_STR/Human_retest_pairs.py'
    else:
        if process_DTI:
            base_path = '/desk/workspace_DTI_humans/Image_comparison_workflow'
            pairs_file = '/Users/arno/Data/Human_DTI_Updated20130123/Human_retest_pairs.py'

            pairs_type1 = [[[0,1],[0,2],[0,3],[0,5],[0,7]],
                           [[8,9]],
                           [[10,11],[10,12],[10,13],[10,14],[10,15]],
                           []]
            pairs_type2 = [[],[],[],[]]
            pairs_type3 = [[[3,4],[5,6]],[],[],[]]

        else:
            base_path = '/desk/workspace_humans/Image_comparison_workflow'
            pairs_file = '/Users/arno/Data/Human_STR_Updated20130123/Human_retest_pairs.py'

            pairs_type1 = [[[0,2],[0,4],[0,7],[0,10],[0,13],[0,15],
                            [0,17],[0,19],[0,21]],
                           [[23,24],[23,26],[23,28],[23,30],[23,32]],
                           [[33,35],[33,37],[33,39],[33,41],[33,43],
                            [33,45],[33,48],[33,51],[33,53]],
                           [[55,58],[55,61],[55,64],[55,65],[55,66]]]
            pairs_type2 = [[[0,1],[2,3],[4,5],[7,8],[10,11],[13,14],
                            [15,16],[17,18],[19,20],[21,22]],
                           [[24,25],[26,27],[28,29],[30,31]],
                           [[33,34],[35,36],[37,38],[39,40],[41,42],
                            [43,44],[45,46],[48,49],[51,52],[53,54]],
                           [[55,56],[58,59],[61,62],[66,67]]]
            pairs_type3 = [[[4,6],[7,9],[10,12]],
                [],
                           [[45,47],[48,50]],
                           [[55,57],[58,60],[61,63]]]

    if pairs_type == 1:
        pairs = pairs_type1
    elif pairs_type == 2:
        pairs = pairs_type2
    elif pairs_type == 3:
        pairs = pairs_type3

    if do_similarity:
        table = os.path.join(base_path, 'Similarity', 'pairwise_similarities.txt')
        plot_ylabel = 'Image similarity (correlation coefficient)'
        plot_title = 'Similarities per site'
    else:
        table = os.path.join(base_path, 'Compare_histograms', 'vector_distances.txt')
        plot_ylabel = 'Normalized histogram distance'
        plot_title = 'Histogram distances per site'


    dot_plot_no_overlaps(table, pairs, site_names, no_ones=True, no_zeros=True,
        plot_xlabel='Sites', plot_ylabel=plot_ylabel, plot_title=plot_title,
        figsize=4, dpi=150)

