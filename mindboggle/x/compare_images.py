#!/usr/bin/env python

"""
Functions for comparing images.


Authors:
    - Arno Klein  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

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
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> infile = os.path.join(path, 'arno', 'mri', 't1weighted.nii.gz')
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

def compute_image_histograms(infiles, nbins=100, threshold=0.0):
    """
    Compute histogram values from multiple nibabel-readable images.

    Parameters
    ----------
    infiles : list of strings
        input file names
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
    >>> from mindboggle.evaluate.compare_images import compute_image_histograms
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> infiles = [os.path.join(path, 'arno', 'mri', 't1weighted.nii.gz'),
    >>>            os.path.join(path, 'arno', 'labels', 'labels.DKT25.manual.nii.gz')]
    >>> compute_image_histograms(infiles, nbins=100, threshold=0.1)

    """
    import os
    from mindboggle.evaluate.compare_images import compute_image_histogram

    histogram_values_list = []

    #---------------------------------------------------------------------------
    # Compute histograms
    #---------------------------------------------------------------------------
    for infile in infiles:

        histogram_values = compute_image_histogram(infile, nbins, threshold)
        histogram_values_list.append(histogram_values)

    return histogram_values_list

def register_images_to_ref_images(files, ref_file_index=1, max_angle=90,
                                  flirt_command='flirt'):
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

    Returns
    -------
    outfiles : list of strings
        output transform file names

    """
    import os

    outfiles = []
    target_file = files[ref_file_index]
    for isource, source_file in enumerate(files):  #(files[1::]):

        # Save transformation matrix
        prefix = 'registered' + str(isource) + '_'
        out_prefix = os.path.join(os.getcwd(), prefix)
        outfile = out_prefix + 'Affine.txt'
        outfiles.append(outfile)
        print('Save registration transform: {0}'.format(outfile))

        min_angle = '-' + str(max_angle)
        max_angle = str(max_angle)
        cmd = ' '.join([flirt_command, '-in', source_file,
                        '-ref', target_file,
                        '-dof 7',
                        '-searchrx', min_angle, max_angle,
                        '-searchry', min_angle, max_angle,
                        '-searchrz', min_angle, max_angle,
                        '-omat', outfile])
        print(cmd)
        os.system(cmd)

    return outfiles

def apply_transforms(files, ref_file_index, transform_files,
                     interp='trilinear', flirt_command='flirt'):
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
    target_file = files[ref_file_index]
    for isource, source_file in enumerate(files):  #(files[1::]):

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
            cmd = ' '.join([flirt_command, '-in', source_file,
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

def compute_image_similarities(files, intersect_masks=False,
                               metric='cc', save_file=False):
    """
    Measure similarity between coregistered nibabel-readable images.

    Parameters
    ----------
    files : list of strings
        file names of coregistered files
    intersect_masks : Boolean
        compute similarity just on intersection?
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
    >>> from mindboggle.evaluate.compare_images import compute_image_similarities
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> file1 = os.path.join(path, 'arno', 'mri', 't1weighted.nii.gz')
    >>> #file2 = os.path.join(path, 'arno', 'mri', 't1weighted.nii.gz')
    >>> file2 = os.path.join(path, 'arno', 'labels', 'labels.DKT25.manual.nii.gz')
    >>> compute_image_similarities([file1,file2], False, 'cc', False)

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
        for i2, volume2 in enumerate(files):
            if i2 == i1:
                pairwise_similarities[i1, i2] = 1
            elif i2 > i1:
                volume2 = nb.load(volume2)
                volume2 = volume2.get_data().ravel()

                if intersect_masks:
                    #volume1[volume1 < threshold_value] = 0
                    #volume2[volume2 < threshold_value] = 0
                    mask = volume1 * volume2
                    mask[mask > 0] = 1
                    volume1 = mask * volume1
                    volume2 = mask * volume2

                similarity = np.corrcoef(volume1,volume2)[0,1]
                pairwise_similarities[i1, i2] = similarity

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
                         yrange=[], figsize=4, dpi=150):
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
        linewidth=0.0, frameon=True, subplotpars=None) #, tight_layout=None)

    plt.plot(x, y, 'o')
    plt.xlim((0,len(data_lists) + 1))
    if yrange:
        plt.ylim((yrange[0],yrange[1]))
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
    # Pairs of indices to file lists:
    # type 1: different_subject_pairs_per_site
    # type 2: same_subject_pairs_per_site
    # type 3: same_subject_intersite_pairs_per_site
    #-----------------------------------------------------------------------------
    def pairs_types_mri_phantoms():

        pairs_type1 = [[],[],[],[]]
        pairs_type2 = [[[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7]],
                       [[8,9],[9,10],[10,11],[11,12],[12,13],[13,14],[14,15]],
                       [[16,17],[17,18],[18,19],[19,20],[20,21],[21,22]],
                       [[23,24],[24,25],[25,26],[26,27],[27,28],[28,29],[29,30]]]
        pairs_type3 = [[],[],[],[]]

        return pairs_type1, pairs_type2, pairs_type3


    def pairs_types_dmri_phantoms():

        pairs_type1 = [[],[],[],[]]
        pairs_type2 = [[[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8]],
                       [[9,10],[10,11],[11,12],[12,13],[13,14],[14,15]],
                       [[16,17],[17,18],[18,19],[19,20],[20,21],[21,22]],
                       [[23,24],[24,25],[25,26],[26,27],[27,28],[28,29]]]
        pairs_type3 = [[],[],[],[]]

        return pairs_type1, pairs_type2, pairs_type3

    def pairs_types_mri():

        pairs_type1 = [[[0,2],[0,4],[0,7],[0,10],[0,13],[0,15],
                        [0,17],[0,19],[0,21]],
                       [[23,25],[23,27],[23,29],[23,31],[23,33]],
                       [[34,36],[34,38],[34,40],[34,42],[34,44],
                        [34,45],[34,47],[34,50],[34,53],[34,55]],
                       [[57,60],[57,63],[57,66],[57,68],[57,69],[57,70]]]
        pairs_type2 = [[[0,1],[2,3],[4,5],[7,8],[10,11],[13,14],
                        [15,16],[17,18],[19,20],[21,22]],
                       [[23,24],[25,26],[27,28],[29,30],[31,32]],
                       [[34,35],[36,37],[38,39],[40,41],[42,43],
                        [44,45],[47,48],[50,51],[53,54],[55,56]],
                       [[57,58],[60,61],[63,64],[70,71]]]
        pairs_type3 = [[[4,6],[7,9],[10,12]],
            [],
                       [[44,46],[47,49],[50,52]],
                       [[57,59],[60,62],[63,65]]]

        return pairs_type1, pairs_type2, pairs_type3


    def pairs_types_dmri():

        pairs_type1 = [[[0,1],[0,2],[0,4],[0,6],[0,8]],
                       [[9,10],[9,11],[9,12],[9,13]],
                       [[14,15],[14,16],[14,17],[14,18],[14,19],
                        [14,21],[14,23]],
                       [[25,27],[25,29]]]
        pairs_type2 = [[],[],[],[]]
        pairs_type3 = [[[2,3],[4,5],[6,7]],[],[[19,20],[21,22],
                                               [23,24]],[[25,26],[27,28],[29,30]]]

        return pairs_type1, pairs_type2, pairs_type3

    #-----------------------------------------------------------------------------
    # Data to run
    #-----------------------------------------------------------------------------
    process_phantoms = True  # use phantom (or human) data?
    process_dmri = False  # use diffusion (or structural) data?
    do_similarity = True  # use image similarity (or histogram) table?
    # Choose one:
    pairs_type = 2  # {1,2,3}
                    # 1: different_subject_pairs_per_site
                    # 2: same_subject_pairs_per_site
                    # 3: same_subject_intersite_pairs_per_site
    site_names = ['CU','MG','TX','UM']
    base_path = '/homedir/Data/EMBARC/Data'

    if process_phantoms:
        if process_dmri:
            base_path = os.path.join(base_path, 'Scratch_dmri/Scratch')
            pairs_type1, pairs_type2, pairs_type3 = pairs_types_dmri_phantoms()
        else:
            base_path = os.path.join(base_path, 'Scratch/Scratch')
            pairs_type1, pairs_type2, pairs_type3 = pairs_types_mri_phantoms()
    else:
        if process_dmri:
            base_path = os.path.join(base_path, 'Scratch_dmri/Scratch2')
            pairs_type1, pairs_type2, pairs_type3 = pairs_types_dmri()
        else:
            base_path = os.path.join(base_path, 'Scratch/Scratch2')
            pairs_type1, pairs_type2, pairs_type3 = pairs_types_mri()
    if pairs_type == 1:
        pairs = pairs_type1
    elif pairs_type == 2:
        pairs = pairs_type2
    elif pairs_type == 3:
        pairs = pairs_type3

    if do_similarity:
        yrange = [0.0,1.0]
        table = os.path.join(base_path, 'Similarity', 'pairwise_similarities.txt')
        plot_ylabel = 'Image similarity (correlation coefficient)'
        plot_title = 'Similarities per site'
    else:
        yrange = [0,0.25]
        table = os.path.join(base_path, 'Compare_histograms', 'vector_distances.txt')
        plot_ylabel = 'Normalized histogram distance'
        plot_title = 'Histogram distances per site'


    dot_plot_no_overlaps(table, pairs, site_names, no_ones=True, no_zeros=True,
        plot_xlabel='Sites', plot_ylabel=plot_ylabel, plot_title=plot_title,
        yrange=yrange, figsize=4, dpi=150)

