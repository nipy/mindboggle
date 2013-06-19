"""
def write_shape_stats(label_files, shape_files,
        area_files='', label_name='', shape_name='',
        exclude_labels=[-1], delimiter=','):
    ""
    Make tables of shape statistics for a given feature (such as label regions).

    Parameters
    ----------
    labels_or_file : list or string
        label number for each vertex or name of VTK file with index scalars
    sulci :  list of integers
        indices to sulci, one per vertex, with -1 indicating no sulcus
    fundi :  list of integers
        indices to fundi, one per vertex, with -1 indicating no fundus
    affine_transform_file : string
        affine transform file to standard space
    transform_format : string
        format for transform file
        Ex: 'txt' for text, 'itk' for ITK, and 'mat' for Matlab format
    area_file :  string
        name of VTK file with surface area scalar values
    exclude_labels : list of lists of integers
        indices to be excluded (in addition to -1)
    delimiter : string
        delimiter between columns, such as ','

    Returns
    -------
    table_file :  string
        output table filename for label shapes

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.utils.io_table import write_shape_stats
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> label_file = 
    >>> shape_file = 
    >>> label_file = 
    >>> label_name = 'labels'
    >>> shape_name = 'thickness'
    >>> exclude_labels = [-1]
    >>> area_file = ''
    >>> delimiter = ','
    >>> #
    >>> write_shape_stats(labels_or_file, sulci, fundi,
    >>>     affine_transform_file, transform_format, area_file,
    >>>     mean_curvature_file, travel_depth_file, geodesic_depth_file,
    >>>     convexity_file, thickness_file, labels_spectra,
    >>>     labels_spectra_norm, labels_spectra_IDs, sulci_spectra,
    >>>     sulci_spectra_norm, sulci_spectra_IDs, exclude_labels, delimiter)

    ""
    import os
    import numpy as np
    from mindboggle.shapes.measure import means_per_label, stats_per_label
    from mindboggle.utils.io_vtk import read_scalars, read_vtk
    from mindboggle.utils.io_table import write_columns

    column_names = []
    first_pass = True
    area_array = []
    if os.path.exists(label_file):
        feature_array, name = read_scalars(label_file, True, True)
    if os.path.exists(shape_file):
        shape_array, name = read_scalars(shape_file, True, True)
    if os.path.exists(label_file):
        label_array, name = read_scalars(label_file, True, True)
    if area_file and os.path.exists(area_file):
        area_array, name = read_scalars(area_file, True, True)

    #---------------------------------------------------------------------
    # Construct a table of average shape values:
    #---------------------------------------------------------------------
    table_file = os.path.join(os.getcwd(), label_name+'_'+shape_name+'.csv')
    table_column_names = []
    columns = []

    #-----------------------------------------------------------------
    # Loop through shape measures:
    #-----------------------------------------------------------------
    table_column_names.extend(column_names[:])
    print('  Compute statistics on {0} {1}'.
         format(label_name, shape_name))

    #-------------------------------------------------------------
    # Mean shapes:
    #-------------------------------------------------------------
    medians, mads, means, sdevs, skews, kurts, \
    lower_quarts, upper_quarts, \
    label_list = stats_per_label(shape_array,
        label_array, exclude_labels, area_array, precision=1)

    # Append shape names and values to columns:
    pr = label_name + ": " + shape_name + ": "
    if np.size(area_array):
        po = " (weighted)"
    else:
        po = ""
    table_column_names.append(pr + 'median' + po)
    """
    table_column_names.append(pr + 'median absolute deviation' + po)
    table_column_names.append(pr + 'mean' + po)
    table_column_names.append(pr + 'standard deviation' + po)
    table_column_names.append(pr + 'skew' + po)
    table_column_names.append(pr + 'kurtosis' + po)
    table_column_names.append(pr + 'lower quartile' + po)
    table_column_names.append(pr + 'upper quartile' + po)
    """
    columns.append(medians)
    """
    columns.append(mads)
    columns.append(means)
    columns.append(sdevs)
    columns.append(skews)
    columns.append(kurts)
    columns.append(lower_quarts)
    columns.append(upper_quarts)
    """

    #-----------------------------------------------------------------
    # Write labels and values to table:
    #-----------------------------------------------------------------
    # Write labels to table:
    write_columns(label_list, 'label', table_file, delimiter)

    # Append columns of shape values to table:
    if columns:
        write_columns(columns, table_column_names, table_file,
                              delimiter, quote=True, input_table=table_file)
    else:
        # Write something to table:
        write_columns([], '', table_file, delimiter)

    return table_file
"""

if __name__== '__main__':

    import os
    import numpy as np
    from mindboggle.shapes.measure import means_per_label, stats_per_label
    from mindboggle.utils.io_vtk import read_scalars, read_vtk
    from mindboggle.utils.io_table import write_columns

    #out_path = '/homedir/Data/Mindboggle-101/'
    #x_path = os.path.join(os.environ['MINDBOGGLE'], 'x')
    ctrl = ['4995_50286','44616_50563','4428009_50287','4429032_50331',
     '4996_50262','44811_50567','4984_50254']
    pre = ['19040_50345', '19056_50501', '4419013_50301', '4934_50201',
    '4974_50244', '19053_50438', '4909_50178', '4970_50237']
    post = ['19040_50368', '19056_50535', '4419013_50322', '4934_50225',
    '4974_50288', '19053_50493', '4909_50184', '4970_50333']
    subjects_list = ctrl
    #subjects_list = read_columns(subjects_list_file, 1)[0]

    exclude_labels = [-1]
    delimiter = ','

    label_name = 'labels'
    shape_name = 'thickness'
    hemi = 'rh'

    #---------------------------------------------------------------------
    # Construct a table of average shape values:
    #---------------------------------------------------------------------
    table_file = os.path.join(os.getcwd(), label_name+'_'+shape_name+'.csv')
    table_column_names = []
    column_names = []
    columns = []

    for subject in subjects_list:
        label_file = '/homedir/Desktop/hippo/workspace/Mindboggle/Labels/_hemi_'+\
                     hemi+'_subject_'+subject+'/DKT_annot_to_VTK/'+hemi+'.DKTatlas40.gcs.vtk'
        shape_file = '/homedir/Desktop/hippo/results/shapes/_hemi_'+hemi+\
                     '_subject_'+subject+'/'+hemi+'.thickness.vtk'
        area_file = '/homedir/Desktop/hippo/results/shapes/_hemi_'+hemi+\
                     '_subject_'+subject+'/'+hemi+'.pial.area.vtk'
    
        area_array = []
        if os.path.exists(label_file):
            feature_array, name = read_scalars(label_file, True, True)
        if os.path.exists(shape_file):
            shape_array, name = read_scalars(shape_file, True, True)
        if os.path.exists(label_file):
            label_array, name = read_scalars(label_file, True, True)
        if area_file and os.path.exists(area_file):
            area_array, name = read_scalars(area_file, True, True)

        #-----------------------------------------------------------------
        # Loop through shape measures:
        #-----------------------------------------------------------------
        table_column_names.extend(column_names[:])
        print('  Compute statistics on {0} {1}'.
             format(label_name, shape_name))
    
        #-------------------------------------------------------------
        # Mean shapes:
        #-------------------------------------------------------------
        medians, mads, means, sdevs, skews, kurts, \
        lower_quarts, upper_quarts, \
        label_list = stats_per_label(shape_array,
            label_array, exclude_labels, area_array, precision=1)
    
        # Append shape names and values to columns:
        pr = label_name + ": " + shape_name + ": "
        if np.size(area_array):
            po = " (weighted)"
        else:
            po = ""
        table_column_names.append(pr + 'median' + po)
        """
        table_column_names.append(pr + 'median absolute deviation' + po)
        table_column_names.append(pr + 'mean' + po)
        table_column_names.append(pr + 'standard deviation' + po)
        table_column_names.append(pr + 'skew' + po)
        table_column_names.append(pr + 'kurtosis' + po)
        table_column_names.append(pr + 'lower quartile' + po)
        table_column_names.append(pr + 'upper quartile' + po)
        """
        columns.append(medians)
        """
        columns.append(mads)
        columns.append(means)
        columns.append(sdevs)
        columns.append(skews)
        columns.append(kurts)
        columns.append(lower_quarts)
        columns.append(upper_quarts)
        """
    
        #-----------------------------------------------------------------
        # Write labels and values to table:
        #-----------------------------------------------------------------
        # Write labels to table:
        write_columns(label_list, 'label', table_file, delimiter)
    
        # Append columns of shape values to table:
        if columns:
            write_columns(columns, table_column_names, table_file,
                                  delimiter, quote=True, input_table=table_file)
        else:
            # Write something to table:
            write_columns([], '', table_file, delimiter)

