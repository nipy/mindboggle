#!/usr/bin/env python
"""
Functions for writing tables


Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com
    - Forrest Sheng Bao, 2012  (forrest.bao@gmail.com)  http://fsbao.net

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def read_columns(filename, n_columns=1, trail=False):
    """
    Read n-column text file. Assumes space(s) as delimiter.

    Parameters
    ----------
    filename :  name of text file [string]
    n_columns :  number of columns to extract [integer]
    trail :  combine all remaining columns as a string
             in the final list [Boolean]

    Returns
    -------
    columns :  a list of lists of strings, one list per column of text.

    """
    import os
    import re

    columns = [[] for x in range(n_columns)]
    if os.path.exists(filename):
        Fp = open(filename, 'r')
        lines = Fp.readlines()
        for line in lines:
            if line:
                row = re.findall(r'\S+', line)
                if len(row) >= n_columns:
                    for icolumn in range(n_columns):
                        if trail and icolumn == n_columns - 1:
                            columns[icolumn].append(' '.join(row[icolumn::]))
                        else:
                            columns[icolumn].append(row[icolumn])
                else:
                    import os
                    os.error('The number of columns in {0} is less than {1}.'.format(
                             filename, n_columns))
        Fp.close()

    return columns


def write_columns(columns, column_names, delimiter=',', quote=True,
                  input_table='', output_table=''):
    """
    Write table with columns and column names.  Assumes space(s) as delimiter.

    If there is an input table file to append to, assume a one-line header.

    Parameters
    ----------
    columns :  list of lists of floats or integers
        values (each list is a column of values)
    column_names :  list of strings
        names of columns
    delimiter : string
        delimiter between columns, such as ','
    bracket : string
        string bracketing each element, such as '"'
    input_table : string (default is empty string)
        name of table file to which the columns are to be appended
    output_table : string
        name of output table file (full path)

    Returns
    -------
    output_table : string
        name of output table file (full path)

    Examples
    --------
    >>> from mindboggle.utils.io_table import write_columns
    >>> labels = ['category one', 'category two', 'category three', 'category four']
    >>> values = [0.12, 0.36, 0.75, 0.03]
    >>> values2 = [32, 87, 53, 23]
    >>> columns = [labels, values]
    >>> column_names = ['label', 'value']
    >>> delimiter = ','
    >>> quote = True
    >>> input_table = ''
    >>> output_table = 'write_columns.csv'
    >>> output_table = write_columns(columns, column_names, delimiter, quote,
    >>>                              input_table, output_table)
    >>> write_columns(values2, 'value 2', delimiter, quote,
    >>>               input_table=output_table, output_table=output_table)

    """
    import os
    import sys
    from mindboggle.utils.io_table import read_columns

    if not output_table:
        if input_table:
            s = os.path.basename(input_table).split('.')[0] + '.csv'
            output_table = os.path.join(os.getcwd(), s)
        else:
            if column_names:
                s = '_'.join(column_names)
            else:
                s = 'table'
            output_table = os.path.join(os.getcwd(), s + '.csv')

    if quote:
        q = '"'
    else:
        q = ''

    #-----------------------
    # Check format of inputs
    #-----------------------
    # If the list contains integers or floats, put in a list:
    if columns:
        if isinstance(columns[0], int) or isinstance(columns[0], float) or \
           isinstance(columns[0], str):
            columns = [columns]
        # If the list contains all lists, accept format:
        elif all([isinstance(x, list) for x in columns]):
            pass
        else:
            print("Error: columns contains unacceptable elements.")
            print("columns type is: {0}".format(type(columns)))
            print("columns length is: {0}".format(len(columns)))
            print("columns[0] type is: {0}".format(type(columns[0])))
            sys.exit()
        # If column_names is a string, create a list containing
        # as many of this string as there are columns.
        if isinstance(column_names, str):
            column_names = [column_names for x in columns]
        elif isinstance(column_names, list):
            if len(column_names) < len(columns):
                column_names = [column_names[0] for x in columns]
            else:
                pass
        else:
            print("Error: column_names is neither a list nor a string")
            sys.exit()

        #-----------------------------------
        # Read columns from input table file
        #-----------------------------------
        input_names = ''
        input_columns = ['' for x in columns[0]]
        if input_table and os.path.exists(input_table):
            input_columns = read_columns(input_table, n_columns=1, trail=True)
            if input_columns and len(input_columns[0]) > 0:
                input_names = input_columns[0][0]
                input_columns = input_columns[0][1::]

        #--------------
        # Write to file
        #--------------
        Fp = open(output_table, 'wa')
        if column_names:
            column_names = [q+x+q for x in column_names]
            if input_table:
                Fp.write(delimiter.join([input_names,
                         delimiter.join(column_names) + "\n"]))
            else:
                Fp.write(delimiter.join(column_names) + "\n")
        #else:
        #    Fp.write(input_names + "\n")

        for irow in range(len(columns[0])):
            if input_table:
                Fp.write(input_columns[irow] + delimiter)
            for icolumn, column in enumerate(columns):
                if icolumn < len(columns)-1:
                    Fp.write('{0}{1}{2}{3}'.format(
                        q, column[irow], q, delimiter))
                else:
                    Fp.write('{0}{1}{2}'.format(q, column[irow], q))
            Fp.write("\n")

        Fp.close()

        if not os.path.exists(output_table):
            raise(IOError(output_table + " not found"))

    else:
        print("NOTE: 'columns' is empty. Nothing written.")

    return output_table


def write_rows(filename, list_of_lines, header=""):
    """
    Write a list to a file, one line per list element.

    Parameters
    ----------
    filename : string
        name of output file
    list_of_lines :  list
        each element is written to file as a line
    header : string (default is empty string)
        header to write at the top of the file

    Returns
    -------
    filename : string
        name of output file

    """

    Fp = open(filename, 'w')

    if header:
        Fp.write(header + '\n')

    for element in list_of_lines:
        Fp.write(str(element) + '\n')

    Fp.close()

    return filename


def write_shape_stats(labels_or_file=[], sulci=[], fundi=[],
        affine_transform_file='', transform_format='itk',
        area_file='', mean_curvature_file='', travel_depth_file='',
        geodesic_depth_file='', freesurfer_convexity_file='',
        freesurfer_thickness_file='',
        labels_spectra=[], labels_spectra_IDs=[],
        sulci_spectra=[], sulci_spectra_IDs=[],
        labels_zernike=[], labels_zernike_IDs=[],
        sulci_zernike=[], sulci_zernike_IDs=[],
        exclude_labels=[-1], delimiter=','):
    """
    Make tables of shape statistics per label, sulcus, and/or fundus.

    Note ::
        This function is tailored for Mindboggle outputs.

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
    mean_curvature_file :  string
        name of VTK file with mean curvature scalar values
    travel_depth_file :  string
        name of VTK file with travel depth scalar values
    geodesic_depth_file :  string
        name of VTK file with geodesic depth scalar values
    freesurfer_convexity_file :  string
        name of VTK file with FreeSurfer convexity scalar values
    freesurfer_thickness_file :  string
        name of VTK file with FreeSurfer thickness scalar values
    labels_zernike : list of lists of floats
        Laplace-Beltrami spectra for labeled regions
    labels_spectra_IDs : list of integers
        unique ID numbers (labels) for labels_spectra
    sulci_spectra : list of lists of floats
        Laplace-Beltrami spectra for sulci
    sulci_spectra_IDs : list of integers
        unique ID numbers (labels) for sulci_spectra
    labels_zernike : list of lists of floats
        Zernike moments for labeled regions
    labels_zernike_IDs : list of integers
        unique ID numbers (labels) for labels_zernike
    sulci_zernike : list of lists of floats
        Zernike moments for sulci
    sulci_zernike_IDs : list of integers
        unique ID numbers (labels) for sulci_zernike
    exclude_labels : list of lists of integers
        indices to be excluded (in addition to -1)
    delimiter : string
        delimiter between columns, such as ','

    Returns
    -------
    label_table :  string
        output table filename for label shapes
    sulcus_table :  string
        output table filename for sulcus shapes
    fundus_table :  string
        output table filename for fundus shapes

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.utils.io_table import write_shape_stats
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> labels_or_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> fundi_file = os.path.join(path, 'arno', 'features', 'fundi.vtk')
    >>> sulci, name = read_scalars(sulci_file)
    >>> fundi, name = read_scalars(fundi_file)
    >>> affine_transform_file = os.path.join(path, 'arno', 'mri', 't1weighted_brain.MNI152Affine.txt')
    >>> #transform_format = 'mat'
    >>> transform_format = 'itk'
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> mean_curvature_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> travel_depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> geodesic_depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.geodesic_depth.vtk')
    >>> freesurfer_convexity_file = ''
    >>> freesurfer_thickness_file = ''
    >>> delimiter = ','
    >>> #
    >>> labels, name = read_scalars(labels_or_file)
    >>> labels_spectra = []
    >>> labels_spectra_IDs = []
    >>> sulci_spectra = []
    >>> sulci_spectra_IDs = []
    >>> labels_zernike = []
    >>> labels_zernike_IDs = []
    >>> sulci_zernike = []
    >>> sulci_zernike_IDs = []
    >>> exclude_labels = [-1]
    >>> #
    >>> write_shape_stats(labels_or_file, sulci, fundi,
    >>>     affine_transform_file, transform_format, area_file,
    >>>     mean_curvature_file, travel_depth_file, geodesic_depth_file,
    >>>     freesurfer_convexity_file, freesurfer_thickness_file,
    >>>     labels_spectra, labels_spectra_IDs,
    >>>     sulci_spectra, sulci_spectra_IDs,
    >>>     labels_zernike, labels_zernike_IDs,
    >>>     sulci_zernike, sulci_zernike_IDs,
    >>>     exclude_labels, delimiter)

    """
    import os
    import numpy as np
    from mindboggle.utils.compute import means_per_label, stats_per_label, \
        sum_per_label
    from mindboggle.utils.io_vtk import read_scalars, read_vtk, \
        apply_affine_transform
    from mindboggle.utils.io_table import write_columns

    # Make sure inputs are lists:
    if isinstance(labels_or_file, np.ndarray):
        labels = labels_or_file.tolist()
    elif isinstance(labels_or_file, list):
        labels = labels_or_file
    elif isinstance(labels_or_file, str):
        labels, name = read_scalars(labels_or_file)
    if isinstance(sulci, np.ndarray):
        sulci = sulci.tolist()
    if isinstance(fundi, np.ndarray):
        fundi = fundi.tolist()

    if not labels and not sulci and not fundi:
        import sys
        sys.exit('No feature data to tabulate in write_shape_stats().')

    #-------------------------------------------------------------------------
    # Feature lists, shape names, and shape files:
    #-------------------------------------------------------------------------
    # Feature lists:
    feature_lists = [labels, sulci, fundi]
    feature_names = ['label', 'sulcus', 'fundus']
    spectra_lists = [labels_spectra, sulci_spectra]
    spectra_ID_lists = [labels_spectra_IDs, sulci_spectra_IDs]
    spectra_names = ['label spectrum', 'sulcus spectrum']
    zernike_lists = [labels_zernike, sulci_zernike]
    zernike_ID_lists = [labels_zernike_IDs, sulci_zernike_IDs]
    zernike_names = ['label zernike', 'sulcus zernike']
    table_names = ['label_shapes.csv', 'sulcus_shapes.csv',
                   'fundus_shapes.csv']

    # Shape names corresponding to shape files below:
    shape_names = ['area', 'mean curvature', 'travel depth', 'geodesic depth',
                   'FreeSurfer convexity', 'FreeSurfer thickness']

    # Load shape files as a list of numpy arrays of per-vertex shape values:
    shape_files = [area_file, mean_curvature_file, travel_depth_file,
                   geodesic_depth_file, freesurfer_convexity_file,
                   freesurfer_thickness_file]
    shape_arrays = []
    column_names = []
    first_pass = True
    area_array = []
    for ishape, shape_file in enumerate(shape_files):
        if os.path.exists(shape_file):
            if first_pass:
                faces, lines, indices, points, npoints, scalars_array, name, \
                    input_vtk = read_vtk(shape_file, True, True)
                points = np.array(points)
                first_pass = False
                if affine_transform_file:
                    affine_points, \
                        foo1 = apply_affine_transform(affine_transform_file,
                                    points, transform_format, save_file=False)
                    affine_points = np.array(affine_points)
            else:
                scalars_array, name = read_scalars(shape_file, True, True)
            if scalars_array.size:
                shape_arrays.append(scalars_array)

                # Store area array:
                if ishape == 0:
                    area_array = scalars_array.copy()

    # Initialize table file names:
    label_table = ''
    sulcus_table = ''
    fundus_table = ''

    # Loop through features / tables:
    for itable, feature_list in enumerate(feature_lists):
        table_column_names = []

        #---------------------------------------------------------------------
        # For each feature, construct a table of average shape values:
        #---------------------------------------------------------------------
        if feature_list:
            feature_name = feature_names[itable]
            columns = []

            #-----------------------------------------------------------------
            # Mean positions in the original space:
            #-----------------------------------------------------------------
            # Compute mean position per feature:
            positions, sdevs, label_list, foo = means_per_label(points,
                feature_list, exclude_labels, area_array)

            # Append mean position per feature to columns:
            table_column_names.append('mean position')
            columns.append(positions)

            #-----------------------------------------------------------------
            # Mean positions in standard space:
            #-----------------------------------------------------------------
            if affine_transform_file:
                # Compute standard space mean position per feature:
                standard_positions, sdevs, label_list, \
                foo = means_per_label(affine_points,
                    feature_list, exclude_labels, area_array)

                # Append standard space mean position per feature to columns:
                table_column_names.append('mean position in standard space')
                columns.append(standard_positions)

            #-----------------------------------------------------------------
            # Loop through shape measures:
            #-----------------------------------------------------------------
            table_column_names.extend(column_names[:])
            for ishape, shape_array in enumerate(shape_arrays):
                shape_name = shape_names[ishape]
                print('  Compute statistics on {0} {1}'.
                      format(feature_name, shape_name))

                # Append shape names and values per feature to columns:
                pr = feature_name + ": " + shape_name + ": "
                if np.size(area_array):
                    po = " (weighted)"
                else:
                    po = ""
                #-------------------------------------------------------------
                # Append total feature areas to columns:
                #-------------------------------------------------------------
                if ishape == 0 and np.size(area_array):
                    sums, label_list = sum_per_label(shape_array,
                        feature_list, exclude_labels)
                    table_column_names.append(pr + 'total')
                    columns.append(sums)
                #-------------------------------------------------------------
                # Append feature shape statistics to columns:
                #-------------------------------------------------------------
                else:
                    medians, mads, means, sdevs, skews, kurts, \
                    lower_quarts, upper_quarts, \
                    label_list = stats_per_label(shape_array,
                        feature_list, exclude_labels, area_array, precision=1)

                    table_column_names.append(pr + 'median' + po)
                    table_column_names.append(pr + 'median absolute deviation' + po)
                    table_column_names.append(pr + 'mean' + po)
                    table_column_names.append(pr + 'standard deviation' + po)
                    table_column_names.append(pr + 'skew' + po)
                    table_column_names.append(pr + 'kurtosis' + po)
                    table_column_names.append(pr + 'lower quartile' + po)
                    table_column_names.append(pr + 'upper quartile' + po)
                    columns.append(medians)
                    columns.append(mads)
                    columns.append(means)
                    columns.append(sdevs)
                    columns.append(skews)
                    columns.append(kurts)
                    columns.append(lower_quarts)
                    columns.append(upper_quarts)

            #-----------------------------------------------------------------
            # Laplace-Beltrami spectra:
            #-----------------------------------------------------------------
            if itable in [0, 1]:
                spectra = spectra_lists[itable]
                spectra_name = spectra_names[itable]
                spectra_IDs = spectra_ID_lists[itable]

                # Order spectra into a list:
                spectrum_list = []
                for label in label_list:
                    if label in spectra_IDs:
                        spectrum = spectra[spectra_IDs.index(label)]
                        spectrum_list.append(spectrum)
                    else:
                        spectrum_list.append('')

                # Append spectral shape name and values to relevant columns:
                columns.append(spectrum_list)
                table_column_names.append(spectra_name)

            #-----------------------------------------------------------------
            # Zernike moments:
            #-----------------------------------------------------------------
            if itable in [0, 1]:
                zernike = zernike_lists[itable]
                zernike_name = zernike_names[itable]
                zernike_IDs = zernike_ID_lists[itable]

                # Order zernike into a list:
                spectrum_list = []
                for label in label_list:
                    if label in zernike_IDs:
                        spectrum = zernike[zernike_IDs.index(label)]
                        spectrum_list.append(spectrum)
                    else:
                        spectrum_list.append('')

                # Append Zernike shape name and values to relevant columns:
                columns.append(spectrum_list)
                table_column_names.append(zernike_name)

            #-----------------------------------------------------------------
            # Write labels/IDs and values to table:
            #-----------------------------------------------------------------
            # Write labels/IDs to table:
            output_table = os.path.join(os.getcwd(), table_names[itable])
            output_table = write_columns(label_list, feature_name, delimiter,
                                         quote=True, input_table='',
                                         output_table=output_table)

            # Append columns of shape values to table:
            if columns:
                write_columns(columns, table_column_names, delimiter,
                              quote=True, input_table=output_table,
                              output_table=output_table)

            if not os.path.exists(output_table):
                raise(IOError(output_table + " not found"))

            #-----------------------------------------------------------------
            # Return correct table file name:
            #-----------------------------------------------------------------
            if itable == 0:
                label_table = output_table
            elif itable == 1:
                sulcus_table = output_table
            elif itable == 2:
                fundus_table = output_table

    return label_table, sulcus_table, fundus_table


def write_vertex_measures(output_table, labels_or_file, sulci=[], fundi=[],
        affine_transform_file='', transform_format='itk',
        area_file='', mean_curvature_file='', travel_depth_file='',
        geodesic_depth_file='', freesurfer_convexity_file='',
        freesurfer_thickness_file='', delimiter=','):
    """
    Make a table of shape values per vertex.

    Note ::
        This function is tailored for Mindboggle outputs.

    Parameters
    ----------
    output_table : string
        output file (full path)
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
    mean_curvature_file :  string
        name of VTK file with mean curvature scalar values
    travel_depth_file :  string
        name of VTK file with travel depth scalar values
    geodesic_depth_file :  string
        name of VTK file with geodesic depth scalar values
    freesurfer_convexity_file :  string
        name of VTK file with FreeSurfer convexity scalar values
    freesurfer_thickness_file :  string
        name of VTK file with FreeSurfer thickness scalar values
    delimiter : string
        delimiter between columns, such as ','

    Returns
    -------
    output_table : table file name for vertex shape values

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.utils.io_table import write_vertex_measures
    >>> #
    >>> output_table = ''#vertex_shapes.csv'
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> labels_or_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> fundi_file = os.path.join(path, 'arno', 'features', 'fundi.vtk')
    >>> sulci, name = read_scalars(sulci_file)
    >>> fundi, name = read_scalars(fundi_file)
    >>> affine_transform_file = os.path.join(path, 'arno', 'mri',
    >>>     't1weighted_brain.MNI152Affine.txt')
    >>> transform_format = 'itk'
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> mean_curvature_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> travel_depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> geodesic_depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.geodesic_depth.vtk')
    >>> freesurfer_convexity_file = ''
    >>> freesurfer_thickness_file = ''
    >>> delimiter = ','
    >>> #
    >>> write_vertex_measures(output_table, labels_or_file, sulci, fundi,
    >>>     affine_transform_file, transform_format, area_file,
    >>>     mean_curvature_file, travel_depth_file, geodesic_depth_file,
    >>>     freesurfer_convexity_file, freesurfer_thickness_file, delimiter)

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_scalars, read_vtk, \
        apply_affine_transform
    from mindboggle.utils.io_table import write_columns

    # Make sure inputs are lists:
    if isinstance(labels_or_file, np.ndarray):
        labels = labels_or_file.tolist()
    elif isinstance(labels_or_file, list):
        labels = labels_or_file
    elif isinstance(labels_or_file, str):
        labels, name = read_scalars(labels_or_file)
    if isinstance(sulci, np.ndarray):
        sulci = sulci.tolist()
    if isinstance(fundi, np.ndarray):
        fundi = fundi.tolist()

    if not labels and not sulci and not fundi:
        import sys
        sys.exit('No feature data to tabulate in write_vertex_measures().')

    # Feature names and corresponding feature lists:
    feature_names = ['label', 'sulcus', 'fundus']
    feature_lists = [labels, sulci, fundi]

    # Shape names corresponding to shape files below:
    shape_names = ['area', 'mean curvature', 'travel depth', 'geodesic depth',
                   'FreeSurfer convexity', 'FreeSurfer thickness']

    # Load shape files as a list of numpy arrays of per-vertex shape values:
    shape_files = [area_file, mean_curvature_file, travel_depth_file,
                   geodesic_depth_file, freesurfer_convexity_file,
                   freesurfer_thickness_file]

    # Append columns of per-vertex scalar values:
    columns = []
    column_names = []
    for ifeature, values in enumerate(feature_lists):
        if values:
            columns.append(values)
            column_names.append(feature_names[ifeature])

    first_pass = True
    for ishape, shape_file in enumerate(shape_files):
        if os.path.exists(shape_file):
            if first_pass:
                u1, u2, u3, points, u4, scalars, u5, u6 = read_vtk(shape_file)
                columns.append(points)
                column_names.append('coordinates')
                first_pass = False
                if affine_transform_file:
                    affine_points, \
                        foo1 = apply_affine_transform(affine_transform_file,
                                    points, transform_format)
                    columns.append(affine_points)
                    column_names.append('coordinates in standard space')
            else:
                scalars, name = read_scalars(shape_file)
            if len(scalars):
                columns.append(scalars)
                column_names.append(shape_names[ishape])

    # Prepend with column of indices and write table
    if not output_table:
        output_table = os.path.join(os.getcwd(), 'vertices.csv')
    write_columns(range(len(columns[0])), 'index', delimiter, quote=True,
                  input_table='', output_table=output_table)
    write_columns(columns, column_names, delimiter, quote=True,
                  input_table=output_table, output_table=output_table)

    if not os.path.exists(output_table):
        raise(IOError(output_table + " not found"))

    return output_table


def write_face_vertex_averages(input_file, output_table='',
                               area_file='', delimiter=','):
    """
    Make table of average vertex values per face.

    Parameters
    ----------
    input_file : string
        name of VTK file with scalars to average
    area_file :  string
        name of VTK file with surface area scalar values
    delimiter : string
        delimiter between columns, such as ','
    output_table :  string
        output table filename

    Returns
    -------
    output_table :  string
        output table filename

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_table import write_face_vertex_averages
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> #input_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> #input_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.travel_depth.vtk')
    >>> input_file = os.path.join(path, 'arno', 'shapes', 'lh.thickness.vtk')
    >>> output_table = ''
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> delimiter = ','
    >>> #
    >>> write_face_vertex_averages(input_file, output_table, area_file, delimiter)

    """
    import os
    import numpy as np

    from mindboggle.utils.io_vtk import read_vtk, read_scalars
    from mindboggle.utils.io_table import write_columns

    faces, lines, indices, points, npoints, scalars, name, \
        input_vtk = read_vtk(input_file, True, True)
    if area_file:
        area_scalars, name = read_scalars(area_file, True, True)

    #---------------------------------------------------------------------
    # For each face, average vertex values:
    #---------------------------------------------------------------------
    columns = []
    for face in faces:
        values = []
        for index in face:
            if area_file:
                values.append(scalars[index] / area_scalars[index])
            else:
                values.append(scalars[index])
        columns.append(np.mean(values))

    #-----------------------------------------------------------------
    # Write to table:
    #-----------------------------------------------------------------
    if not output_table:
        output_table = os.path.join(os.getcwd(), 'average_face_values.csv')

    write_columns(columns, '', delimiter, quote=False, input_table='',
                  output_table=output_table)

    if not os.path.exists(output_table):
        raise(IOError(output_table + " not found"))

    return output_table


def write_average_face_values_per_label(input_indices_vtk,
                    input_values_vtk='', area_file='',
                    output_stem='', exclude_values=[-1], background_value=-1):
    """
    Write out a separate VTK file for each integer
    in (the first) scalar list of an input VTK file.
    Optionally write the values drawn from a second VTK file.

    Parameters
    ----------
    input_indices_vtk : string
        path of the input VTK file that contains indices as scalars
    input_values_vtk : string
        path of the input VTK file that contains values as scalars
    output_stem : string
        path and stem of the output VTK file
    exclude_values : list or array
        values to exclude
    background_value : integer or float
        background value in output VTK files
    scalar_name : string
        name of a lookup table of scalars values

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_table import write_average_face_values_per_label
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> input_indices_vtk = os.path.join(path, 'allen', 'labels', 'lh.DKTatlas100.gcs.vtk')
    >>> input_values_vtk = os.path.join(path, 'allen', 'shapes', 'lh.thickness.vtk')
    >>> area_file = os.path.join(path, 'allen', 'shapes', 'lh.pial.area.vtk')
    >>> output_stem = 'labels_thickness'
    >>> exclude_values = [-1]
    >>> background_value = -1
    >>> #
    >>> write_average_face_values_per_label(input_indices_vtk,
    >>>     input_values_vtk, area_file, output_stem, exclude_values, background_value)
    >>> #
    >>> # View:
    >>> #example_vtk = os.path.join(os.getcwd(), output_stem + '0.vtk')
    >>> #from mindboggle.utils.plots import plot_surfaces
    >>> #plot_surfaces(example_vtk)

    """
    import os
    import numpy as np
    from mindboggle.utils.io_vtk import read_scalars, read_vtk, write_vtk
    from mindboggle.utils.io_table import write_columns
    from mindboggle.utils.mesh import remove_faces

    # Load VTK file:
    faces, lines, indices, points, npoints, scalars, scalar_names, \
        foo1 = read_vtk(input_indices_vtk, True, True)
    if area_file:
        area_scalars, name = read_scalars(area_file, True, True)
    print("Explode the scalar list in {0}".
          format(os.path.basename(input_indices_vtk)))
    if input_values_vtk != input_indices_vtk:
        values, name = read_scalars(input_values_vtk, True, True)
        print("Explode the scalar list of values in {0} "
              "with the scalar list of indices in {1}".
              format(os.path.basename(input_values_vtk),
                     os.path.basename(input_indices_vtk)))
    else:
        values = np.copy(scalars)

    # Loop through unique (non-excluded) scalar values:
    unique_scalars = [int(x) for x in np.unique(scalars)
                      if x not in exclude_values]
    for scalar in unique_scalars:

        keep_indices = [x for sublst in faces for x in sublst]
        new_faces = remove_faces(faces, keep_indices)

        # Create array and indices for scalar value:
        select_scalars = np.copy(scalars)
        select_scalars[scalars != scalar] = background_value
        scalar_indices = [i for i,x in enumerate(select_scalars) if x==scalar]
        print("  Scalar {0}: {1} vertices".format(scalar, len(scalar_indices)))

        #---------------------------------------------------------------------
        # For each face, average vertex values:
        #---------------------------------------------------------------------
        output_table = os.path.join(os.getcwd(),
                                    output_stem+str(scalar)+'.csv')
        columns = []
        for face in new_faces:
            values = []
            for index in face:
                if area_file:
                    values.append(scalars[index] / area_scalars[index])
                else:
                    values.append(scalars[index])
            columns.append(np.mean(values))

        #-----------------------------------------------------------------
        # Write to table:
        #-----------------------------------------------------------------
        write_columns(columns, '', delimiter=',', quote=False,
                      input_table='', output_table=output_table)

        # Write VTK file with scalar value:
        #output_vtk = os.path.join(os.getcwd(), output_stem + str(scalar) + '.vtk')
        #write_vtk(output_vtk, points, indices, lines, new_faces,
        #          [select_values.tolist()], [output_scalar_name])

        if not os.path.exists(output_table):
            raise(IOError(output_table + " not found"))


def alternate_columns_from_tables(table_files, write_table=True,
                                  output_table='', delimiter=','):
    """
    Alternate columns from a list of tables and make a new table.

    This program assumes that the tables consist of only numbers.

    Parameters
    ----------
    table_files : list of strings
        table files (full paths)
    write_table : Boolean
        write output table?
    output_table : string
        output table file name
    delimiter : string
        delimiter between output table columns, such as ','

    Returns
    -------
    columns : list of lists of floats or integers
        columns of data
    output_table :  string
        output table file name

    Examples
    --------
    >>> from mindboggle.utils.io_table import alternate_columns_from_tables
    >>> table_files = ['/drop/MB/data/arno/tables/label_shapes.csv',
    >>>                '/drop/MB/data/arno/tables/label_shapes.csv']
    >>> write_table = True
    >>> output_table = ''
    >>> delimiter = ','
    >>> columns, output_table = alternate_columns_from_tables(table_files, write_table, output_table, delimiter)

    """
    import os
    import csv

    from mindboggle.utils.io_table import write_columns

    #-------------------------------------------------------------------------
    # Construct a list of all tables:
    #-------------------------------------------------------------------------
    tables = []
    for table_file in table_files:
        if not os.path.exists(table_file):
            raise(IOError(table_file + " not found"))
        else:
            reader = csv.reader(open(table_file, 'rb'),
                                delimiter=',', quotechar='"')
            tables.append([list(x) for x in zip(*reader)])

    #-------------------------------------------------------------------------
    # Alternate columns:
    #-------------------------------------------------------------------------
    columns = []
    for icolumn, column in enumerate(tables[0]):
        for table in tables:
            columns.append(table[icolumn])

    #-------------------------------------------------------------------------
    # Write table:
    #-------------------------------------------------------------------------
    if write_table and columns:
        if all([len(x) == len(columns[0]) for x in columns]):
            if not output_table:
                output_table = os.path.join(os.getcwd(),
                                            'alternating_columns.csv')
            write_columns(columns, column_names='', delimiter=delimiter,
                          quote=True, input_table='',
                          output_table=output_table)
        else:
            print('Columns do not have the same length.')

    return columns, output_table


def select_column_from_tables(tables, column_name, label_name='',
                              write_table=True, output_table='',
                              delimiter=','):
    """
    Select column from list of tables and make a new table.

    Parameters
    ----------
    tables : list of strings
        table files (full paths)
    column_name :  string
        column name to select
    label_name :  string
        column name for column with labels (if empty, no label column added)
    write_table : Boolean
        write output table?
    output_table : string
        output table file name
    delimiter : string
        delimiter between output table columns, such as ','

    Returns
    -------
    tables : list of strings
        input table files (full paths)
    columns : list of lists of floats or integers
        columns of data
    column_name :  string
        column name to select
    row_names : list of strings
        row labels (common strings in the label column of tables)
    row_names_title : string
        row_names column header
    output_table :  string
        output table file name

    Examples
    --------
    >>> from mindboggle.utils.io_table import select_column_from_tables
    >>> tables = ['/home/arno/mindboggled/tables/left/UM0029UMMR1R1_antsCorticalThickness/label_shapes.csv']
    >>> column_name = 'label: FreeSurfer thickness: mean (weighted)'
    >>> #tables = ['/drop/tables/left/UM0029UMMR2R1_FS11212/vertices.csv',
    >>> #          '/drop/tables/left/UM0029UMMR2R1_repositioned/vertices.csv']
    >>> #column_name = "FreeSurfer thickness"
    >>> label_name = 'label'
    >>> write_table = True
    >>> output_table = ''
    >>> delimiter = ','
    >>> #
    >>> select_column_from_tables(tables, column_name, label_name,
    >>>                           write_table, output_table, delimiter)

    """
    import os
    import sys
    import csv

    from mindboggle.utils.io_table import write_columns

    #-------------------------------------------------------------------------
    # Construct a table:
    #-------------------------------------------------------------------------
    columns = []
    row_names = []
    row_names_title = ''
    first = True
    for input_table in tables:

        #---------------------------------------------------------------------
        # Extract column from the table for each subject:
        #---------------------------------------------------------------------
        if not os.path.exists(input_table):
            raise(IOError(input_table + " not found"))
        else:
            reader = csv.reader(open(input_table, 'rb'),
                                delimiter=',', quotechar='"')
            input_columns = [list(x) for x in zip(*reader)]

            #-----------------------------------------------------------------
            # Extract column with table_name as its header:
            #-----------------------------------------------------------------
            icolumn_name = -1
            icolumn_label = -1
            for icolumn, column in enumerate(input_columns):
                # new_column = []
                # rms = ['"', "'"]
                # for c in column:
                #     for s in rms:
                #         c = c.strip()
                #         c = c.strip(s)
                #         c = c.strip()
                #     new_column.append(c)
                # column = new_column
                hdr = column[0]
                if hdr == column_name:
                    icolumn_name = icolumn
                elif hdr == label_name:
                    icolumn_label = icolumn
            if icolumn_name >= 0:
                columns.append(input_columns[icolumn_name][1::])
            else:
                sys.exit('No column name "{0}".'.format(column_name))

            #-----------------------------------------------------------------
            # Don't use the labels if label columns are unequal across tables:
            #-----------------------------------------------------------------
            if icolumn_label >= 0:
                if first:
                    row_names = input_columns[icolumn_label][1::]
                    row_names_title = input_columns[icolumn_label][0].strip()
                    first = False
                elif input_columns[icolumn_label][1::] != row_names:
                    print('Label columns are not the same across tables.')
                    label_name = False

    #-------------------------------------------------------------------------
    # Write table:
    #-------------------------------------------------------------------------
    if write_table and columns:
        if all([len(x) == len(columns[0]) for x in columns]):

            labels = [os.path.basename(x) for x in tables]

            if not output_table:
                rms = ['"', '`', ',', '/', '<', '>', '?', ':', ';',
                       '(', ')', '[', ']', '{', '}', '|', '\\',
                       '!', '@', '#', '$', '%', '^', '&', '*', '=', '_']
                table_name = column_name
                for s in rms:
                    table_name = table_name.replace(s, '_')
                table_name = table_name.replace(' ', '')
                output_table = os.path.join(os.getcwd(), table_name + '.csv')
                output_table = output_table.replace('_.', '.')

            if label_name:
                write_columns(row_names, row_names_title, delimiter, quote=True,
                              input_table='', output_table=output_table)
                write_columns(columns, labels, delimiter, quote=True,
                              input_table=output_table, output_table=output_table)
            else:
                write_columns(columns, labels, delimiter, quote=True,
                              input_table='', output_table=output_table)

        else:
            print('Not saving table.')

    return tables, columns, column_name, row_names, row_names_title, output_table


def select_column_from_mindboggle_tables(subjects, hemi, tables_dir,
        table_name, column_name, label_name='label',
        write_table=True, output_table='', delimiter=','):
    """
    Select column from Mindboggle shape tables and make a new table.

    For example, extract the median travel depth column for a given feature
    (label region or sulcus ID) across a set of subjects, and make a table.

    Expects: <tables_dir>/['left','right']/<subject>/<table_name>

    Parameters
    ----------
    subjects :  list of strings
        names of subjects processed by Mindboggle
    hemi :  string
        hemisphere in {'left', 'right}
    tables_dir : string
        name of Mindboggle tables directory
    table_name : string
        name of Mindboggle table file
    column_name :  string
        column name to select
    label_name :  string
        column name for column with labels (if empty, no label column added)
    write_table : Boolean
        write output table?
    output_table : string
        output table file name
    delimiter : string
        delimiter between output table columns, such as ','

    Returns
    -------
    tables : list of strings
        input table files (full paths)
    columns : list of lists of floats or integers
        columns of data
    column_name :  string
        column name to select
    row_names : list of strings
        row labels (common strings in the label column of tables)
    row_names_title : string
        row_names column header
    output_table :  string
        output table file name

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_table import select_column_from_mindboggle_tables
    >>> subjects = ['Twins-2-1']
    >>> hemi = 'left'
    >>> tables_dir = os.path.join(os.environ['HOME'], 'mindboggled', 'tables')
    >>> table_name = "label_shapes.csv"
    >>> column_name = "label: FreeSurfer thickness: median (weighted)"
    >>> label_name = 'label'
    >>> write_table = True
    >>> output_table = ''
    >>> delimiter = ','
    >>> #
    >>> select_column_from_mindboggle_tables(subjects, hemi, tables_dir,
    >>>     table_name, column_name, label_name,
    >>>     write_table, output_table, delimiter)

    """
    import os

    from mindboggle.utils.io_table import select_column_from_tables

    #-------------------------------------------------------------------------
    # Construct list of Mindboggle shape table file names:
    #-------------------------------------------------------------------------
    tables = []
    for subject in subjects:
        tables.append(os.path.join(tables_dir, hemi, subject, table_name))

    #-------------------------------------------------------------------------
    # Extract columns and construct new table:
    #-------------------------------------------------------------------------
    tables, columns, column_name, row_names, row_names_title, \
    output_table = select_column_from_tables(tables, column_name, label_name,
                                             write_table, output_table,
                                             delimiter)

    return tables, columns, column_name, row_names, row_names_title, output_table


def xyz2nii(input_xyz_file, output_nii_file='', origin=[], pad=10):
    """
    Convert [x,y,z] coordinate file to nifti (nii.gz) volume file.

    Parameters
    ----------
    input_xyz_file : string
        input [x,y,z] coordinate text file
    output_nii_file : string
        output nifti (nii.gz) volume file

    Returns
    -------
    output_nii_file : string
        output nifti (nii.gz) volume file

    Examples
    --------
    >>> from mindboggle.utils.io_table import xyz2nii
    >>> input_xyz_file = '/Users/arno/Dropbox/MSSM/Nebojsa/face.xyz.txt'
    >>> origin = []
    >>> pad = 10
    >>> output_nii_file = ''
    >>> xyz2nii(input_xyz_file)

    """
    import os
    import numpy as np
    import nibabel as nb

    # Load coordinates and scalars:
    XYZscalars = np.loadtxt(input_xyz_file)
    XYZ = np.round(XYZscalars[:, 0:3])
    #scalars = XYZscalars[:, 3::]

    if origin:
        XYZ -= origin

    XYZ += np.abs(np.min(XYZ, axis=0)) + [pad, pad, pad]
    XYZ = np.round(XYZ)
    dims = np.max(XYZ, axis=0) + [pad, pad, pad]
    data = np.zeros(dims)

    # Loop through rows or array and write 1s in image volume:
    for irow, xyz in enumerate(XYZ):
        data[xyz[0], xyz[1], xyz[2]] = 1

    # Write output image volume:
    if not output_nii_file:
        output_nii_file = os.path.join(os.getcwd(), 'xyz.nii.gz')
    img = nb.Nifti1Image(data, affine=np.eye(4,4))
    img.to_filename(output_nii_file)

    return output_nii_file
