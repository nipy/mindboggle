#!/usr/bin/env python
"""
Functions for creating tables.


Authors:
    - Arno Klein, 2012-2016  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2016,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""


def write_shape_stats(labels_or_file=[], sulci=[], fundi=[],
        affine_transform_files=[], inverse_booleans=[], transform_format='itk',
        area_file='', normalize_by_area=False, mean_curvature_file='',
        travel_depth_file='', geodesic_depth_file='',
        freesurfer_thickness_file='', freesurfer_curvature_file='',
        freesurfer_sulc_file='',
        labels_spectra=[], labels_spectra_IDs=[],
        sulci_spectra=[], sulci_spectra_IDs=[],
        labels_zernike=[], labels_zernike_IDs=[],
        sulci_zernike=[], sulci_zernike_IDs=[],
        exclude_labels=[-1], verbose=False):
    """
    Make tables of shape statistics per label, sulcus, and/or fundus.

    Note ::
        This function is tailored for Mindboggle outputs.

    Parameters
    ----------
    labels_or_file : list or string
        label number for each vertex or name of VTK file with index scalars
    sulci : list of integers
        indices to sulci, one per vertex, with -1 indicating no sulcus
    fundi : list of integers
        indices to fundi, one per vertex, with -1 indicating no fundus
    affine_transform_files : list of strings
        affine transform files to standard space
    inverse_booleans : list of of zeros and ones
        for each transform, 1 to take the inverse, else 0
    transform_format : string
        format for transform file
        Ex: 'txt' for text, 'itk' for ITK, and 'mat' for Matlab format
    area_file :  string
        name of VTK file with surface area scalar values
    normalize_by_area : bool
        normalize all shape measures by area of label/feature? (UNTESTED)
    mean_curvature_file :  string
        name of VTK file with mean curvature scalar values
    travel_depth_file :  string
        name of VTK file with travel depth scalar values
    geodesic_depth_file :  string
        name of VTK file with geodesic depth scalar values
    freesurfer_thickness_file :  string
        name of VTK file with FreeSurfer thickness scalar values
    freesurfer_curvature_file :  string
        name of VTK file with FreeSurfer curvature (curv) scalar values
    freesurfer_sulc_file :  string
        name of VTK file with FreeSurfer convexity (sulc) scalar values
    labels_spectra : list of lists of floats
        Laplace-Beltrami spectra for each labeled region
    labels_spectra_IDs : list of integers
        unique labels for labels_spectra
    sulci_spectra : list of lists of floats
        Laplace-Beltrami spectra for each sulcus
    sulci_spectra_IDs : list of integers
        unique sulcus IDs for sulci_spectra
    labels_zernike : list of lists of floats
        Zernike moments for each labeled region
    labels_zernike_IDs : list of integers
        unique labels for labels_zernike
    sulci_zernike : list of lists of floats
        Zernike moments for each sulcus
    sulci_zernike_IDs : list of integers
        unique sulcus IDs for sulci_zernike
    exclude_labels : list of lists of integers
        indices to be excluded (in addition to -1)
    verbose : bool
        print statements?

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
    >>> from mindboggle.mio.tables import write_shape_stats
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> label_file = fetch_data(urls['left_freesurfer_labels'])
    >>> sulci_file = fetch_data(urls['left_sulci'])
    >>> fundi_file = fetch_data(urls['left_fundus_per_sulcus'])
    >>> mean_curvature_file = fetch_data(urls['left_mean_curvature'])
    >>> travel_depth_file = fetch_data(urls['left_travel_depth'])
    >>> geodesic_depth_file = fetch_data(urls['left_geodesic_depth'])
    >>> area_file = fetch_data(urls['left_area'])
    >>> freesurfer_thickness_file = ''
    >>> freesurfer_curvature_file = ''
    >>> freesurfer_sulc_file = ''
    >>> sulci, name = read_scalars(sulci_file)
    >>> fundi, name = read_scalars(fundi_file)
    >>> affine_transform_files = []
    >>> inverse_booleans = []
    >>> transform_format = 'itk'
    >>> normalize_by_area = False
    >>> labels, name = read_scalars(label_file)
    >>> labels_spectra = []
    >>> labels_spectra_IDs = []
    >>> sulci_spectra = []
    >>> sulci_spectra_IDs = []
    >>> labels_zernike = []
    >>> labels_zernike_IDs = []
    >>> sulci_zernike = []
    >>> sulci_zernike_IDs = []
    >>> exclude_labels = [-1]
    >>> verbose = False
    >>> label_table, sulcus_table, fundus_table = write_shape_stats(label_file,
    ...     sulci, fundi, affine_transform_files, inverse_booleans,
    ...     transform_format, area_file, normalize_by_area,
    ...     mean_curvature_file, travel_depth_file, geodesic_depth_file,
    ...     freesurfer_thickness_file, freesurfer_curvature_file,
    ...     freesurfer_sulc_file, labels_spectra, labels_spectra_IDs,
    ...     sulci_spectra, sulci_spectra_IDs, labels_zernike,
    ...     labels_zernike_IDs, sulci_zernike, sulci_zernike_IDs,
    ...     exclude_labels, verbose)

    """
    import os
    import numpy as np
    import pandas as pd

    from mindboggle.guts.compute import stats_per_label
    from mindboggle.guts.compute import means_per_label
    from mindboggle.guts.compute import sum_per_label
    from mindboggle.mio.vtks import read_scalars, read_vtk
    from mindboggle.mio.vtks import apply_affine_transforms
    from mindboggle.mio.labels import DKTprotocol

    dkt = DKTprotocol()

    # Make sure inputs are lists:
    if isinstance(labels_or_file, np.ndarray):
        labels = [int(x) for x in labels_or_file]
    elif isinstance(labels_or_file, list):
        labels = labels_or_file
    elif isinstance(labels_or_file, str):
        labels, name = read_scalars(labels_or_file)
    if isinstance(sulci, np.ndarray):
        sulci = [int(x) for x in sulci]
    if isinstance(fundi, np.ndarray):
        fundi = [int(x) for x in fundi]

    if not labels and not sulci and not fundi:
        raise IOError('No feature data to tabulate in write_shape_stats().')

    spectrum_start = 1  # Store all columns of spectral components (0),
                        # or start from higher frequency components (>=1)

    # ------------------------------------------------------------------------
    # Feature lists, shape names, and shape files:
    # ------------------------------------------------------------------------
    # Feature lists:
    feature_lists = [labels, sulci, fundi]
    feature_names = ['label', 'sulcus', 'fundus']
    spectra_lists = [labels_spectra, sulci_spectra]
    spectra_ID_lists = [labels_spectra_IDs, sulci_spectra_IDs]
    zernike_lists = [labels_zernike, sulci_zernike]
    zernike_ID_lists = [labels_zernike_IDs, sulci_zernike_IDs]
    table_names = ['label_shapes.csv', 'sulcus_shapes.csv',
                   'fundus_shapes.csv']

    # Shape names corresponding to shape files below:
    shape_names = ['area', 'travel depth', 'geodesic depth',
                   'mean curvature', 'freesurfer curvature',
                   'freesurfer thickness', 'freesurfer convexity (sulc)']

    # Load shape files as a list of numpy arrays of per-vertex shape values:
    shape_files = [area_file, travel_depth_file, geodesic_depth_file,
                   mean_curvature_file, freesurfer_curvature_file,
                   freesurfer_thickness_file, freesurfer_sulc_file]
    shape_arrays = []
    first_pass = True
    area_array = []

    for ishape, shape_file in enumerate(shape_files):
        if os.path.exists(shape_file):
            if first_pass:
                points, indices, lines, faces, scalars_array, scalar_names, \
                    npoints, input_vtk = read_vtk(shape_file, True, True)
                points = np.array(points)
                first_pass = False
                if affine_transform_files and transform_format:
                    affine_points, \
                        foo1 = apply_affine_transforms(affine_transform_files,
                                    inverse_booleans, transform_format,
                                    points, vtk_file_stem='')
            else:
                scalars_array, name = read_scalars(shape_file, True, True)
            if scalars_array.size:
                shape_arrays.append(scalars_array)

                # Store area array:
                if ishape == 0:
                    area_array = scalars_array.copy()

    if normalize_by_area:
        use_area = area_array
    else:
        use_area = []

    # Initialize table file names:
    label_table = ''
    sulcus_table = ''
    fundus_table = ''

    # Loop through features / tables:
    for itable, feature_list in enumerate(feature_lists):
        column_names = []

        # ----------------------------------------------------------------
        # Label names:
        # ----------------------------------------------------------------
        label_title = 'name'
        if itable == 0:
            label_numbers = dkt.cerebrum_cortex_DKT31_numbers
            label_names = dkt.cerebrum_cortex_DKT31_names
        elif itable in [1, 2]:
            label_numbers = dkt.sulcus_numbers
            label_names = dkt.sulcus_names
        else:
            label_numbers = []
            label_names = []
        include_labels = label_numbers
        nlabels = len(label_numbers)

        # --------------------------------------------------------------------
        # For each feature, construct a table of average shape values:
        # --------------------------------------------------------------------
        if feature_list:
            feature_name = feature_names[itable]
            columns = []

            # ----------------------------------------------------------------
            # Loop through shape measures:
            # ----------------------------------------------------------------
            column_names.extend(column_names[:])
            for ishape, shape_array in enumerate(shape_arrays):
                shape = shape_names[ishape]
                if verbose:
                    print('  Compute statistics on {0} {1}...'.
                        format(feature_name, shape))
                # ------------------------------------------------------------
                # Append feature areas to columns:
                # ------------------------------------------------------------
                if ishape == 0 and np.size(area_array):
                    sums, label_list = sum_per_label(shape_array,
                        feature_list, include_labels, exclude_labels)
                    column_names.append(shape)
                    columns.append(sums)
                # ------------------------------------------------------------
                # Append feature shape statistics to columns:
                # ------------------------------------------------------------
                else:
                    medians, mads, means, sdevs, skews, kurts, \
                    lower_quarts, upper_quarts, \
                    label_list = stats_per_label(shape_array, feature_list,
                                        include_labels, exclude_labels,
                                        area_array, precision=1)

                    column_names.append(shape + ': median')
                    column_names.append(shape + ': MAD')
                    column_names.append(shape + ': mean')
                    column_names.append(shape + ': SD')
                    column_names.append(shape + ': skew')
                    column_names.append(shape + ': kurtosis')
                    column_names.append(shape + ': 25%')
                    column_names.append(shape + ': 75%')
                    columns.append(medians)
                    columns.append(mads)
                    columns.append(means)
                    columns.append(sdevs)
                    columns.append(skews)
                    columns.append(kurts)
                    columns.append(lower_quarts)
                    columns.append(upper_quarts)

            # ----------------------------------------------------------------
            # Mean positions in the original space:
            # ----------------------------------------------------------------
            # Compute mean position per feature:
            positions, sdevs, label_list, foo = means_per_label(points,
                feature_list, include_labels, exclude_labels, use_area)

            # Append mean x,y,z position per feature to columns:
            xyz_positions = np.asarray(positions)
            for ixyz, xyz in enumerate(['x','y','z']):
                column_names.append('mean position: {0}'.format(xyz))
                columns.append(xyz_positions[:, ixyz].tolist())

            # ----------------------------------------------------------------
            # Mean positions in standard space:
            # ----------------------------------------------------------------
            if affine_transform_files and transform_format:
                # Compute standard space mean position per feature:
                standard_positions, sdevs, label_list, \
                foo = means_per_label(affine_points,
                    feature_list, include_labels, exclude_labels, use_area)

                # Append standard space x,y,z position per feature to columns:
                xyz_std_positions = np.asarray(standard_positions)
                for ixyz, xyz in enumerate(['x','y','z']):
                    column_names.append('mean position in standard space:'
                                        ' {0}'.format(xyz))
                    columns.append(xyz_std_positions[:, ixyz].tolist())

            # ----------------------------------------------------------------
            # Laplace-Beltrami spectra:
            # ----------------------------------------------------------------
            if itable in [0, 1]:
                spectra = spectra_lists[itable]
                if spectra:
                    spectra_IDs = spectra_ID_lists[itable]

                    # Construct a matrix of spectra:
                    len_spectrum = len(spectra[0])
                    spectrum_matrix = np.zeros((nlabels, len_spectrum))
                    for ilabel, label in enumerate(include_labels):
                        if label in spectra_IDs:
                            spectrum = spectra[spectra_IDs.index(label)]
                            spectrum_matrix[ilabel, 0:len_spectrum] = spectrum

                    # Append spectral shape name and values to columns:
                    for ispec in range(spectrum_start, len_spectrum):
                        columns.append(spectrum_matrix[:, ispec].tolist())
                        column_names.append('Laplace-Beltrami spectrum:'
                                            ' component {0}'.format(ispec+1))

            # ----------------------------------------------------------------
            # Zernike moments:
            # ----------------------------------------------------------------
            if itable in [0, 1]:
                zernike = zernike_lists[itable]
                if zernike:
                    zernike_IDs = zernike_ID_lists[itable]

                    # Construct a matrix of Zernike moments:
                    len_moments = len(zernike[0])
                    moments_matrix = np.zeros((nlabels, len_moments))
                    for ilabel, label in enumerate(include_labels):
                        if label in zernike_IDs:
                            moments = zernike[zernike_IDs.index(label)]
                            moments_matrix[ilabel, 0:len_moments] = moments

                    # Append Zernike shape name and values to columns:
                    for imoment in range(0, len_moments):
                        columns.append(moments_matrix[:, imoment].tolist())
                        column_names.append('Zernike moments: component {0}'.
                                            format(imoment+1))

            # ----------------------------------------------------------------
            # Write labels/IDs and values to table:
            # ----------------------------------------------------------------
            # Write labels/IDs to table:
            output_table = os.path.join(os.getcwd(), table_names[itable])

            if columns:
                df1 = pd.DataFrame({'ID': label_numbers})
                df2 = pd.DataFrame(np.transpose(columns),
                                   columns = column_names)
                df = pd.concat([df1, df2], axis=1)
                if label_names:
                    df0 = pd.DataFrame({'name': label_names})
                    df = pd.concat([df0, df], axis=1)
                df.to_csv(output_table, index=False)

            if not os.path.exists(output_table):
                raise IOError(output_table + " not found")

            # ----------------------------------------------------------------
            # Return correct table file name:
            # ----------------------------------------------------------------
            if itable == 0:
                label_table = output_table
            elif itable == 1:
                sulcus_table = output_table
            elif itable == 2:
                fundus_table = output_table

    return label_table, sulcus_table, fundus_table


def write_vertex_measures(output_table, labels_or_file, sulci=[], fundi=[],
        affine_transform_files=[], inverse_booleans=[],
        transform_format='itk', area_file='', mean_curvature_file='',
        travel_depth_file='', geodesic_depth_file='',
        freesurfer_thickness_file='', freesurfer_curvature_file='',
        freesurfer_sulc_file=''):
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
    affine_transform_files : list of strings
        affine transform files to standard space
    inverse_booleans : list of of zeros and ones
        for each transform, 1 to take the inverse, else 0
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
    freesurfer_thickness_file :  string
        name of VTK file with FreeSurfer thickness scalar values
    freesurfer_curvature_file :  string
        name of VTK file with FreeSurfer curvature (curv) scalar values
    freesurfer_sulc_file :  string
        name of VTK file with FreeSurfer convexity (sulc) scalar values

    Returns
    -------
    output_table : table file name for vertex shape values

    Examples
    --------
    >>> import os
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.tables import write_vertex_measures
    >>> output_table = '' #vertex_shapes.csv'
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> labels_or_file = fetch_data(urls['left_freesurfer_labels'])
    >>> sulci_file = fetch_data(urls['left_sulci'])
    >>> fundi_file = sulci_file #''#fetch_data(urls['left_fundus_per_sulcus'])
    >>> mean_curvature_file = fetch_data(urls['left_mean_curvature'])
    >>> travel_depth_file = fetch_data(urls['left_travel_depth'])
    >>> geodesic_depth_file = fetch_data(urls['left_geodesic_depth'])
    >>> area_file = fetch_data(urls['left_area'])
    >>> freesurfer_thickness_file = fetch_data(urls['left_freesurfer_thickness'])
    >>> freesurfer_curvature_file = fetch_data(urls['left_freesurfer_curvature'])
    >>> freesurfer_sulc_file = fetch_data(urls['left_freesurfer_sulc'])
    >>> sulci, name = read_scalars(sulci_file)
    >>> if fundi_file:
    ...     fundi, name = read_scalars(fundi_file)
    ... else:
    ...     fundi = []
    >>> affine_transform_file = fetch_data(urls['affine_mni_transform'])
    >>> inverse_booleans = [1]
    >>> transform_format = 'itk'
    >>> swap_xy = True
    >>> affine_rename = affine_transform_file + '.txt'
    >>> os.rename(affine_transform_file, affine_rename)
    >>> os.rename(labels_or_file, labels_or_file + '.vtk')
    >>> os.rename(area_file, area_file + '.vtk')
    >>> os.rename(mean_curvature_file, mean_curvature_file + '.vtk')
    >>> os.rename(travel_depth_file, travel_depth_file + '.vtk')
    >>> os.rename(geodesic_depth_file, geodesic_depth_file + '.vtk')
    >>> os.rename(freesurfer_thickness_file, freesurfer_thickness_file + '.vtk')
    >>> os.rename(freesurfer_curvature_file, freesurfer_curvature_file + '.vtk')
    >>> os.rename(freesurfer_sulc_file, freesurfer_sulc_file + '.vtk')
    >>> labels_or_file = labels_or_file + '.vtk'
    >>> area_file = area_file + '.vtk'
    >>> mean_curvature_file = mean_curvature_file + '.vtk'
    >>> travel_depth_file = travel_depth_file + '.vtk'
    >>> geodesic_depth_file = geodesic_depth_file + '.vtk'
    >>> freesurfer_thickness_file = freesurfer_thickness_file + '.vtk'
    >>> freesurfer_curvature_file = freesurfer_curvature_file + '.vtk'
    >>> freesurfer_sulc_file = freesurfer_sulc_file + '.vtk'
    >>> affine_transform_files = [] # [affine_rename] # requires ANTs to test
    >>> output_table = write_vertex_measures(output_table, labels_or_file,
    ...     sulci, fundi, affine_transform_files, inverse_booleans,
    ...     transform_format, area_file, mean_curvature_file,
    ...     travel_depth_file, geodesic_depth_file, freesurfer_thickness_file,
    ...     freesurfer_curvature_file, freesurfer_sulc_file)

    """
    import os
    import numpy as np
    import pandas as pd

    from mindboggle.mio.vtks import read_scalars, read_vtk, \
        apply_affine_transforms

    # Make sure inputs are lists:
    if isinstance(labels_or_file, np.ndarray):
        labels = [int(x) for x in labels_or_file]
    elif isinstance(labels_or_file, list):
        labels = labels_or_file
    elif isinstance(labels_or_file, str):
        labels, name = read_scalars(labels_or_file)
    if isinstance(sulci, np.ndarray):
        sulci = [int(x) for x in sulci]
    if isinstance(fundi, np.ndarray):
        fundi = [int(x) for x in fundi]

    if not labels and not sulci and not fundi:
        raise IOError('No feature data to tabulate in write_vertex_measures().')

    # Feature names and corresponding feature lists:
    feature_names = ['label ID', 'sulcus ID', 'fundus ID']
    feature_lists = [labels, sulci, fundi]

    # Shape names corresponding to shape files below:
    shape_names = ['area', 'travel depth', 'geodesic depth',
                   'mean curvature', 'freesurfer curvature',
                   'freesurfer thickness', 'freesurfer convexity (sulc)']

    # Load shape files as a list of numpy arrays of per-vertex shape values:
    shape_files = [area_file, travel_depth_file, geodesic_depth_file,
                   mean_curvature_file, freesurfer_curvature_file,
                   freesurfer_thickness_file, freesurfer_sulc_file]

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

                # Append x,y,z position per vertex to columns:
                points, indices, lines, faces, scalars, scalar_names, \
                    npoints, input_vtk = read_vtk(shape_file)
                xyz_positions = np.asarray(points)
                for ixyz, xyz in enumerate(['x','y','z']):
                    column_names.append('position: {0}'.format(xyz))
                    columns.append(xyz_positions[:, ixyz].tolist())
                first_pass = False

                # Append standard space x,y,z position to columns:
                if affine_transform_files and transform_format:
                    affine_points, \
                        foo1 = apply_affine_transforms(affine_transform_files,
                                    inverse_booleans, transform_format,
                                    points, vtk_file_stem='')
                    xyz_std_positions = affine_points
                    for ixyz, xyz in enumerate(['x','y','z']):
                        column_names.append('position in standard space:'
                                            ' {0}'.format(xyz))
                        columns.append(xyz_std_positions[:, ixyz].tolist())
            else:
                scalars, name = read_scalars(shape_file)
            if len(scalars):
                columns.append(scalars)
                column_names.append(shape_names[ishape])

    # Prepend with column of indices and write table
    if not output_table:
        output_table = os.path.join(os.getcwd(), 'vertices.csv')

    df = pd.DataFrame(np.transpose(columns), columns = column_names)
    df.to_csv(output_table, index=False)

    if not os.path.exists(output_table):
        raise IOError(output_table + " not found")

    return output_table


def write_face_vertex_averages(input_file, output_table='', area_file=''):
    """
    Make table of average vertex values per face
    (divided by face area if area_file provided).

    Parameters
    ----------
    input_file : string
        name of VTK file with scalars to average
    area_file :  string
        name of VTK file with surface area scalar values
    output_table :  string
        output table filename

    Returns
    -------
    output_table :  string
        output table filename

    Examples
    --------
    >>> from mindboggle.mio.tables import write_face_vertex_averages
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> input_file = fetch_data(urls['left_travel_depth'])
    >>> area_file = fetch_data(urls['left_area'])
    >>> output_table = ''
    >>> output_table = write_face_vertex_averages(input_file, output_table,
    ...                                           area_file)

    """
    import os
    import numpy as np
    import pandas as pd

    from mindboggle.mio.vtks import read_vtk, read_scalars

    points, indices, lines, faces, scalars, scalar_names, \
        npoints, input_vtk = read_vtk(input_file, True, True)
    if area_file:
        area_scalars, name = read_scalars(area_file, True, True)

    # --------------------------------------------------------------------
    # For each face, average vertex values:
    # --------------------------------------------------------------------
    columns = []
    for face in faces:
        values = []
        for index in face:
            if area_file:
                values.append(scalars[index] / area_scalars[index])
            else:
                values.append(scalars[index])
        columns.append(np.mean(values))

    # ----------------------------------------------------------------
    # Write to table:
    # ----------------------------------------------------------------
    if not output_table:
        output_table = os.path.join(os.getcwd(), 'average_face_values.csv')

    df = pd.DataFrame({'': columns})
    df.to_csv(output_table, index=False)

    if not os.path.exists(output_table):
        raise IOError(output_table + " not found")

    return output_table


def write_average_face_values_per_label(input_indices_vtk,
        input_values_vtk='', area_file='', output_stem='',
        exclude_values=[-1], background_value=-1, verbose=False):
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
    verbose : bool
        print statements?

    Examples
    --------
    >>> import os
    >>> from mindboggle.mio.tables import write_average_face_values_per_label
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> input_indices_vtk = fetch_data(urls['left_freesurfer_labels'])
    >>> input_values_vtk = fetch_data(urls['left_mean_curvature'])
    >>> area_file = fetch_data(urls['left_area'])
    >>> output_stem = 'labels_thickness'
    >>> exclude_values = [-1]
    >>> background_value = -1
    >>> verbose = False
    >>> write_average_face_values_per_label(input_indices_vtk,
    ...     input_values_vtk, area_file, output_stem, exclude_values,
    ...     background_value, verbose)

    View vtk file (skip test):

    >>> from mindboggle.mio.plots import plot_surfaces
    >>> example_vtk = os.path.join(os.getcwd(), output_stem + '0.vtk')
    >>> plot_surfaces(example_vtk) # doctest: +SKIP

    """
    import os
    import numpy as np
    import pandas as pd

    from mindboggle.mio.vtks import read_scalars, read_vtk, write_vtk
    from mindboggle.guts.mesh import keep_faces

    # Load VTK file:
    points, indices, lines, faces, scalars, scalar_names, npoints, \
        input_vtk = read_vtk(input_indices_vtk, True, True)
    if area_file:
        area_scalars, name = read_scalars(area_file, True, True)
    if verbose:
        print("Explode the scalar list in {0}".
            format(os.path.basename(input_indices_vtk)))
    if input_values_vtk != input_indices_vtk:
        if verbose:
            print("Explode the scalar list of values in {0} "
                  "with the scalar list of indices in {1}".
                format(os.path.basename(input_values_vtk),
                       os.path.basename(input_indices_vtk)))

    # Loop through unique (non-excluded) scalar values:
    unique_scalars = [int(x) for x in np.unique(scalars)
                      if x not in exclude_values]
    for scalar in unique_scalars:

        keep_indices = [x for sublst in faces for x in sublst]
        new_faces = keep_faces(faces, keep_indices)

        # Create array and indices for scalar value:
        select_scalars = np.copy(scalars)
        select_scalars[scalars != scalar] = background_value
        scalar_indices = [i for i,x in enumerate(select_scalars) if x==scalar]
        if verbose:
            print("  Scalar {0}: {1} vertices".format(scalar,
                                                      len(scalar_indices)))

        # --------------------------------------------------------------------
        # For each face, average vertex values:
        # --------------------------------------------------------------------
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

        # ----------------------------------------------------------------
        # Write to table:
        # ----------------------------------------------------------------
        df = pd.DataFrame({'': columns})
        df.to_csv(output_table, index=False)
        if not os.path.exists(output_table):
            raise IOError(output_table + " not found")


def select_column_from_tables(tables, index=0, write_table=True,
                              output_table=''):
    """
    Select column from list of tables, make a new table.

    Note: If more than one table, column must be of the same length.

    Parameters
    ----------
    tables : list of strings
        table files (full paths)
    index : integer
        index for column to select (from each table)
    write_table : bool
        write output table?
    output_table : string
        output table file name

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
    row_stats :  list
        row statistics
    row_stats_names :  list
        names of column statistics
    output_table :  string
        output table file name

    Examples
    --------
    >>> from mindboggle.mio.tables import select_column_from_tables
    >>> from mindboggle.mio.fetch_data import prep_tests
    >>> urls, fetch_data = prep_tests()
    >>> tables = [fetch_data(urls['thickinthehead_freesurfer_labels_table']),
    ...           fetch_data(urls['thickinthehead_freesurfer_labels_table'])]
    >>> index = 2
    >>> write_table = True
    >>> output_table = ''
    >>> output = select_column_from_tables(tables, index, write_table,
    ...                                    output_table)
    >>> columns = output[1][0]
    >>> columns[0]
    2.8010000000000002
    >>> columns[1]
    3.9430000000000001
    >>> columns[2]
    4.0280000000000005

    """
    import os
    import pandas as pd

    # ------------------------------------------------------------------------
    # Construct a table:
    # ------------------------------------------------------------------------
    columns = []
    for input_table in tables:

        # --------------------------------------------------------------------
        # Extract column from the table for each subject:
        # --------------------------------------------------------------------
        if not os.path.exists(input_table):
            raise IOError(input_table + " not found")
        else:
            input_columns = pd.read_csv(input_table)
            columns.append(input_columns.iloc[:, index])

    # ------------------------------------------------------------------------
    # Write tables:
    # ------------------------------------------------------------------------
    if write_table and columns:
        if all([len(x) == len(columns[0]) for x in columns]):
            if not output_table:
                output_table = os.path.join(os.getcwd(),
                                            'select_column_from_tables.csv')
            df = pd.DataFrame({'': columns})
            df.to_csv(output_table, index=False)
        else:
            raise IOError('Not saving table.')

    return tables, columns, output_table


def select_column_from_mindboggle_tables(subjects, hemi, index, tables_dir,
        table_name, is_surface_table=True, write_table=True, output_table=''):
    """
    Select column from Mindboggle shape tables and make a new table.

    For example, extract the median travel depth column for the label regions
    across a set of subjects, and make a new table.

    Expects::
        <tables_dir>/<subject>/tables/['left','right']_cortical_surface/<table_name>

    Parameters
    ----------
    subjects :  list of strings
        names of subjects processed by Mindboggle
    hemi :  string
        hemisphere in {'left', 'right}
    index : integer
        index for column to select
    tables_dir : string
        name of Mindboggle tables directory
    table_name : string
        name of Mindboggle table file
    is_surface_table : bool
        if True, use path to surface tables
    write_table : bool
        write output table?
    output_table : string
        output table file name

    Returns
    -------
    tables : list of strings
        input table files (full paths)
    columns : list of lists of floats or integers
        columns of data
    output_table :  string
        output table file name

    Examples
    --------
    >>> import os
    >>> from mindboggle.mio.tables import select_column_from_mindboggle_tables
    >>> path = os.environ['MINDBOGGLE_DATA'] # doctest: +SKIP
    >>> subject1 = os.path.basename(path) # doctest: +SKIP
    >>> subject2 = os.path.basename(path) # doctest: +SKIP
    >>> subjects = [subject1, subject2] # doctest: +SKIP
    >>> hemi = 'left'
    >>> index = 2
    >>> tables_dir = os.path.dirname(path) # doctest: +SKIP
    >>> table_name = "label_shapes.csv"
    >>> label_name = 'Label name'
    >>> is_surface_table = True
    >>> write_table = True
    >>> output_table = ''
    >>> tables, cols, output = select_column_from_mindboggle_tables(subjects,
    ...     hemi, index, tables_dir, table_name, is_surface_table,
    ...     write_table, output_table) # doctest: +SKIP
    >>> cols[0][0] # doctest: +SKIP
    878.03969839999979
    >>> cols[0][1] # doctest: +SKIP
    3085.6236725000008
    >>> cols[0][2] # doctest: +SKIP
    1761.2330760000002

    """
    import os

    from mindboggle.mio.tables import select_column_from_tables

    # ------------------------------------------------------------------------
    # Construct list of Mindboggle shape table file names:
    # ------------------------------------------------------------------------
    tables = []
    for subject in subjects:
        if is_surface_table:
            table = os.path.join(tables_dir, subject, 'tables',
                                 hemi+'_cortical_surface', table_name)
        else:
            table = os.path.join(tables_dir, subject, 'tables', table_name)
        tables.append(table)

    # ------------------------------------------------------------------------
    # Extract columns and construct new table:
    # ------------------------------------------------------------------------
    tables, columns, output_table = select_column_from_tables(tables, index,
        write_table, output_table)

    return tables, columns, output_table


def explode_mindboggle_tables(subject_path='', output_path='',
                              break_column='label ID', verbose=False):
    """
    Given a subject name corresponding to Mindboggle outputs, break up each
    hemisphere's vertices.csv file into a separate table file for each label.

    Parameters
    ----------
    subject_path : string
        path to subject directory
    output_path : string
        output path/directory
    break_column : string
        column header that contains the integers to break up into tables
    verbose : bool
        print statements?

    Examples
    --------
    >>> import os
    >>> from mindboggle.mio.tables import explode_mindboggle_tables
    >>> subject_path = os.environ['MINDBOGGLE_DATA'] # doctest: +SKIP
    >>> output_path = '.'
    >>> break_column = 'label ID'
    >>> verbose = False
    >>> explode_mindboggle_tables(subject_path, output_path, break_column,
    ...                           verbose) # doctest: +SKIP
    """
    import os
    import numpy as np
    import pandas as pd

    if not os.path.exists(output_path):
        if verbose:
            print("{0} does not exist".format(output_path))
    else:

        for side in ['left', 'right']:

            output_dir = os.path.join(output_path, side + '_exploded_tables')
            if not os.path.exists(output_dir):
                if verbose:
                    print("Create missing output directory: {0}".
                          format(output_dir))
                os.mkdir(output_dir)
            if os.path.exists(output_dir):

                vertices_table = os.path.join(subject_path, 'tables',
                                              side + '_cortical_surface',
                                              'vertices.csv')
                shape_column_names = ['travel depth',
                                      'geodesic depth',
                                      'mean curvature',
                                      'freesurfer curvature',
                                      'freesurfer thickness',
                                      'freesurfer convexity (sulc)']
                shape_column_names2 = ['travel_depth',
                                       'geodesic_depth',
                                       'mean_curvature',
                                       'freesurfer_curvature',
                                       'freesurfer_thickness',
                                       'freesurfer_sulc']

                if verbose:
                    print("Explode {0} by {1} values".
                        format(vertices_table, break_column))

                df = pd.read_csv(vertices_table,
                        header=0,
                        index_col=break_column)
                        #usecols=shape_column_names)

                df1 = df[shape_column_names]
                unique_labels = np.unique(df1.index)

                df1.columns = shape_column_names2

                for label in unique_labels:
                    label_table = df1.loc[label]

                    out_file = os.path.join(output_dir, str(label) + '.csv')
                    label_table.to_csv(out_file, index=False)

                    if not os.path.exists(out_file):
                        raise IOError(out_file + " not found")
            else:
                raise IOError('Unable to make directory {0}'.
                              format(output_dir))


# ============================================================================
# Doctests
# ============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)  # py.test --doctest-modules
