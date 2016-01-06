#!/usr/bin/env python
"""
Functions for creating tables.


Authors:
    - Arno Klein, 2012-2015  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2015,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

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
        exclude_labels=[-1]):
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
    affine_transform_files : list of strings
        affine transform files to standard space
    inverse_booleans : list of of zeros and ones
        for each transform, 1 to take the inverse, else 0
    transform_format : string
        format for transform file
        Ex: 'txt' for text, 'itk' for ITK, and 'mat' for Matlab format
    area_file :  string
        name of VTK file with surface area scalar values
    normalize_by_area : Boolean
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
    >>> from mindboggle.mio.vtks import read_scalars
    >>> from mindboggle.mio.tables import write_shape_stats
    >>> path = '/homedir/mindboggled/Twins-2-1'
    >>> labels_or_file = os.path.join(path, 'labels', 'left_cortical_surface', 'freesurfer_cortex_labels.vtk')
    >>> sulci_file = os.path.join(path, 'features', 'left_cortical_surface', 'sulci.vtk')
    >>> fundi_file = os.path.join(path, 'features', 'left_cortical_surface', 'fundus_per_sulcus.vtk')
    >>> sulci, name = read_scalars(sulci_file)
    >>> fundi, name = read_scalars(fundi_file)
    >>> affine_transform_files = [] #os.path.join(path, 'mri', 't1weighted_brain.MNI152Affine.txt')
    >>> inverse_booleans = []
    >>> #transform_format = 'mat'
    >>> transform_format = 'itk'
    >>> area_file = os.path.join(path, 'shapes', 'left_cortical_surface', 'area.vtk')
    >>> normalize_by_area = False
    >>> mean_curvature_file = os.path.join(path, 'shapes', 'left_cortical_surface', 'mean_curvature.vtk')
    >>> travel_depth_file = os.path.join(path, 'shapes', 'left_cortical_surface', 'travel_depth.vtk')
    >>> geodesic_depth_file = os.path.join(path, 'shapes', 'left_cortical_surface', 'geodesic_depth.vtk')
    >>> freesurfer_thickness_file = ''
    >>> freesurfer_curvature_file = ''
    >>> freesurfer_sulc_file = ''
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
    >>>     affine_transform_files, inverse_booleans, transform_format,
    >>>     area_file, normalize_by_area,
    >>>     mean_curvature_file, travel_depth_file, geodesic_depth_file,
    >>>     freesurfer_thickness_file, freesurfer_curvature_file,
    >>>     freesurfer_sulc_file,
    >>>     labels_spectra, labels_spectra_IDs,
    >>>     sulci_spectra, sulci_spectra_IDs,
    >>>     labels_zernike, labels_zernike_IDs,
    >>>     sulci_zernike, sulci_zernike_IDs,
    >>>     exclude_labels)

    """
    import os
    import numpy as np
    import pandas as pd

    from mindboggle.guts.compute import means_per_label, stats_per_label, \
        sum_per_label
    from mindboggle.mio.vtks import read_scalars, read_vtk, \
        apply_affine_transforms
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
        import sys
        sys.exit('No feature data to tabulate in write_shape_stats().')

    spectrum_start = 1  # Store all columns of spectral components (0),
                        # or start from higher frequency components (>=1)

    #-------------------------------------------------------------------------
    # Feature lists, shape names, and shape files:
    #-------------------------------------------------------------------------
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

        #-----------------------------------------------------------------
        # Label names:
        #-----------------------------------------------------------------
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

        #---------------------------------------------------------------------
        # For each feature, construct a table of average shape values:
        #---------------------------------------------------------------------
        if feature_list:
            feature_name = feature_names[itable]
            columns = []

            #-----------------------------------------------------------------
            # Loop through shape measures:
            #-----------------------------------------------------------------
            column_names.extend(column_names[:])
            for ishape, shape_array in enumerate(shape_arrays):
                shape = shape_names[ishape]
                print('  Compute statistics on {0} {1}...'.
                      format(feature_name, shape))
                #-------------------------------------------------------------
                # Append feature areas to columns:
                #-------------------------------------------------------------
                if ishape == 0 and np.size(area_array):
                    sums, label_list = sum_per_label(shape_array,
                        feature_list, include_labels, exclude_labels)
                    column_names.append(shape)
                    columns.append(sums)
                #-------------------------------------------------------------
                # Append feature shape statistics to columns:
                #-------------------------------------------------------------
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

            #-----------------------------------------------------------------
            # Mean positions in the original space:
            #-----------------------------------------------------------------
            # Compute mean position per feature:
            positions, sdevs, label_list, foo = means_per_label(points,
                feature_list, include_labels, exclude_labels, use_area)

            # Append mean x,y,z position per feature to columns:
            xyz_positions = np.asarray(positions)
            for ixyz, xyz in enumerate(['x','y','z']):
                column_names.append('mean position: {0}'.format(xyz))
                columns.append(xyz_positions[:, ixyz].tolist())

            #-----------------------------------------------------------------
            # Mean positions in standard space:
            #-----------------------------------------------------------------
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

            #-----------------------------------------------------------------
            # Laplace-Beltrami spectra:
            #-----------------------------------------------------------------
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

            #-----------------------------------------------------------------
            # Zernike moments:
            #-----------------------------------------------------------------
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

            #-----------------------------------------------------------------
            # Write labels/IDs and values to table:
            #-----------------------------------------------------------------
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
        affine_transform_files=[], inverse_booleans=[],
        transform_format='itk',
        area_file='', mean_curvature_file='', travel_depth_file='',
        geodesic_depth_file='', freesurfer_thickness_file='',
        freesurfer_curvature_file='', freesurfer_sulc_file=''):
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
    >>> #
    >>> output_table = ''#vertex_shapes.csv'
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> labels_or_file = os.path.join(path, 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> sulci_file = os.path.join(path, 'features', 'sulci.vtk')
    >>> fundi_file = os.path.join(path, 'features', 'fundi.vtk')
    >>> sulci, name = read_scalars(sulci_file)
    >>> fundi, name = read_scalars(fundi_file)
    >>> affine_transform_files = [os.path.join(path, 'mri',
    >>>     't1weighted_brain.MNI152Affine.txt')]
    >>> inverse_booleans = [1]
    >>> transform_format = 'itk'
    >>> swap_xy = True
    >>> area_file = os.path.join(path, 'shapes', 'lh.pial.area.vtk')
    >>> mean_curvature_file = os.path.join(path, 'shapes', 'lh.pial.mean_curvature.vtk')
    >>> travel_depth_file = os.path.join(path, 'shapes', 'lh.pial.travel_depth.vtk')
    >>> geodesic_depth_file = os.path.join(path, 'shapes', 'lh.pial.geodesic_depth.vtk')
    >>> freesurfer_thickness_file = ''
    >>> freesurfer_curvature_file = ''
    >>> freesurfer_sulc_file = ''
    >>> #
    >>> write_vertex_measures(output_table, labels_or_file, sulci, fundi,
    >>>     affine_transform_files, inverse_booleans, transform_format, area_file,
    >>>     mean_curvature_file, travel_depth_file, geodesic_depth_file,
    >>>     freesurfer_thickness_file, freesurfer_curvature_file, freesurfer_sulc_file)

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
        import sys
        sys.exit('No feature data to tabulate in write_vertex_measures().')

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
        raise(IOError(output_table + " not found"))

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
    >>> import os
    >>> from mindboggle.mio.tables import write_face_vertex_averages
    >>> path = '/homedir/mindboggled'
    >>> input_file = os.path.join(path, 'Twins-2-1', 'shapes', 'left_cortical_surface', 'freesurfer_thickness.vtk')
    >>> area_file = os.path.join(path, 'Twins-2-1', 'shapes', 'left_cortical_surface', 'area.vtk')
    >>> output_table = ''
    >>> #
    >>> write_face_vertex_averages(input_file, output_table, area_file)

    """
    import os
    import numpy as np
    import pandas as pd

    from mindboggle.mio.vtks import read_vtk, read_scalars

    points, indices, lines, faces, scalars, scalar_names, \
        npoints, input_vtk = read_vtk(input_file, True, True)
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

    df = pd.DataFrame({'': columns})
    df.to_csv(output_table, index=False)

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
    >>> from mindboggle.mio.tables import write_average_face_values_per_label
    >>> path = '/homedir/mindboggled'
    >>> input_indices_vtk = os.path.join(path, 'Twins-2-1', 'labels', 'left_cortical_surface', 'freesurfer_cortex_labels.vtk')
    >>> input_values_vtk = os.path.join(path, 'Twins-2-1', 'shapes', 'left_cortical_surface', 'freesurfer_thickness.vtk')
    >>> area_file = os.path.join(path, 'Twins-2-1', 'shapes', 'left_cortical_surface', 'area.vtk')
    >>> output_stem = 'labels_thickness'
    >>> exclude_values = [-1]
    >>> background_value = -1
    >>> #
    >>> write_average_face_values_per_label(input_indices_vtk,
    >>>     input_values_vtk, area_file, output_stem, exclude_values, background_value)
    >>> #
    >>> # View:
    >>> #example_vtk = os.path.join(os.getcwd(), output_stem + '0.vtk')
    >>> #from mindboggle.mio.plots import plot_surfaces
    >>> #plot_surfaces(example_vtk)

    """
    import os
    import numpy as np
    import pandas as pd

    from mindboggle.mio.vtks import read_scalars, read_vtk, write_vtk
    from mindboggle.guts.mesh import remove_faces

    # Load VTK file:
    points, indices, lines, faces, scalars, scalar_names, npoints, \
        input_vtk = read_vtk(input_indices_vtk, True, True)
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
        df = pd.DataFrame({'': columns})
        df.to_csv(output_table, index=False)

        # Write VTK file with scalar value:
        #output_vtk = os.path.join(os.getcwd(), output_stem + str(scalar) + '.vtk')
        #write_vtk(output_vtk, points, indices, lines, new_faces,
        #          [select_values.tolist()], [output_scalar_name])

        if not os.path.exists(output_table):
            raise(IOError(output_table + " not found"))


def select_column_from_tables(tables, index=0, write_table=True,
                              output_table=''):
    """
    Select column from list of tables, make a new table.

    Parameters
    ----------
    tables : list of strings
        table files (full paths)
    index : integer
        index for column to select
    write_table : Boolean
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
    >>> path = '/homedir/mindboggled'
    >>> tables = [os.path.join(path, 'Twins-2-1', 'tables', 'thickinthehead_per_freesurfer_cortex_label.csv'),
    >>>           os.path.join(path, 'Twins-2-1', 'tables', 'thickinthehead_per_freesurfer_cortex_label.csv')]
    >>> index = '2'
    >>> write_table = True
    >>> output_table = ''
    >>> select_column_from_tables(tables, index, write_table, output_table)

    """
    import os
    import pandas as pd
    import numpy as np

    #-------------------------------------------------------------------------
    # Construct a table:
    #-------------------------------------------------------------------------
    columns = []
    for input_table in tables:

        #---------------------------------------------------------------------
        # Extract column from the table for each subject:
        #---------------------------------------------------------------------
        if not os.path.exists(input_table):
            raise(IOError(input_table + " not found"))
        else:
            input_columns = pd.read_csv(input_table)
            columns.append(input_columns.iloc[:, index])

    #-------------------------------------------------------------------------
    # Write tables:
    #-------------------------------------------------------------------------
    if write_table and columns:
        if all([len(x) == len(columns[0]) for x in columns]):
            if not output_table:
                output_table = os.path.join(os.getcwd(),
                                            'select_column_from_tables.csv')
            df = pd.DataFrame({'': columns})
            df.to_csv(output_table, index=False)
        else:
            print('Not saving table.')

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
    is_surface_table : Boolean
        if True, use path to surface tables
    write_table : Boolean
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
    >>> subjects = ['Twins-2-1', 'Colin27-1']
    >>> hemi = 'left'
    >>> index = 2
    >>> tables_dir = os.path.join(os.environ['HOME'], 'mindboggled')
    >>> table_name = "label_shapes.csv"
    >>> label_name = 'Label name'
    >>> is_surface_table = True
    >>> write_table = True
    >>> output_table = ''
    >>> select_column_from_mindboggle_tables(subjects, hemi, index, tables_dir,
    >>>     table_name, is_surface_table, write_table, output_table)

    """
    import os

    from mindboggle.mio.tables import select_column_from_tables

    #-------------------------------------------------------------------------
    # Construct list of Mindboggle shape table file names:
    #-------------------------------------------------------------------------
    tables = []
    for subject in subjects:
        if is_surface_table:
            table = os.path.join(tables_dir, subject, 'tables',
                                 hemi+'_cortical_surface', table_name)
        else:
            table = os.path.join(tables_dir, subject, 'tables', table_name)
        tables.append(table)

    #-------------------------------------------------------------------------
    # Extract columns and construct new table:
    #-------------------------------------------------------------------------
    tables, columns, output_table = select_column_from_tables(tables, index,
        write_table, output_table)

    return tables, columns, output_table

def explode_mindboggle_tables(subject, subject_path='', output_path='',
                              break_column='label ID'):
    """
    Given a subject name corresponding to Mindboggle outputs, break up each
    hemisphere's vertices.csv file into a separate table file for each label.

    Parameters
    ----------
    subject : string
        name of subject run through Mindboggle
    subject_path : string
        path to subject's parent directory
    output_path : string
        output path/directory
    break_column : string
        column header that contains the integers to break up into tables

    Examples
    --------
    >>> from mindboggle.mio.tables import explode_mindboggle_tables
    >>> subject = 'Twins-2-1'
    >>> subject_path = '/Users/arno/mindboggled'
    >>> output_path = '/desk'
    >>> break_column = 'label ID'
    >>> explode_mindboggle_tables(subject, subject_path, output_path, break_column)
    """
    import os
    import numpy as np
    import pandas as pd

    if not os.path.exists(output_path):
        print("{0} does not exist".format(output_path))
    else:

        for side in ['left', 'right']:

            output_dir = os.path.join(output_path, side + '_exploded_tables')
            if not os.path.exists(output_dir):
                print("Create missing output directory: {0}".format(output_dir))
                os.mkdir(output_dir)
            if os.path.exists(output_dir):

                vertices_table = os.path.join(subject_path, subject, 'tables',
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
                        raise(IOError(out_file + " not found"))
            else:
                print('Unable to make directory {0}'.format(output_dir))


#=============================================================================
# Doctests
#=============================================================================
if __name__ == "__main__":
    import doctest
    doctest.testmod()
