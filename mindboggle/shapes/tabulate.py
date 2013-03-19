#!/usr/bin/env python
"""
Functions for writing tables


Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def write_mean_shapes_tables(labels, subfolds=[], fundi=[], sulci=[],
                             area_file='', depth_file='',
                             mean_curvature_file='', gauss_curvature_file='',
                             max_curvature_file='', min_curvature_file='',
                             thickness_file='', convexity_file='',
                             exclude_labels=[-1]):
    """
    Make tables of mean shape values per label, subfold, fundus, and/or sulcus.

    Parameters
    ----------
    labels :  list of integers
        indices to labels, one per vertex, with -1 indicating no label
    subfolds :  list of integers
        indices to subfolds, one per vertex, with -1 indicating no subfold
    fundi :  list of integers
        indices to fundi, one per vertex, with -1 indicating no fundus
    sulci :  list of integers
        indices to sulci, one per vertex, with -1 indicating no sulcus
    area_file :  string
        name of VTK file with scalar surface area values
    depth_file :  string
        name of VTK file with scalar depth values
    mean_curvature_file :  string
        name of VTK file with scalar mean curvature values
    gauss_curvature_file :  string
        name of VTK file with scalar Gaussian curvature values
    max_curvature_file :  string
        name of VTK file with scalar maximum curvature values
    min_curvature_file :  string
        name of VTK file with scalar minimum curvature values
    thickness_file :  string
        name of VTK file with scalar thickness values
    convexity_file :  string
        name of VTK file with scalar convexity values
    exclude_labels : list of integers
        indices to be excluded (in addition to -1)

    Returns
    -------
    label_table :  string
        output table filename for label shapes
    subfold_table :  string
        output table filename for subfold shapes
    fundus_table :  string
        output table filename for fundus shapes
    sulcus_table :  string
        output table filename for sulcus shapes
    norm_label_table :  string
        output table filename for label shapes normalized by area
    norm_subfold_table :  string
        output table filename for subfold shapes normalized by area
    norm_fundus_table :  string
        output table filename for fundus shapes normalized by area
    norm_sulcus_table :  string
        output table filename for sulcus shapes normalized by area

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.shapes.tabulate import write_mean_shapes_tables
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> subfolds_file = os.path.join(path, 'arno', 'features', 'subfolds.vtk')
    >>> fundi_file = os.path.join(path, 'arno', 'features', 'fundi.vtk')
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> labels, name = read_scalars(labels_file)
    >>> subfolds, name = read_scalars(subfolds_file)
    >>> fundi, name = read_scalars(fundi_file)
    >>> sulci, name = read_scalars(sulci_file)
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> mean_curvature_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.avg.vtk')
    >>> gauss_curvature_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.gauss.vtk')
    >>> max_curvature_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.max.vtk')
    >>> min_curvature_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.min.vtk')
    >>> thickness_file = ''
    >>> convexity_file = ''
    >>> exclude_labels = [-1]
    >>> #
    >>> write_mean_shapes_tables(labels, subfolds, fundi, sulci, \
    >>>     area_file, depth_file, mean_curvature_file, \
    >>>     gauss_curvature_file, max_curvature_file, min_curvature_file, \
    >>>     thickness_file, convexity_file, exclude_labels)

    """
    import os
    from mindboggle.shapes.measure import mean_value_per_label
    from mindboggle.utils.io_vtk import read_scalars
    from mindboggle.utils.io_file import write_columns

    # Feature lists:
    feature_lists = [labels, subfolds, fundi, sulci]
    table_names = ['label_shapes.csv', 'subfold_shapes.csv',
                   'fundus_shapes.csv', 'sulcus_shapes.csv']

    # Shape names corresponding to shape files below:
    shape_names = ['area', 'depth', 'mean_curvature', 'gauss_curvature',
                   'max_curvature', 'min_curvature', 'thickness', 'convexity']

    # Load shape files as a list of numpy arrays of per-vertex shape values:
    shape_files = [area_file, depth_file, mean_curvature_file,
                   gauss_curvature_file, max_curvature_file, min_curvature_file,
                   thickness_file, convexity_file]
    shape_arrays = []
    column_names = []
    normalize_by_area = False
    for ishape, shape_file in enumerate(shape_files):
        if os.path.exists(shape_file):
            scalars_array, name = read_scalars(shape_file, True, True)
            if scalars_array.size:
                shape_arrays.append(scalars_array)
                column_names.append(shape_names[ishape])

                # Store area array:
                if ishape == 0:
                    normalize_by_area = True
                    area_array = scalars_array.copy()

    # Initialize table file names:
    label_table = None
    subfold_table = None
    fundus_table = None
    sulcus_table = None
    norm_label_table = None
    norm_subfold_table = None
    norm_fundus_table = None
    norm_sulcus_table = None

    # Loop through features / tables:
    for itable, feature_list in enumerate(feature_lists):
        print('itbl: {0} ftrlist: {1}'.format(itable, type(feature_list)))
        if feature_list:

            # Loop through shape measures:
            columns = []
            norm_columns = []
            for shape_array in shape_arrays:

                # Compute mean shape value per feature:
                mean_values, label_list, surface_areas, \
                norm_mean_values = mean_value_per_label(shape_array,
                    feature_list, exclude_labels, normalize_by_area, area_array)

                columns.append(mean_values)
                if normalize_by_area:
                    norm_columns.append(norm_mean_values)

            # For each feature, construct a table of average shape values:
            if columns:

                # Define output table:
                table_file = os.path.join(os.getcwd(), table_names[itable])

                # Write labels to table:
                write_columns(label_list, 'label', table_file)

                # Append columns of average shape values to table:
                write_columns(columns, column_names, table_file, table_file)

                # Write labels and normalized average shape values to table:
                if normalize_by_area:
                    norm_table_file = os.path.join(os.getcwd(),
                        'norm_' + table_names[itable])
                    write_columns(label_list, 'label', norm_table_file)
                    write_columns(norm_columns, column_names,
                                  norm_table_file, norm_table_file)

                # Return correct table file name:
                if itable == 0:
                    label_table = table_file
                elif itable == 1:
                    subfold_table = table_file
                elif itable == 2:
                    fundus_table = table_file
                elif itable == 3:
                    sulcus_table = table_file
                if normalize_by_area:
                    if itable == 0:
                        norm_label_table = norm_table_file
                    elif itable == 1:
                        norm_subfold_table = norm_table_file
                    elif itable == 2:
                        norm_fundus_table = norm_table_file
                    elif itable == 3:
                        norm_sulcus_table = norm_table_file

    return label_table, subfold_table, fundus_table, sulcus_table, \
           norm_label_table, norm_subfold_table, \
           norm_fundus_table, norm_sulcus_table


def write_vertex_shapes_table(table_file,
                              labels, subfolds=[], fundi=[], sulci=[],
                              area_file='', depth_file='', depth_rescaled_file='',
                              mean_curvature_file='', gauss_curvature_file='',
                              max_curvature_file='', min_curvature_file='',
                              thickness_file='', convexity_file=''):
    """
    Make a table of shape values per vertex.

    Parameters
    ----------
    table_file : output filename (without path)
    labels :  list of integers
        indices to labels, one per vertex, with -1 indicating no label
    subfolds :  list of integers
        indices to subfolds, one per vertex, with -1 indicating no subfold
    fundi :  list of integers
        indices to fundi, one per vertex, with -1 indicating no fundus
    sulci :  list of integers
        indices to sulci, one per vertex, with -1 indicating no sulcus
    area_file :  string
        name of VTK file with scalar surface area values
    depth_file :  string
        name of VTK file with scalar depth values
    depth_rescaled_file :  string
        name of VTK file with scalar rescaled depth values
    mean_curvature_file :  string
        name of VTK file with scalar mean curvature values
    gauss_curvature_file :  string
        name of VTK file with scalar Gaussian curvature values
    max_curvature_file :  string
        name of VTK file with scalar maximum curvature values
    min_curvature_file :  string
        name of VTK file with scalar minimum curvature values
    thickness_file :  string
        name of VTK file with scalar thickness values
    convexity_file :  string
        name of VTK file with scalar convexity values

    Returns
    -------
    shape_table : table file name for vertex shape values

    Examples
    --------
    >>> import os
    >>> from mindboggle.utils.io_vtk import read_scalars
    >>> from mindboggle.shapes.tabulate import write_vertex_shapes_table
    >>> #
    >>> table_file = 'vertex_shapes.csv'
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> subfolds_file = os.path.join(path, 'arno', 'features', 'subfolds.vtk')
    >>> fundi_file = os.path.join(path, 'arno', 'features', 'fundi.vtk')
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> labels, name = read_scalars(labels_file)
    >>> subfolds, name = read_scalars(subfolds_file)
    >>> fundi, name = read_scalars(fundi_file)
    >>> sulci, name = read_scalars(sulci_file)
    >>> area_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.area.vtk')
    >>> depth_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.depth.vtk')
    >>> depth_rescaled_file = os.path.join(path, 'arno', 'shapes', 'depth_rescaled.vtk')
    >>> mean_curvature_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.avg.vtk')
    >>> gauss_curvature_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.gauss.vtk')
    >>> max_curvature_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.max.vtk')
    >>> min_curvature_file = os.path.join(path, 'arno', 'shapes', 'lh.pial.curv.min.vtk')
    >>> thickness_file = ''
    >>> convexity_file = ''
    >>> #
    >>> write_vertex_shapes_table(table_file, labels, subfolds, fundi, sulci, \
    >>>     area_file, depth_file, depth_rescaled_file, mean_curvature_file, \
    >>>     gauss_curvature_file, max_curvature_file, min_curvature_file, \
    >>>     thickness_file, convexity_file)

    """
    import os
    from mindboggle.utils.io_vtk import read_scalars
    from mindboggle.utils.io_file import write_columns

    # Feature names and corresponding feature lists:
    feature_names = ['label', 'subfold', 'fundus', 'sulcus']
    feature_lists = [labels, subfolds, fundi, sulci]

    # Shape names and corresponding shape files:
    shape_names = ['area', 'depth', 'depth_rescaled', 'mean_curvature',
                   'gauss_curvature', 'max_curvature', 'min_curvature',
                   'thickness', 'convexity']
    shape_files = [area_file, depth_file, depth_rescaled_file,
                   mean_curvature_file, gauss_curvature_file,
                   max_curvature_file, min_curvature_file,
                   thickness_file, convexity_file]

    # Append columns of per-vertex scalar values:
    columns = []
    column_names = []
    for ifeature, values in enumerate(feature_lists):
        print('iftr: {0} values: {1} columns: {2}'.format(ifeature, type(values), type(columns)))
        if values:
            if not columns:
                indices = range(len(values))
            columns.append(values)
            column_names.append(feature_names[ifeature])

    for ishape, shape_file in enumerate(shape_files):
        if os.path.exists(shape_file):
            scalars, name = read_scalars(shape_file)
            if len(scalars):
                if not columns:
                    indices = range(len(scalars))
                columns.append(scalars)
                column_names.append(shape_names[ishape])

    # Prepend with column of indices and write table
    shapes_table = os.path.join(os.getcwd(), table_file)
    write_columns(indices, 'index', shapes_table)
    write_columns(columns, column_names, shapes_table, shapes_table)

    return shapes_table

