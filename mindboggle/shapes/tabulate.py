#!/usr/bin/env python
"""
Functions for writing tables


Authors:
    - Arno Klein, 2012-2013  (arno@mindboggle.info)  http://binarybottle.com

Copyright 2013,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

"""

def write_mean_shapes_table(table_file, column_names, labels, depth_file,
                            depth_rescaled_file, mean_curvature_file,
                            gauss_curvature_file,
                            max_curvature_file, min_curvature_file,
                            thickness_file='', convexity_file='',
                            norm_vtk_file='', exclude_labels=[]):
    """
    Make a table of mean values per label per measure.

    Parameters
    ----------
    filename :  output filename (without path)
    column_names :  names of columns [list of strings]
    labels :  name of label file or list of labels (same length as values)
    *shape_files :  vtk files with scalar values
    norm_vtk_file :  name of file containing per-vertex normalization values
                     (e.g., surface areas)
    exclude_labels : list of integer labels to be excluded

    Returns
    -------
    means_file : table file name for mean shape values
    norm_means_file : table file name for mean shape values normalized by area

    Examples
    --------
    >>> import os
    >>> from mindboggle.shapes.tabulate import write_mean_shapes_table
    >>> table_file = 'test_write_mean_shapes_table.txt'
    >>> column_names = ['labels', 'area', 'depth', 'depth_rescaled', 'mean_curvature',
    >>>                 'gauss_curvature', 'max_curvature', 'min_curvature',
    >>>                 'thickness', 'convexity']
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> labels_file = os.path.join(path, 'arno', 'labels', 'lh.labels.DKT25.manual.vtk')
    >>> exclude_values = [-1]
    >>> area_file = os.path.join(path, 'arno', 'measures', 'lh.pial.area.vtk')
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> depth_rescaled_file = os.path.join(path, 'arno', 'measures', 'depth_rescaled.vtk')
    >>> mean_curv_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.avg.vtk')
    >>> gauss_curv_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.gauss.vtk')
    >>> max_curv_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.max.vtk')
    >>> min_curv_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.min.vtk')
    >>> #
    >>> write_mean_shapes_table(table_file, column_names, labels_file,
    >>>                         depth_file, depth_rescaled_file, mean_curv_file,
    >>>                         gauss_curv_file, max_curv_file, min_curv_file,
    >>>                         thickness_file='', convexity_file='',
    >>>                         norm_vtk_file=area_file, exclude_labels=exclude_values)

    """
    import os
    import numpy as np
    from mindboggle.shapes.measure import mean_value_per_label
    from mindboggle.utils.io_vtk import read_scalars
    from mindboggle.utils.io_file import write_columns

    # Load per-vertex labels and normalization vtk files
    if type(labels) == str:
        labels, name = read_scalars(labels, return_first=True, return_array=True)
    if norm_vtk_file:
        norms, name = read_scalars(norm_vtk_file, return_first=True, return_array=True)
    else:
        norms = np.ones(len(labels))

    # List files
    vtk_files = [depth_file, depth_rescaled_file, mean_curvature_file,
                 gauss_curvature_file, max_curvature_file, min_curvature_file,
                 thickness_file, convexity_file]

    # Append columns of values to table
    columns = []
    norm_columns = []
    for i, vtk_file in enumerate(vtk_files):
        values, name = read_scalars(vtk_file, return_first=True, return_array=True)
        if values.size:
            mean_values, norm_mean_values, norm_values, \
                label_list = mean_value_per_label(values, norms, labels,
                                                  exclude_labels)
            columns.append(mean_values)
            norm_columns.append(norm_mean_values)
        else:
            column_names[i] = ''

    if columns:

        # Prepend with column of normalization values
        columns.insert(0, norm_values)
        norm_columns.insert(0, norm_values)
        column_names.insert(0, 'area')

        means_file = os.path.join(os.getcwd(), table_file)
        write_columns(label_list, 'label', means_file)
        write_columns(columns, column_names, means_file, means_file)

        norm_means_file = os.path.join(os.getcwd(), 'norm_' + table_file)
        write_columns(label_list, 'label', norm_means_file)
        write_columns(norm_columns, column_names, norm_means_file, norm_means_file)

    else:
        print("NOTE: 'columns' is empty. Nothing written.")

    return means_file, norm_means_file

def write_vertex_shapes_table(table_file, column_names,
                              labels_file, subfolds_file, sulci_file, fundi_file,
                              area_file, depth_file, depth_rescaled_file,
                              mean_curvature_file, gauss_curvature_file,
                              max_curvature_file, min_curvature_file,
                              thickness_file='', convexity_file=''):
    """
    Make a table of shape values per vertex.

    Parameters
    ----------
    table_file : output filename (without path)
    column_names : names of columns [list of strings]
    *vtk_files : arbitrary number of vtk files with per-vertex scalar values
                 (set each missing file to an empty string)

    Returns
    -------
    shape_table : table file name for vertex shape values

    Examples
    --------
    >>> import os
    >>> from mindboggle.shapes.tabulate import write_vertex_shapes_table
    >>>
    >>> path = os.environ['MINDBOGGLE_DATA']
    >>> table_file = 'test_write_vertex_shape_table.txt'
    >>> column_names = ['labels', 'subfolds', 'sulcus', 'fundus', 'area', 'depth',
    >>>                 'depth_rescaled', 'mean_curvature', 'gauss_curvature',
    >>>                 'max_curvature', 'min_curvature', '', '']
    >>> labels_file = ''
    >>> fundi_file = ''
    >>> subfolds_file = os.path.join(path, 'arno', 'features', 'subfolds.vtk')
    >>> sulci_file = os.path.join(path, 'arno', 'features', 'sulci.vtk')
    >>> area_file = os.path.join(path, 'arno', 'measures', 'lh.pial.area.vtk')
    >>> depth_file = os.path.join(path, 'arno', 'measures', 'lh.pial.depth.vtk')
    >>> depth_rescaled_file = os.path.join(path, 'arno', 'measures', 'depth_rescaled.vtk')
    >>> mean_curv_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.avg.vtk')
    >>> gauss_curv_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.gauss.vtk')
    >>> max_curv_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.max.vtk')
    >>> min_curv_file = os.path.join(path, 'arno', 'measures', 'lh.pial.curv.min.vtk')
    >>> #
    >>> write_vertex_shapes_table(table_file, column_names,
    >>>     labels_file, sulci_file, fundi_file, area_file, depth_file,
    >>>     depth_rescaled_file, mean_curv_file, gauss_curv_file,
    >>>     max_curv_file, min_curv_file, '', '')

    """
    import os
    from mindboggle.utils.io_vtk import read_scalars
    from mindboggle.utils.io_file import write_columns

    # List files
    vtk_files = [labels_file, subfolds_file, sulci_file, fundi_file,
                 area_file, depth_file, depth_rescaled_file,
                 mean_curvature_file, gauss_curvature_file,
                 max_curvature_file, min_curvature_file,
                 thickness_file, convexity_file]

    # Append columns of values to table
    columns = []
    for i, vtk_file in enumerate(vtk_files):
        if vtk_file:
            values, name = read_scalars(vtk_file)
            if values:
                if not columns:
                    indices = range(len(values))
                columns.append(values)
        else:
            column_names[i] = ''

    # Prepend with column of indices and write table
    shape_table = os.path.join(os.getcwd(), table_file)
    write_columns(indices, 'index', shape_table)
    write_columns(columns, column_names, shape_table, shape_table)

    return shape_table

