#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import argparse
import numpy as np

from mindboggle.mio.colors import distinguishable_colors, label_adjacency_matrix


if __name__ == "__main__":

    description = ('calculate colormap for labeled image;'
                   'calculated result is stored in output_dirname/colors.npy')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('label_filename', help='path to the label image')
    parser.add_argument('output_dirname', help='path to the folder storing '
                                               'temporary files and result')
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    args = parser.parse_args()

    if not os.path.isdir(args.output_dirname):
        os.makedirs(args.output_dirname)
        
    matrix_filename = os.path.join(args.output_dirname, 'matrix.npy')
    colormap_filename = os.path.join(args.output_dirname, 'colormap.npy')
    labels_filename = os.path.join(args.output_dirname, 'labels.npy')
    colors_filename = os.path.join(args.output_dirname, 'colors.npy')

    if args.verbose:
        print('finding adjacency maps...')

    if not os.path.isfile(matrix_filename) or \
            not os.path.isfile(labels_filename):
        labels, matrix = label_adjacency_matrix(args.label_filename,
                                                out_dir=args.output_dirname)[:2]
        matrix = matrix.as_matrix()[:, 1:]
        np.save(matrix_filename, matrix)
        np.save(labels_filename, labels)
    else:
        labels = np.load(labels_filename)
        matrix = np.load(matrix_filename)

    if args.verbose:
        print('finding colormap...')

    if not os.path.isfile(colormap_filename):
        num_colors = len(labels)
        colormap = distinguishable_colors(ncolors=num_colors,
                                          plot_colormap=False,
                                          save_csv=False,
                                          out_dir=args.output_dirname)
        np.save(colormap_filename, colormap)
    else:
        colormap = np.load(colormap_filename)

    if args.verbose:
        print('finding label colors')

    if not os.path.isfile(colors_filename):
        label_colors = colors.group_colors(colormap, 
                                           args.label_filename, 
                                           IDs=labels,
                                           adjacency_matrix=matrix, 
                                           out_dir=args.output_dirname,
                                           plot_colors=False,
                                           plot_graphs=False)
        np.save(colors_filename, label_colors)
