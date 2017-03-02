#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
sys.path.append('../mindboggle/mio')
import colors
import numpy as np
import os
import argparse
import shutil


def parse_inputs():

    des = 'calculate colormap for labeled image, store the colors in the file'\
        + ' colors.npy'
    parser = argparse.ArgumentParser(description=des)
    parser.add_argument('label_filename', help='the path of the label image file')
    help = 'path of folder to store intermediate files and final results'
    parser.add_argument('output_dir', help=help)
    args = parser.parse_args()
    return args


def main(args):

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
        
    matrix_filename = os.path.join(args.output_dir, 'matrix.npy')
    colormap_filename = os.path.join(args.output_dir, 'colormap.npy')
    IDs_filename = os.path.join(args.output_dir, 'IDs.npy')
    colors_filename = os.path.join(args.output_dir, 'colors.npy')

    print('finding adjacency maps')
    if not os.path.isfile(matrix_filename) and not os.path.isfile(IDs_filename):
        IDs, matrix, output = colors.label_adjacency_matrix(args.label_filename,
                                                            out_dir=args.output_dir)
        matrix = matrix.as_matrix()[:, 1:]
        np.save(matrix_filename, matrix)
        np.save(IDs_filename, IDs)
    else:
        IDs = np.load(IDs_filename)
        matrix = np.load(matrix_filename)

    print('finding colormap')
    if not os.path.isfile(colormap_filename):
        ncolors = len(IDs)
        colormap = colors.distinguishable_colors(ncolors=ncolors,
                                                 plot_colormap=False,
                                                 save_csv=False,
                                                 out_dir=args.output_dir)
        np.save(colormap_filename, colormap)
    else:
        colormap = np.load(colormap_filename)

    print('finding label colors')
    label_colors = colors.group_colors(colormap, args.label_filename, IDs=IDs,
                                       adjacency_matrix=matrix, 
                                       out_dir=args.output_dir,
                                       plot_colors=True,
                                       plot_graphs=False)
    np.save(colors_filename, label_colors)


if __name__ == "__main__":
    args = parse_inputs()
    main(args)
