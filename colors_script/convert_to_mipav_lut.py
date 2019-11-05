#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse

def parse_inputs():
    des = 'convert colormap to mipav lut file.'
    parser = argparse.ArgumentParser(description=des)
    parser.add_argument('colormap_filename')
    parser.add_argument('ids_filename')
    parser.add_argument('output_filename')
    args = parser.parse_args()
    return args

def main(args):

    labels = np.load(args.ids_filename)
    colors = np.load(args.colormap_filename) * 255

    contents = list()
    contents.append('<LUT>')
    contents.append('256\t# Size of LUT Arrays')

    for i in range(256):
        if i in labels:
            idx = np.where(labels == i)[0][0]
            c = colors[idx].tolist()
        else:
            c = [0.0, 0.0, 0.0]
        c_str = [str(cc) for cc in c]
        line = '\t'.join([str(i), '1.0', *c_str])
        contents.append(line)

    with open(args.output_filename, 'w') as file:
        file.write('\n'.join(contents))


if __name__ == "__main__":
    args = parse_inputs()
    main(args)
