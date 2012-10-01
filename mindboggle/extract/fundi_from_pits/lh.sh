#!/bin/bash

# Usage example:  sh lh.sh ~/test_subject/surf
time python vtk_extract.py --thick $1/lh.thickness --convex $1/lh.sulc --fundi $1/lh.fundi.vtk --pits $1/lh.pits.vtk --sulci $1/lh.sulci.vtk  --sulciThld 0.15  --meancurv $1/lh.mean_curv.vtk --gausscurv $1/lh.gauss_curv.vtk $1/lh.depth.vtk  --clouchoux
