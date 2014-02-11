#!/bin/bash

# Usage example:  sh rh.sh ~/test_subject/surf

# for vtk_extract.py using Clouchoux-style pits
time python vtk_extract.py --thick $1/rh.thickness --convex $1/rh.sulc --fundi $1/rh.fundi.vtk --pits $1/rh.pits.vtk --sulci $1/rh.sulci.vtk  --sulciThld 0.15  --meancurv $1/rh.mean_curv.vtk --gausscurv $1/rh.gauss_curv.vtk --clouchoux $1/rh.depth.vtk

# for vtk_extract.py without Clouchoux-style pits
time python vtk_extract.py --thick $1/rh.thickness --convex $1/rh.sulc --fundi $1/rh.fundi.vtk --pits $1/rh.pits.vtk --sulci $1/rh.sulci.vtk  --sulciThld 0.15  --meancurv $1/rh.mean_curv.vtk --gausscurv $1/rh.gauss_curv.vtk $1/rh.depth.vtk

# for fs_extract.py
time python fs_extract.py --use conv --thick $1/rh.thickness --conv $1/rh.sulc --curv $1/rh.curv --fundi $1/rh.fundi.vtk --pits $1/rh.pits.vtk --sulci $1/rh.sulci.vtk --sulciThld 0.15 $1/rh.pial
