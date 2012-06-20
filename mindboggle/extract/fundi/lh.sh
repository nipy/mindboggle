#!/bin/bash

#time python extract.py -S $1/lh.inflated -C $1/lh.curv -T $1/lh.thickness $1/lh.sulc $1/lh.pial
time python vtk_extract.py -S $1/lh.inflated -C $1/lh.sulc -T $1/lh.thickness $1/lh.euclidean.vtk
#time python vtk_extract.py $1/surf/lh.euclidean.vtk
