#!/bin/bash

#time python extract.py -S $1/lh.inflated -C $1/lh.curv $1/lh.sulc $1/lh.pial
#time python extract.py -S $1/lh.inflated   $1/lh.sulc $1/lh.pial
time python extract.py -S $1/lh.inflated -C $1/lh.curv -T $1/lh.thickness $1/lh.sulc $1/lh.pial