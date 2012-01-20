#!/bin/bash

time python extract.py -S $1/rh.inflated -T $1/rh.thickness -C $1/rh.curv $1/rh.sulc $1/rh.pial
#time python extract.py -S $1/rh.inflated  -T $1/rh.thickness $1/rh.sulc $1/rh.pial 
