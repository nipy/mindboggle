Developed by Kiho Im
Kiho.Im@childrens.harvard.edu
kiho.sky@gmail.com


Sulcal pit extraction using a FreeSurfer surface

Usage: sulcal_pit_fs.py <-diffusion_fwhm> <-merge_area> <-merge_distance> lh.smoothwm value.txt(depth, lh.sulc) <final_value>
ex) sulcal_pit_fs.py 10 20 5 lh.smoothwm lh.sulc -5


Parameters - (Im et al., Cereb. Cortex 2010)
-diffuseion_fwhm: fwhm for diffusion smoothing on surface
-merge_area: threshold of area in merging process 
-merge_distance: threshold of distance in merging process, CAUTION!: If you want use the threshold of 10 mm, put '5' (half value) value when running software.

final_value: Watershed runs until the depth of the vertex in the list is less or more than a threshold of final_value.


In sulcal_pit_fs.py
mris_convert, mris_distance_transform: FreeSurfer tool
SurfSmooth: AFNI tool
fstex2mnitex, mesh_watershed_fs, merge_distance_seeds_fsdistmap: Kiho's tool

mris_convert lh.smoothwm lh.smoothwm.asc
> convert binary form to ascii form

#mris_convert -c lh.sulc lh.smoothwm lh.sulc.asc
> format convert 

#fstex2mnitex lh.sulc.asc lh.sulc.dset
> format convert

#SurfSmooth -i_fs lh.smoothwm -input lh.sulc.dset -met HEAT_07 -fwhm 10 -output lh.sulc.f10.1D.dset
> heat kernel smoothing (Moo K. Chung) of depth map
> You can also use other surface smoothing tool, such as matlab software distributed by Chung (http://www.stat.wisc.edu/~mchung/softwares/software.html)

#mesh_watershed_fs -from_max -merge 20 -all lh.smoothwm.asc lh.sulc.f10.1D.dset -5 lh.sulc.f10ma20.wts.dset lh.sulc.f10ma20.wss.dset lh.sulc.f10ma20.wss.label
> Watershed with merging process
> More detail: mesh_watershed_fs help

#mris_distance_transform lh.smoothwm lh.sulc.f10ma20.wss.label unsigned lh.sulc.f10ma20.wss.vv
> For distance merging, calculate distance map. Sulcal pits are seeds.
> output lh.sulc.f10ma20.wss.vv.w

#mris_convert lh.sulc.f10ma20.wss.vv.w lh.sulc.f10ma20.wss.vv.asc
> format convert

#merge_distance_seeds_fsdistmap -max -wss lh.smoothwm.asc lh.sulc.f10.1D.dset lh.sulc.f10ma20.wss.vv.asc lh.sulc.f10ma20.wss.dset 5 lh.sulc.f5ma5md5.wss.dset
> Distance merging
> More detail: merge_distance_seeds_fsdistmap help 


Final sulcal pit map: lh.sulc.f5ma5md5.wss.dset

lh.sulc.f10ma20.wss.dset: sulcal pits before running distance merging
lh.sulc.f10ma20.wts.dset: sulcal catchment basin before running distance merging
