#!/usr/bin/python

#mris_convert lh.smoothwm lh.smoothwm.asc
#mris_convert -c lh.sulc lh.smoothwm lh.sulc.asc
#fstex2mnitex lh.sulc.asc lh.sulc.dset
#SurfSmooth -i_fs lh.smoothwm -input lh.sulc.dset -met HEAT_07 -fwhm 10 -output lh.sulc.f10.1D.dset
#mesh_watershed_fs -from_max -merge 20 -all lh.smoothwm.asc lh.sulc.f10.1D.dset -5 lh.sulc.f10ma20.wts.dset lh.sulc.f10ma20.wss.dset lh.sulc.f10ma20.wss.label
#mris_distance_transform lh.smoothwm lh.sulc.f10ma20.wss.label unsigned lh.sulc.f10ma20.wss.vv
#	>> output lh.sulc.f10ma20.wss.vv.w
#mris_convert lh.sulc.f10ma20.wss.vv.w lh.sulc.f10ma20.wss.vv.asc
#merge_distance_seeds_fsdistmap -max -wss lh.smoothwm.asc lh.sulc.f10.1D.dset lh.sulc.f10ma20.wss.vv.asc lh.sulc.f10ma20.wss.dset 5 lh.sulc.f5ma5md5.wss.dset


import sys
from subprocess import *

if sys.argv[1]=='help':
        print 'Usage: sulcal_pit_fs.py <-diffusion_fwhm> <-merge_area> <-merge_distance> lh.smoothwm value.txt(curv or depth, lh.sulc) <final_value>'
        sys.exit(0)

diffusion_opt=sys.argv[1]
merge_area_opt=sys.argv[2]
merge_dist_opt=sys.argv[3]

input_surf=sys.argv[4]
valuefile=sys.argv[5]
final_value=sys.argv[6]

surf_ascii=input_surf+'.asc'
value_ascii=valuefile+'.asc'
mnivalue=value_ascii.replace('asc','dset')
smooth_value=mnivalue.replace('dset', 'f'+diffusion_opt+'.1D.dset')
watershed_seeds=smooth_value.replace('.1D','ma'+merge_area_opt+'.wss')
merge_seeds=smooth_value.replace('.1D','ma'+merge_area_opt+'md'+merge_dist_opt+'.wss')
watershed_regions=watershed_seeds.replace('wss','wts')
seeds_label=watershed_seeds.replace('dset','label')
seeds_distance=watershed_seeds.replace('dset','vv')
seeds_distance_out=seeds_distance+'.w'
seeds_distance_ascii=seeds_distance_out.replace('vv.w','vv.asc')


Popen(['mris_convert',input_surf,surf_ascii]).wait()
Popen(['mris_convert','-c', valuefile, input_surf, value_ascii]).wait()
Popen(['fstex2mnitex', value_ascii, mnivalue]).wait()
Popen(['SurfSmooth','-i_fs',input_surf, '-input', mnivalue, '-met', 'HEAT_07', '-fwhm', diffusion_opt, '-output',smooth_value]).wait()
Popen(['mesh_watershed_fs', '-from_max','-merge', merge_area_opt,'-all',surf_ascii,smooth_value,final_value,watershed_regions, watershed_seeds, seeds_label]).wait()
print 'Extract watershed regions and seeds\n'
Popen(['mris_distance_transform',input_surf,seeds_label,'unsigned',seeds_distance]).wait()
Popen(['mris_convert',seeds_distance_out,seeds_distance_ascii]).wait()
Popen(['merge_distance_seeds_fsdistmap','-max', '-wss', surf_ascii, smooth_value, seeds_distance_ascii, watershed_seeds, merge_dist_opt,merge_seeds]).wait()
print 'Merging with distance map on the surface\n'