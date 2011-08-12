Usage : mesh_sulcalparcel_usingpit -from_max/-from_min -fs/-mni a.asc/a.obj field.map pit.txt final_value save_region.txt

If you use freesurfer surface, you should do '-fs  lh.smoothwm.asc' The freesurfer surface should be ascii format.
field.map is the feature map file (maybe, depth file, lh.sulc.f10.1D.dset) that was used for sulcal pit extraction.
watershed progress from_max: max value to final_value (decrease), -from_min: min value to final_value (increase)
final_value: do watershed until final_value is met
