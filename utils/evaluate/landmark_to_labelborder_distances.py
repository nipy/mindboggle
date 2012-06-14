#!/usr/bin/python
"""
Measure the Hausdorff distance between landmarks and label boundaries for each subject.

(c) 2011, @rno klein
"""

import os
from os.path import exists
from subprocess import call

verbose = 1
dim = 3
source_files  = ['m1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m11','m12']
landmark_type = 'fundi_forrest_bao'
#landmark_type = 'fundi_gang_li'
#landmark_type = 'fundi_brain_visa'
#landmark_type = 'ribbons_brain_visa'

out_path    = '/hd2/Archive/registration_evaluation_2011_output/'
in_dir      = os.path.join( out_path, 'Transformed_Landmarks/')
ext         = '.nii.gz'
results_dir = os.path.join( out_path, 'Results/')

atlas_borders_dir  = '/hd2/Brains/CUMC12/Atlas_label_shells_aggregated_binary/'

results_file = results_dir+'avg_hausdorff_distance_to_label_borders_'+landmark_type+'.txt'
f = open(results_file, 'w');
count = 0
average_avg = 0
average_std = 0
average_max = 0
append = in_dir + 'distance_'
temp_file1 = results_dir+'temp2_hausdorff1.txt'
temp_file2 = results_dir+'temp2_hausdorff2.txt'
print(landmark_type)
for file in source_files:
    source = file+'_to_template_'+landmark_type+ext
    source_file = in_dir + source
    target = file+ext
    target_file = atlas_borders_dir + target
    if os.path.exists(source_file) and \
       os.path.exists(target_file):

        count += 1
        average_dist = 0
        
        """
        c3d m1.nii.gz -sdt -o dm1.nii.gz
        c3d m2.nii.gz -sdt -o dm2.nii.gz
        c3d dm2.nii.gz m1.nii.gz -multiply -o d2m1.nii.gz
        c3d dm1.nii.gz m2.nii.gz -multiply -o d1m2.nii.gz
        c3d d2m1.nii.gz m1.nii.gz -lstat
        LabelID        Mean        StdD         Max         Min       Count
            0       0.00000     0.00000     0.00000     0.00000     8119765
            1       3.77377     2.05970    13.85825     0.00000        6699
        c3d d1m2.nii.gz m2.nii.gz -lstat
        LabelID        Mean        StdD         Max         Min       Count
            0       0.00000     0.00000     0.00000     0.00000     8120521
            1       3.45281     1.81034    14.03562     0.00000        5943
        """
        args1 = " ".join(['c3d', source_file, '-sdt -o', append+source])
        args2 = " ".join(['c3d', target_file, '-sdt -o', append+target])
        args3 = " ".join(['c3d', append+target, source_file, \
                         '-multiply -o', append+'from_'+source])
        args4 = " ".join(['c3d', append+source, target_file, \
                         '-multiply -o', append+'from_'+target])
        args5 = " ".join(['c3d', append+'from_'+target, \
                                 target_file, '-label-statistics > ', temp_file1])
        args6 = " ".join(['c3d', append+'from_'+source, \
                                 source_file, '-label-statistics > ', temp_file2])
        """
        print(args1); 
        print(args2); 
        print(args3); 
        print(args4); 
        print(args5); 
        print(args6); 
        """
#        p = call(args1, shell="True")
#        p = call(args2, shell="True")
        p = call(args3, shell="True")
        p = call(args4, shell="True")
        p = call(args5, shell="True")
        p = call(args6, shell="True")

        f1 = open(temp_file1,'r')
        f2 = open(temp_file2,'r')
        temp1 = f1.read()
        temp2 = f2.read()
        results1 = [0,0,0,0]
        results2 = [0,0,0,0]
        if temp1 != '':
            results1 = temp1.split("\n")[2].split()
        if temp2 != '':
            results2 = temp2.split("\n")[2].split()
        avg_results = (float(results1[1]) + float(results2[1]))/2.0
        std_results = (float(results1[2]) + float(results2[2]))/2.0
        max_results = (float(results1[3]) + float(results2[3]))/2.0
        print_out = " ".join([file, '(mean, std, max):', \
                              str(avg_results), str(std_results), str(max_results)])
        print(print_out)
        f.close()
        f = open(results_file, 'a')
        f.write(print_out + '\n')
        average_avg += avg_results
        average_std += std_results
        average_max += max_results
    else:
        if not os.path.exists(source_landmarks):
            raise NameError('Check input file ' + source_landmarks)
        elif not os.path.exists(target_landmarks):
            raise NameError('Check input file ' + target_landmarks)

average_avg = average_avg/count
average_std = average_std/count
average_max = average_max/count

print_out = " ".join(['Average (mean, std, max) Hausdorff distance to label borders:', \
            str(average_avg), str(average_std), str(average_max)])
print(print_out);
f.close()
f = open(results_file, 'a');
f.write(print_out)
f.close()
