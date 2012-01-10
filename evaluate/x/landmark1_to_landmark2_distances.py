#!/usr/bin/python
"""
Measure the average minimum distance from landmarks of one type (e.g., pits)
to landmarks of another type (e.g., fundi) within a subject.

(c) 2011, @rno klein
"""

import os
from os.path import exists
from subprocess import call

verbose = 1
dim = 3
files  = ['m1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m11','m12']
landmark_types1 = ['fundi_forrest_bao','fundi_gang_li','fundi_brain_visa']
landmark_types2 = ['pits_forrest_bao','pits_yrjo_hame','pits_kiho_im']
landmark_dirs1 = ['/hd2/Brains/CUMC12/Landmarks/fundi_forrest_bao_binary/',
                  '/hd2/Brains/CUMC12/Landmarks/fundi_gang_li_binary/',
                  '/hd2/Brains/CUMC12/Landmarks/fundi_brain_visa_binary/']
landmark_dirs2 = ['/hd2/Brains/CUMC12/Landmarks/pits_forrest_bao_binary/',
                  '/hd2/Brains/CUMC12/Landmarks/pits_yrjo_hame_binary/',
                  '/hd2/Brains/CUMC12/Landmarks/pits_kiho_im_binary/']
out_path = '/hd2/Archive/registration_evaluation_2011_output/'
dist_dir = out_path + 'Landmark_distances/'
results_dir = os.path.join( out_path, 'Results/')
ext = '.nii.gz'

count = len(files)
temp_file = results_dir+'temp_mindistance.txt'
for idir1,dir1 in enumerate(landmark_dirs1):
    landmark_type1 = landmark_types1[idir1]
    for idir2,dir2 in enumerate(landmark_dirs2):
        landmark_type2 = landmark_types2[idir2]
        results_file = results_dir+'avg_min_distance_'+landmark_type1+'_to_'+landmark_type2+'.txt'
        f = open(results_file, 'w')
        average_avg = 0
        average_std = 0
        average_max = 0
        for ifile,file in enumerate(files):
            source_file = dir1+file+ext
            dfile = dist_dir + 'distance_' + file + '_' + landmark_type1 + ext
            args = " ".join(['c3d', source_file, '-sdt -o', dfile])
            p = call(args, shell="True")
            target_file = dir2+file+ext
            if os.path.exists(source_file) and \
               os.path.exists(target_file):

                maskdist_file = dist_dir+'distance_from_'+file+'_'+landmark_type1+'_to_'+landmark_type2+ext
    
                """
                c3d m1.nii.gz -sdt -o dm1.nii.gz
                c3d dm1.nii.gz m2.nii.gz -multiply -o d1m2.nii.gz
                c3d d1m2.nii.gz m2.nii.gz -lstat
                LabelID        Mean        StdD         Max         Min       Count
                    0       0.00000     0.00000     0.00000     0.00000     8119765
                    1       3.77377     2.05970    13.85825     0.00000        6699
                """
                args2 = " ".join(['c3d', dfile, target_file, '-multiply -o', maskdist_file])
                args3 = " ".join(['c3d', maskdist_file, \
                                         target_file, '-label-statistics > ', temp_file])
                p = call(args2, shell="True")
                p = call(args3, shell="True")
    
                ftemp = open(temp_file,'r')
                temp = ftemp.read()
                results = [0,0,0,0]
                if temp != '':
                    results = temp.split("\n")[2].split()
                avg_results = float(results[1])
                std_results = float(results[2])
                max_results = float(results[3])
                print_out = " ".join([file,
                                      'avg_min_distance_'+landmark_type1+'_to_'+landmark_type2+'.txt', 
                                      '(mean, std, max):',
                                      str(avg_results), str(std_results), str(max_results)])
                print(print_out)
                f.close()
                f = open(results_file, 'a')
                f.write(print_out + '\n')
    
                average_avg += avg_results
                average_std += std_results
                average_max += max_results

            else:
                if not os.path.exists(source_file):
                    raise NameError('Check input file ' + source_file)
                elif not os.path.exists(target_file):
                    raise NameError('Check input file ' + target_file)

        average_avg = average_avg/count
        average_std = average_std/count
        average_max = average_max/count
        
        print_out = " ".join(['Average (mean, std, max) minimum distance:', \
                    str(average_avg), str(average_std), str(average_max)])
        print(print_out);
        f.close()
        f = open(results_file, 'a');
        f.write(print_out)
        f.close()
