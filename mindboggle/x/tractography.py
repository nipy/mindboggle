#!/usr/bin/env python
"""
Tractography starting from each in a list of labeled regions.

(c) 2013  Arno Klein  .  arno@mindboggle.info  .  www.binarybottle.com
    under Apache License Version 2.0
"""
import os, sys

#-------------------------------------------------------------------------------
# Paths:
#-------------------------------------------------------------------------------
mindboggle_dir = 'results'
subjects_dir = 'pre-post-treatment_predictors_fromTJ_MZ_20130327'
out_path = os.path.join(subjects_dir, 'tracks')
reg_dir = os.path.join(subjects_dir, 'registered')
dti_dir = os.path.join(subjects_dir, 'dti')
bedpostx_dir = os.path.join(subjects_dir, 'bedpostx')
brain_dir = os.path.join(subjects_dir, 'atropos')
dticoreg_dir = os.path.join(subjects_dir, 'dti_coreg')
subject_list = 'lists/pre-post-treatment_predictors.txt'
fid = open(subject_list, 'r')
subjects = fid.readlines()
subjects = [x.strip() for x in subjects]
subject_brains = []
for subject in subjects:
    new_file = os.path.join(brain_dir, subject+'_atroposoutput/brain_n3.nii.gz')
    subject_brains.append(new_file)

#-------------------------------------------------------------------------------
# Labels:
#-------------------------------------------------------------------------------
#labels = [1,2,4,5,6,101,102,104,105,106]
#label_names = ['R CA1', 'R subiculum', 'R CA4/dentate gyrus',
#    'R CA2/CA3', 'R stratum radiatum/lacunosum/moleculare',
#    'L CA1', 'L subiculum', 'L CA4/dentate gyrus',
#    'L CA2/CA3', 'L stratum radiatum/lacunosum/moleculare']
labels = [18,54]
label_names = ['L Amygdala','R Amygdala']
hippo = True
whole_hippo = False
if hippo:
    if whole_hippo:
        labels = [1,2]
        label_names = ['Hippocampus']
    else:
        labels = [4,104]
        label_names = ['R CA4/dentate gyrus','L CA4/dentate gyrus']

#-------------------------------------------------------------------------------
# Main program:
#-------------------------------------------------------------------------------
transform_labels = False
extract_label = False
transform_label = False
extract_dti_label = False
run_probtrackx = False
transform_probtrackx = True
run_camino_prep = False
run_camino_track = False
do_these = True
if do_these:
  for isubject,subject in enumerate(subjects):
    out_dir = os.path.join(out_path, subject)
    os.system('mkdir ' + out_dir)

    #-------------------------------------------------------------------------
    # Files:
    #-------------------------------------------------------------------------
    brain_file = subject_brains[isubject]
    brain_to_dti = os.path.join(out_dir, 'brain_to_dti.nii.gz')
    dti_brain = os.path.join(dticoreg_dir, subject, 'eroded_DTI_Tensor_FA.nii.gz')
    xfm = os.path.join(out_dir, 'brain_to_dti')
    labels_to_dti = os.path.join(out_dir, 'labels_to_dti.nii.gz')
    if hippo:
        labels_file = os.path.join(reg_dir, subject + '_labels1.nii.gz')
    else:
        labels_file = os.path.join(mindboggle_dir, 'labels_volume',
            '_subject_'+subject, 'labels.DKT31.DKTatlas.native.nii.gz')
        #labels_file = os.path.join(os.environ['SUBJECTS_DIR'],
        #    subject, 'mri', 'aparc+aseg.nii.gz')
    dti_brain_mask = os.path.join(out_dir, 'dti_brain_mask.nii.gz')
    if run_probtrackx:
        bedpostx_dir_subj = os.path.join(bedpostx_dir, subject + '_bedpostxoutput')
    if run_camino_prep or run_camino_track:
        tensor_file = os.path.join(dti_dir, subject, 's006a001.nii.gz')
        bval = os.path.join(dti_dir, 's006a001.bval')
        bvec = os.path.join(dti_dir, 's006a001.bvec')

    #-------------------------------------------------------------------------
    # Transform to DTI space:
    #-------------------------------------------------------------------------
    if transform_labels:
        # Register T1 to DTI:
        args = ['ANTS 3 -m CC[' + dti_brain +','+ brain_file +',1,2]',
                '-o', xfm,
                '-r Gauss[2,0] -t SyN[0.5] -i 30x11x1 --use-Histogram-Matching',
                '--number-of-affine-iterations 10000x10000x10000x10000x10000']
        #args = ['flirt -in', brain_file, '-ref', dti_brain,
        #        '-o', brain_to_dti,'-omat', xfm] 
        cmd = ' '.join(args)
        print(cmd); os.system(cmd)
        args = ['c3d', dti_brain, '-binarize -o', dti_brain_mask]
        cmd = ' '.join(args)
        print(cmd); os.system(cmd)

        # Transform labels to DTI:
        args = ['WarpImageMultiTransform 3', labels_file, labels_to_dti,
                '-R', dti_brain, xfm + 'Warp.nii.gz', xfm + 'Affine.txt --use-NN'] 
        #args = ['flirt -in', labels_file, '-applyxfm -init', xfm,
        #        '-ref', dti_brain, '-interp nearestneighbour -o', labels_to_dti] 
        cmd = ' '.join(args)
        print(cmd); os.system(cmd)

    #-------------------------------------------------------------------------
    # Loop through labels:
    #-------------------------------------------------------------------------
    for label in labels:
        out_dir_label = os.path.join(out_dir, 'label'+str(label))
        os.system('mkdir ' + out_dir_label)
        label_file = os.path.join(out_dir_label, 'label.nii.gz')
        label_to_dti = os.path.join(out_dir_label, 'label_to_dti.nii.gz')
        if transform_probtrackx:
            probtrackx_file = os.path.join(out_dir_label, 'fdt_paths.nii.gz')
            probtrackx_to_orig = os.path.join(out_dir_label, 'fdt_paths_to_orig.nii.gz')

        #---------------------------------------------------------------------
        # Extract / transform label regions:
        #---------------------------------------------------------------------
        if extract_label:
            if hippo:
                if label == 1:
                    args = ['c3d', labels_file,
                            '-threshold', '1 6 1 0 -o', label_file]
                elif label == 2:
                    args = ['c3d', labels_file,
                            '-threshold', '101 106 1 0 -o', label_file]
            else:
                args = ['c3d', labels_file,
                        '-threshold', str(label), str(label), '1 0 -o', label_file]
            cmd = ' '.join(args)
            print(cmd); os.system(cmd)

            if transform_label == True:
                args = ['WarpImageMultiTransform 3', label_file, label_to_dti,
                        '-R', dti_brain, xfm + 'Warp.nii.gz', xfm + 'Affine.txt --use-NN']
                #args = ['flirt -in', label_file, '-applyxfm -init', xfm,
                #        '-ref', dti_brain, '-interp nearestneighbour -o', label_to_dti]
                cmd = ' '.join(args)
                print(cmd); os.system(cmd)

        elif extract_dti_label:
            args = ['c3d', labels_to_dti,
                    '-threshold', str(label), str(label), '1 0 -o', label_to_dti]
            cmd = ' '.join(args)
            print(cmd); os.system(cmd)

        #---------------------------------------------------------------------
        # FSL tractography (assumes bedpostx already run):
        #---------------------------------------------------------------------
        if run_probtrackx:
            args = ['probtrackx --mode=seedmask',
                    '-x', label_to_dti,
                    '-l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd',
                    '-s', os.path.join(bedpostx_dir_subj, 'merged'),
                    '-m', dti_brain_mask, '--dir=' + out_dir_label]
            cmd = ' '.join(args)
            print(cmd)
            os.system(cmd)

        if transform_probtrackx:
            args = ['WarpImageMultiTransform 3', probtrackx_file,
                    probtrackx_to_orig,
                    '-R', brain_file, 
                    '-i', xfm + 'Affine.txt', 
                    xfm + 'InverseWarp.nii.gz --use-NN']
            cmd = ' '.join(args)
            print(cmd); os.system(cmd)

        #---------------------------------------------------------------------
        # Camino DTI processing:
        #---------------------------------------------------------------------
        if run_camino_prep:
            from run_camino import run_camino
            run_camino(out_dir, tensor_file, bval, bvec, dti_brain_mask,
               roi_file=label_to_dti, annotation='_'+str(label), track_type=False)

        #---------------------------------------------------------------------
        # Camino tractography:
        #---------------------------------------------------------------------
        if run_camino_track:
            from run_camino import run_camino
            run_camino(out_dir, tensor_file, bval, bvec, dti_brain_mask,
               roi_file=label_to_dti, annotation='_'+str(label), track_type='fact')
