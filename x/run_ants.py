"""
Run antsCorticalThickness.sh on Mindboggle-101 brains
"""
import os

run_all = True

if run_all:

    names = ['OASIS-TRT-20', 'MMRR-21', 'NKI-RS-22', 'NKI-TRT-20',
             'Afterthought', 'Colin27', 'Twins-2', 'MMRR-3T7T-2', 'HLN-12']
    numbers = [20,21,22,20,1,1,2,2,12]

    i1 = 0
    names = [names[i1]]
    numbers = [numbers[i1]]

    path1 = '/homedir/Data/Brains/Mindboggle101/subjects/'
    end1a = '/mri/orig/001.mgz'
    end1b = '/mri/orig/001.nii.gz'
    path2 = '/data/Brains/Atropos_templates/OASIS-30_Atropos_template/'
    end2a = 'T_template0.nii.gz'
    end2b = 'T_template0_BrainCerebellum.nii.gz'
    end2c = 'T_template0_BrainCerebellumProbabilityMask.nii.gz'
    end2d = 'T_template0_BrainCerebellumExtractionMask.nii.gz'
    end2e = 'Priors2/priors%d.nii.gz'
    convert = False

    for i,name in enumerate(names):
        number = numbers[i]
        for n in range(1,number+1):
            if convert:
                s = 'mri_convert {0}{1}-{2}{3} {0}{1}-{2}{4} ' \
                    .format(path1, name, n, end1a, end1b)
                print(s)
                os.system(s)

            prefix = 'antsCorticalThickness/{0}-{1}/ants'.format(name, n)

            s = 'antsCorticalThickness.sh -d 3 -n 3 -w 0.25 ' \
                '-a {0}{1}-{2}{3} ' \
                '-o {4} ' \
                '-e {5}/{6} ' \
                '-t {5}/{7} ' \
                '-m {5}/{8} ' \
                '-f {5}/{9} ' \
                '-p {5}/{10} ' \
                .format(path1, name, n, end1b, prefix, path2, end2a, end2b,
                        end2c, end2d, end2e)
            print(s)
            os.system(s)
