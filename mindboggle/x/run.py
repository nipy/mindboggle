
import os

names = ['OASIS-TRT-20', 'MMRR-21', 'NKI-RS-22', 'NKI-TRT-20']
numbers = [20,21,22,20]
prefix = 'ants'
ants_dir = '/data/Brains/Mindboggle101/antsCorticalThickness' #OASIS-TRT-20_volumes/OASIS-TRT-20-11

for i,name in enumerate(names):
    number = numbers[i]
    for n in range(1,number+1):
        s = 'mindboggle -n 8 {0}-{1} --ants_prefix {2}/{0}_volumes/{0}-{1}/ants --sulci --fundi --vertices --spectra 10 --thickness'.format(name, n, ants_dir, prefix)
        print(s)
        os.system(s)

names = ['Afterthought', 'Colin27', 'Twins-2', 'MMRR-3T7T-2', 'HLN-12']
numbers = [1,1,2,2,12]
prefix = 'ants'
ants_dir = '/data/Brains/Mindboggle101/antsCorticalThickness/Extra-18_volumes'

for i,name in enumerate(names):
    number = numbers[i]
    for n in range(1,number+1):
        s = 'mindboggle -n 8 {0}-{1} --ants_prefix {2}/{0}-{1}/ants --sulci --fundi --vertices --spectra 10 --thickness'.format(name, n, ants_dir, prefix)
        print(s)
        os.system(s)
