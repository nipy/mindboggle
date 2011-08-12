import libfundivtx, libfundifc, libbasin, fileio, libfundipit, libvtk # this line imports my own library
import sys, getopt # this line imports pulic libraries

try:
    opts, args = getopt.getopt(sys.argv[1:],'S:T')
except getopt.GetoptError, err:
    print str(err) # will print something like "option -a not recognized"
    print "\n\t Usage: python extract.py [-T] [-S Second Surface File] CurvFile SurfFile"
    sys.exit(2)

if len(args)<2:
    print "\n\t Please give a CurvFile and a SurfFile"
    print "\t Usage: python extract.py [-T] [-S Second Surface File] CurvFile SurfFile"
    sys.exit()

[CurvFile, SurfFile] = args

#print args
#print opts

print "\tRead in FreeSurfer curvature-format per-vertex file:", CurvFile
print "\tRead in FreeSurfer surfer file:", SurfFile

ToVTK = False
SurfFile2 = ''

#print opts

for o,p in opts:
    if o in ['-T']:
        ToVTK = True
        print "\tResults will be mapped onto a surface and saved in VTK format"
        libvtk.surf2VTK(SurfFile)
    elif o in ['-S']:
        SurfFile2 = p
        print "\tExtraction result be mapped onto the second surface file:", SurfFile2
        libvtk.surf2VTK(SurfFile2)


# extract basin and pits
libbasin.getBasin(CurvFile, SurfFile, ToVTK, SurfFile2 = SurfFile2)

# extract fundus curves from vertexes directly as Gang Li's
libfundivtx.getFundi(CurvFile, SurfFile, ToVTK=ToVTK, SurfFile2=SurfFile2)

# extract fundus curves from fundus face strips
libfundifc.getFundi(CurvFile, SurfFile, ToVTK=ToVTK, SurfFile2=SurfFile2)

# extract fundus curves from pits
libfundipit.getFundi(CurvFile, SurfFile, ToVTK=ToVTK, SurfFile2=SurfFile2)


'''
import libfundivtx, libfundifc, libbasin, fileio, libfundipit# this line imports my own library
import sys, getopt # this line imports pulic libraries

def printUsage():
    print "\n\t Please give a CurvFile and a SurfFile"
    print "\t Usage: python extract.py [-T] [-F V/S/P] [-S Second Surface File] CurvFile SurfFile"
    print "-T: Save results into VTK format"
    print "-F: the source of fundi. Fundi can be extracted from vertexes directly(V, default), or skeleton of fundus face strips (S), or pits (P)"
    print "-S: secondary surface file where results will be mapped onto and saved in VTK format. For example, if primiary surface is the pial surface, the secondary surface can be inflated surface"

try:
    opts, args = getopt.getopt(sys.argv[1:],'S:F:T')
except getopt.GetoptError, err:
    print str(err) # will print something like "option -a not recognized"
    printUsage()
    sys.exit(2)

#print opts

if len(args)<2:
    print "\n\t Please give a CurvFile and a SurfFile"
    printUsage()
    sys.exit()

[CurvFile, SurfFile] = args

#print args
#print opts

print "\tRead in FreeSurfer curvature-format per-vertex file:", CurvFile
print "\tRead in FreeSurfer surfer file:", SurfFile

ToVTK = False
SurfFile2 = ''
FundusFrom = 'P'

for o,p in opts:
    if o in ['-T']:
        ToVTK = True
        print "\tResults will be mapped onto a surface and saved in VTK format"
    elif o in ['-S']:
        SurfFile2 = p
        print "\tExtraction result be mapped onto the second surface file:", SurfFile2
    elif o in ['-F']:
        if p == 'V':
            print "\t Fundi will be extracted from vertexes directly"
            FundusFrom = p
        elif p == 'S':
            print "\t Fundi will be extracted from skeleton of fundus face strips"
            FundusFrom = p
        elif p == 'P':
            print "\t Fundi will be extracted from Pits"
            FundusFrom = p

# extrac basin
#libbasin.getBasin(CurvFile, SurfFile, ToVTK, SurfFile2 = SurfFile2)

if FundusFrom == 'V':
# extract fundi directly from vertex cloud/curve segments
    libfundivtx.getFundi(CurvFile, SurfFile, ToVTK=ToVTK, SurfFile2=SurfFile2)
elif FundusFrom == 'S':
# extract fundi from face strips
    libfundifc.getFundi(CurvFile, SurfFile, ToVTK=ToVTK, SurfFile2=SurfFile2)
elif FundusFrom == 'P':
# extract fundi from face strips
    libfundipit.getFundi(CurvFile, SurfFile, ToVTK=ToVTK, SurfFile2=SurfFile2)
'''
