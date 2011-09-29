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
#libfundifc.getFundi(CurvFile, SurfFile, ToVTK=ToVTK, SurfFile2=SurfFile2)

# extract fundus curves from pits
#libfundipit.getFundi(CurvFile, SurfFile, ToVTK=ToVTK, SurfFile2=SurfFile2)