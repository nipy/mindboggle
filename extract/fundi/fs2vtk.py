import fileio,  libvtk # this line imports my own library
import sys  # this line imports pulic libraries

InputFile = sys.argv[1]
OutputFile = sys.argv[2]

Vrtx, Face = fileio.readSurf(InputFile)

OutFP = open(OutputFile, 'w')

libvtk.writeHeader(OutFP)
libvtk.writePoint(OutFP, Vrtx)
libvtk.writeFace(OutFP, Face)

OutFP.close()