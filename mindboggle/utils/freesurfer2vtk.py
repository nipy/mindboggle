#!/usr/bin/python

"""
Surface conversion from FreeSurfer to VTK polydata format.


Authors:  Forrest Sheng Bao http://fsbao.net
Version:  0.2, last update on 2012-06-19

(c) 2012  Mindbogglers (www.mindboggle.info), under Apache License Version 2.0

"""

def read_surface(filename):
    import struct, os
    f = open(filename, "rb")
    f.seek(3)  # skip the first 3 Bytes "Magic" number

    s = f.read(50)   # the second field is string of creation information of variable length
    End2 = s.find('\n\n',0)  # end of the second field is a '\n\n'

    f.seek(3+End2+2)  # jump to immediate Byte after the creating information

    s = f.read(8)
    vertex_count, face_count = struct.unpack(">ii", s)
    # print("This hemisphere has", vertex_count, "vertices and", face_count, "Faces")

    vertices, faces = [], []

    for i in xrange(0, vertex_count):
        s = f.read(8)
        R, A = struct.unpack(">ff", s)
        f.seek(-4, os.SEEK_CUR)
        s = f.read(8)
        A, S = struct.unpack(">ff", s)
        vertices.append([R,A,S]) # R, A, S are the coordinates of vertices
        #print i

    for i in xrange(0, face_count):
        s = f.read(8)
        V0, V1 = struct.unpack(">ii", s)
        f.seek(-4, os.SEEK_CUR)
        s = f.read(8)
        V1, V2 = struct.unpack(">ii", s)
        faces.append([V0, V1, V2])
        #print i, V0, V1, V2

    return vertices, faces

def write_header(Fp, Title='', Header='# vtk DataFile Version 2.0', FileType='ASCII', DataType='POLYDATA'):
    '''Write the all non-data information for a VTK-format file
    
    This part matches three things in VTK 4.2 File Formats doc
    
    Part 1: Header
    Part 2: Title (256 characters maximum, terminated with newline \n character)
    Part 3: Data type, either ASCII or BINARY
    Part 4: Geometry/topology. Type is one of:
        STRUCTURED_POINTS
        STRUCTURED_GRID
        UNSTRUCTURED_GRID
        POLYDATA
        RECTILINEAR_GRID
        FIELD

    '''
    
    Fp.write(Header)
    Fp.write("\n")
    Fp.write(Title)
    Fp.write("\n")
    Fp.write(FileType)
    Fp.write("\n")
    Fp.write("DATASET ")
    Fp.write(DataType)
    Fp.write("\n")
            
def write_points(Fp, PointList, Type="float"):
    """Print coordinates of points, the POINTS section in DATASET POLYDATA section 
    """
    Fp.write("POINTS " + str(len(PointList)) + " " + Type + "\n")
    for i in xrange(0, len(PointList)):
        [R, A, S] = PointList[i]
        Fp.write(str(R) + " " + str(A) + " " + str(S) + "\n")
    
def write_faces(f, face_list, vertices_per_face=3):
    """Print vertices forming triangular meshes, the POLYGONS section in DATASET POLYDATA section 
    """
    f.write("POLYGONS " + str(len(face_list)) + " " + str( (vertices_per_face + 1) * len(face_list)  )  + '\n' )
    for i in xrange(0, len(face_list)):
        [V0, V1, V2] = face_list[i]
        f.write( str(vertices_per_face) + " " + str(V0) + " " + str(V1) + " " + str(V2) + "\n")

def freesurfer2vtk(in_file):
    """
    Convert FreeSurfer surface file to vtk format
    """

    from os import path, getcwd, error
    from freesurfer2vtk import read_surface, write_header, write_points, write_faces

    # Check type:
    if type(in_file) == str:
        pass
    elif type(in_file) == list:
        in_file = in_file[0]
    else:
        error("Check format of " + in_file)

    vertices, faces = read_surface(in_file)

    out_file = path.join(getcwd(), in_file + '.vtk')
    f = open(out_file, 'w')

    write_header(f)
    write_points(f, vertices)
    write_faces(f, faces)

    f.close()

    return out_file

