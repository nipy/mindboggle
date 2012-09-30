#!/usr/bin/python
"""
Copy structures and maps from one surface to another, vertex by vertex    

Authors:
    - Forrest Sheng Bao  (forrest.bao@gmail.com)  http://fsbao.net

Copyright 2012,  Mindboggle team (http://mindboggle.info), Apache v2.0 License

To make this a stand alone script, some functions are copied to here from mindboggle.utils.io_vtk

Usage: 
     python maps_to_second_surface.py file_containing_structures_and_maps file_containing_2nd_surface output_file
     
Examples:
     python maps_to_second_surface.py lh.gaussian_curv.pial.vtk lh.inflated.vtk lh.gaussian_curv.inflated.vtk

Dependencies:
    python-vtk: vtk's official Python binding
"""

import vtk
import sys

def load_maps(Filename):
    """Load all maps from the first VTK files 
    
    parameters
    --------------
    Filename: string
        the path to the first VTK file
        
    returns
    ---------
    MapsInFile: list of lists of floats/integers
        Each element (size: 1 by #vertexes) is a list represent a map on the surface
    
    Note
    ------
    
        1. This function differs from load_scalar in io_vtk which loads only one  scalar map
        2. For structures, we only support copying lines, vertexes (indexes of points) and triangular faces from one surface to another.  
        3. We assume that all vertexes are written in one line in VERTICES segment. 
           If your VERTICES structure is different, please modify the code accordingly  
    
    """
    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(Filename)
    Reader.ReadAllScalarsOn()  # Activate the reading of all scalars
    Reader.Update()
    
    Data = Reader.GetOutput()
    PointData = Data.GetPointData()
    
    if Data.GetNumberOfPolys() > 0:
        Faces = [[Data.GetPolys().GetData().GetValue(j) for j in xrange(i*4 + 1, i*4 + 4)]
                 for i in xrange(Data.GetPolys().GetNumberOfCells())]
    else:
        Faces = []
    
    if Data.GetNumberOfLines() > 0:
        Lines  = [[Data.GetLines().GetData().GetValue(j) for j in xrange(i*3+1, i*3+3) ]
                  for i in xrange(Data.GetNumberOfLines())]
    else:
        Lines = []

    if Data.GetNumberOfVerts() > 0:
        Vertices = [Data.GetVerts().GetData().GetValue(i) for i in xrange(1, Data.GetVerts().GetSize() )]
        # The reason the reading starts from 1 is because we need to avoid the 
    else:
        Vertices = []
    
    MapsInFile = []
    MapNamesInFile = []
    
    if Reader.GetNumberOfScalarsInFile() > 0:
        for Scalar_Index in range(Reader.GetNumberOfScalarsInFile()):
            Scalar_Name = Reader.GetScalarsNameInFile(Scalar_Index)
            print("Loading {0} (named \"{1}\") of {2} scalars in file {3}...".
              format(Scalar_Index + 1,
                     Reader.GetScalarsNameInFile(Scalar_Index), Reader.GetNumberOfScalarsInFile(), Filename))
            ScalarArray = PointData.GetArray(Scalar_Name)
            if ScalarArray:
                Scalar = [ScalarArray.GetValue(i) for i in xrange(ScalarArray.GetSize())]
            else:
                print "An empty scalar map was read. Please check the integrity of the source VTK"
                sys.exit()
            MapsInFile.append(Scalar)
            MapNamesInFile.append(Scalar_Name) 
        
    return Faces, Lines, Vertices, MapsInFile, MapNamesInFile

def load_2nd_surface(Filename):
    """
    Load points of the 2nd surface.

    Parameters
    ----------
    filename : string
        The path/filename of a VTK format file.

    Returns
    -------
    points : list of lists of floats (see return_arrays)
        Each element is a list of 3-D coordinates of a vertex on a surface mesh

    """
    Reader = vtk.vtkDataSetReader()
    Reader.SetFileName(Filename)
    Reader.ReadAllScalarsOn()  # Activate the reading of all scalars
    Reader.Update()

    Data = Reader.GetOutput()
    Points = [list(Data.GetPoint(point_id))
              for point_id in xrange(Data.GetNumberOfPoints())]

    return Points

def write_2nd_surface_with_maps(vtk_file, Points, Faces, Lines, Indices, LUTs, LUT_names, Title="Created by MindBoggle"):
    """Write the second surface along with maps from the 1st file into the target file 
   
    Scalar definition includes specification of a lookup table.
    The definition of a lookup table is optional. If not specified,
    the default VTK table will be used (and tableName should be "default").

    SCALARS dataName dataType numComp
    LOOKUP_TABLE tableName

    Parameters
    ----------
    vtk_file : string
        path of the output VTK file
    Points :  list of 3-tuples of floats
        each element has 3 numbers representing the coordinates of the points
    Indices : list of integers
        indices of vertices. they are called VERTICES in VTK format though
    Faces : list of 3-tuples of integers
        indices to the three vertices of a face on the mesh
    Lines : list of 2-tuples of integers
        Each element contains two indices  of the two mesh vertexes forming a line on the mesh
    LUTs : list of lists of floats
        each list contains values assigned to the vertices
    LUT_names : list of strings
        each element is the name of a LUT
    Title : string
        Title for the target VTK file

    """    
    
    def write_vtk_header(Fp, Header='# vtk DataFile Version 2.0',
                         Title='Generated by Mindboggle (www.mindboggle.info)',
                         fileType='ASCII', dataType='POLYDATA'):
        """
        Write header information for a VTK-format file::
    
            # vtk DataFile Version 2.0
            Generated by Mindboggle (www.mindboggle.info)
            ASCII
            DATASET POLYDATA
    
        This part matches three things in the VTK 4.2 File Formats doc:
          - Part 1: Header
          - Part 2: Title (256 characters maximum, terminated with newline character)
          - Part 3: Data type, either ASCII or BINARY
          - Part 4: Geometry/topology. dataType is one of:
              - STRUCTURED_POINTS
              - STRUCTURED_GRID
              - UNSTRUCTURED_GRID
              - POLYDATA
              - RECTILINEAR_GRID
              - FIELD
    
        """
    
        Fp.write('{0}\n{1}\n{2}\nDATASET {3}\n'.format(Header, Title, fileType, dataType))
    
    def write_vtk_points(Fp, points, dataType="float"):
        """
        Write coordinates of points, the POINTS section in DATASET POLYDATA section::
    
            POINTS 150991 float
            -7.62268877029 -81.2403945923 -1.44539153576
            ...
    
        Indices are 0-offset. Thus the first point is point id 0::
    
            POINTS n dataType
            p0x p0y p0z
            ...
            p(n-1)x p(n-1)y p(n-1)z
    
        """
        import numpy as np
    
        Fp.write('POINTS {0} {1}\n'.format(len(points), dataType))
    
        n = np.shape(points)[1]
        for point in points:
            if n == 3:
                [R, A, S] = point
                Fp.write('{0} {1} {2}\n'.format(R, A, S))
            elif n == 2:
                [R, A] = point
                Fp.write('{0} {1}\n'.format(R, A))
            else:
                print('ERROR: Unrecognized number of coordinates per point')
    
    def write_vtk_faces(Fp, faces):
        """
        Write indices to vertices forming triangular meshes,
        the POLYGONS section in DATASET POLYDATA section::
    
            POLYGONS 301978 1207912
            3 0 1 4
            ...
    
        """
        import numpy as np
    
        n = np.shape(faces)[1]
        if n == 3:
            face_name = 'POLYGONS '
            Fp.write('{0} {1} {2}\n'.format(face_name, len(faces),
                     len(faces) * (n + 1)))
        elif n == 2:
            face_name = 'LINES '
            Fp.write('{0} {1} {2}\n'.format(face_name, len(faces),
                     len(faces) * (n + 1)))
        else:
            print('ERROR: Unrecognized number of vertices per face')
    
        for face in faces:
            if n == 3:
                [V0, V1, V2] = face
                Fp.write('{0} {1} {2} {3}\n'.format(n, V0, V1, V2))
            elif n == 2:
                [V0, V1] = face
                Fp.write('{0} {1} {2}\n'.format(n, V0, V1))

    def write_vtk_lines(Fp, lines):
        """
        Save connected line segments to a VTK file.
    
        Parameters
        ----------
        Fp: pointer to a file
            pointer to the file to write lines
        lines : list of 2-tuples of integers
            Each element is an edge on the mesh, consisting of 2 integers
            representing the 2 vertices of the edge
        """

        write_vtk_faces(Fp, lines)
    
    def write_vtk_vertices(Fp, indices):
        """
        Write indices to vertices, the VERTICES section
        in the DATASET POLYDATA section::
    
            VERTICES 140200 420600
            3 130239 2779 10523
            ...
    
        Indices are 0-offset. Thus the first point is point id 0::
    
            VERTICES n size
            numPoints0 i0 j0 k0
            ...
            numPoints_[n-1] i_[n-1] j_[n-1] k_[n-1]
    
        """
    
        Fp.write('VERTICES {0} {1}\n{2} '.format(
                 len(indices), len(indices) + 1, len(indices)))
        [Fp.write('{0} '.format(i)) for i in indices]
        Fp.write('\n')
    
    def write_vtk_LUT(Fp, LUT, LUTName, at_LUT_begin=True):
        """
        Write per-VERTEX values as a scalar LUT into a VTK file::
    
            POINT_DATA 150991
            SCALARS Max_(majority_labels) int 1
            LOOKUP_TABLE default
            11 11 11 11 11 11 11 11 11 11 ...
    
        Parameters
        ----------
        LUT :  list of floats
        at_LUT_begin : [Boolean] True if the first vertex LUT in a VTK file
    
        """
    
        if at_LUT_begin:
            Fp.write('POINT_DATA {0}\n'.format(len(LUT)))
        Fp.write('SCALARS {0} float\n'.format(LUTName))
        Fp.write('LOOKUP_TABLE {0}\n'.format(LUTName))
        for Value in LUT:
            Fp.write('{0}\n'.format(Value))
        Fp.write('\n')
    
    # End of nested functions
    import os

    vtk_file = os.path.join(os.getcwd(), vtk_file)

    Fp = open(vtk_file,'w')
    write_vtk_header(Fp, Title=Title)
    write_vtk_points(Fp, Points)
    if Faces != []:
        write_vtk_faces(Fp, Faces)
    if Lines !=[]:
        write_vtk_lines(Fp, Lines)
        print "Lines written"
    if Indices != []:
        write_vtk_vertices(Fp, Indices)
        print "Pits written"
    if len(LUTs) > 0:
        # Make sure that LUTs is a list of lists
        if type(LUTs[0]) != list:
            LUTs = [LUTs]
        for i, LUT in enumerate(LUTs):
            if i == 0:
                if len(LUT_names) == 0:
                    LUT_name = 'scalars'
                else:
                    LUT_name = LUT_names[i]
                write_vtk_LUT(Fp, LUT, LUT_name)
            else:
                if len(LUT_names) < i + 1:
                    LUT_name = 'scalars'
                else:
                    LUT_name  = LUT_names[i]
                write_vtk_LUT(Fp, LUT, LUT_name, at_LUT_begin=False)
    Fp.close()

    return vtk_file

if __name__ == "__main__":
    Faces, Lines, Vertices, Maps, MapNamesInFile = load_maps(sys.argv[1])
    Points = load_2nd_surface(sys.argv[2])
    if len(sys.argv) > 4:
        write_2nd_surface_with_maps(sys.argv[3], Points, Faces, Lines, Vertices, Maps, MapNamesInFile, Title=sys.argv[4])
    else:
        write_2nd_surface_with_maps(sys.argv[3], Points, Faces, Lines, Vertices, Maps, MapNamesInFile)