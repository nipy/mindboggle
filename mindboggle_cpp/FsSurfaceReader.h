/* ********************************************************************
 * FsSurfaceReader
 *
 * Copyright 2011 Joachim Giard
 * Université catholique de Louvain, Belgium
 * Apache License, Version 2.0
 *
 * Class containing methods to read a freesurfer mesh file and to
 * convert it to vtk data.
 *
 * *******************************************************************/

#ifndef FS_SURFACE_READER_H_
#define FS_SURFACE_READER_H_

#include <vtkPolyData.h>

#include <iostream>
#include <fstream>

using namespace std;

class FsSurfaceReader
{
public:

        //constructor and destructor
        FsSurfaceReader(char* fileName);
        ~FsSurfaceReader();

        vtkPolyData* GetVTKData();

private:

        vtkPolyData* m_mesh;
        int readInt3(ifstream &myfile);
        std::string readTo2eol(ifstream &myfile);

};


#endif // FS_SURFACE_READER_H_
