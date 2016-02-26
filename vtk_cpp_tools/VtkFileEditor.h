/* ********************************************************************
 * VtkFileEditor
 *
 * Copyright 2011 Joachim Giard
 * Université catholique de Louvain, Belgium
 * Apache License, Version 2.0
 *
 * Class containing methods to edit a vtk file
 *
 * *******************************************************************/

#ifndef VTKFILEEDITOR_H
#define VTKFILEEDITOR_H

#include <vtkDoubleArray.h>

#include <iostream>
#include <fstream>

using namespace std;

class VtkFileEditor
{
public:
    //constructor and destructor
    VtkFileEditor(char* fileName);
    ~VtkFileEditor();

    void CreateField(char* fieldName, vtkDoubleArray* values);

private:
    ifstream m_inFile;
    ofstream m_outFile;

};

#endif // VTKFILEEDITOR_H
