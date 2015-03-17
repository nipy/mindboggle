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

#include "VtkFileEditor.h"

VtkFileEditor::VtkFileEditor(char *fileName)
{
    m_inFile.open (fileName, ios::in);
    m_outFile.open (fileName, ios::out | ios::app);
}


VtkFileEditor::~VtkFileEditor()
{
    m_inFile.close();
    m_outFile.close();
}

void VtkFileEditor::CreateField(char *fieldName, vtkDoubleArray *values)
{
    m_outFile<<endl;
    m_outFile<<"SCALARS "<<fieldName<<" double"<<endl;
    m_outFile<<"LOOKUP_TABLE default"<<endl;

    for(int i = 0; i<values->GetNumberOfTuples() ; i++)
    {
        m_outFile<<values->GetValue(i)<<endl;
    }

}
