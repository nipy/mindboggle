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

#include "FsSurfaceReader.h"



#include <vtkPoints.h>

#include <vtkCellArray.h>
#include <vtkIdList.h>

#include <vtkPolyDataWriter.h>

#include <vtkPolyDataNormals.h>

//using namespace std;

FsSurfaceReader::FsSurfaceReader(char *fileName)
{
    ifstream myfile;
    myfile.open (fileName, ios::in | ios::binary);

    vtkPolyData* pd = vtkPolyData::New();

    int fileIdType = readInt3(myfile);

    std::string creator = readTo2eol(myfile);

    unsigned char tmp[4];
    int x;

    m_mesh = vtkPolyData::New();

    myfile.read((char*)tmp, 4);
    x = (((int)(tmp[0])) << 24) |
        (((int)(tmp[1])) << 16) |
        (((int)(tmp[2])) << 8) |
            (int)tmp[3];
    int vertexCount = x;
    cout<<"vertex count "<<x<<endl;

    myfile.read((char*)tmp, 4);
    x = (((int)(tmp[0])) << 24) |
        (((int)(tmp[1])) << 16) |
        (((int)(tmp[2])) << 8) |
            (int)tmp[3];
    int faceCount = x;
    cout<<"face count "<<x<<endl;


    vtkPoints* points = vtkPoints::New();
    float pos[3];
    for(int i = 0 ; i<vertexCount ; i++){
        for(int j = 0; j<3 ; j++){
            myfile.read((char*)tmp, 4);
            x = (((int)(tmp[0])) << 24) |
                (((int)(tmp[1])) << 16) |
                (((int)(tmp[2])) << 8) |
                    (int)tmp[3];

            float *f = (float*)&x;
            float y = *f;
            pos[j] = y;
        }
        points->InsertNextPoint(pos);
    }

    cout<<"points Inserted"<<endl;


    vtkIndent indent;
    points->PrintSelf(cout, indent);

    pd->SetPoints(points);

    vtkCellArray* cells = vtkCellArray::New();
    for(int i = 0 ; i<faceCount ; i++){
        vtkIdList* pointList = vtkIdList::New();
        for(int j = 0; j<3 ; j++){
            myfile.read((char*)tmp, 4);
            x = (((int)(tmp[0])) << 24) |
                (((int)(tmp[1])) << 16) |
                (((int)(tmp[2])) << 8) |
                    (int)tmp[3];

            pointList->InsertNextId(x);
        }

        cells->InsertNextCell(pointList);

    }

    cout<<"faces Inserted"<<endl;

//    cells->PrintSelf(cout, indent);


    pd->SetPolys(cells);
    pd->Update();

    cout<<endl;

    pd->PrintSelf(cout,indent);

    vtkPolyDataNormals *pdn = vtkPolyDataNormals::New();
    pdn->SetInputData(pd);
    pdn->SetFeatureAngle(90);
    pdn->SplittingOff();
    pdn->Update();

    cout<<endl;

    pdn->GetOutput()->PrintSelf(cout,indent);

    m_mesh->DeepCopy(pdn->GetOutput());

//    vtkPolyDataWriter* pdw = vtkPolyDataWriter::New();
//    pdw->SetInputData(m_mesh);
//    pdw->SetFileName("test.vtk");
//    pdw->Write();
//    pdw->Update();
//    pdw->Delete();

    cells->Delete();
    pdn->Delete();
    points->Delete();

}

FsSurfaceReader::~FsSurfaceReader()
{
    m_mesh->Delete();
}


vtkPolyData* FsSurfaceReader::GetVTKData()
{
    return m_mesh;
}

int FsSurfaceReader::readInt3(ifstream &myfile)
{
    unsigned char tmp[4];
    int x;
    myfile.read((char*)tmp, 3);
    x = (((int)(tmp[0])) << 16) |
        (((int)(tmp[1])) << 8) |
            (int)tmp[2];
//    cout<<"File Type Id "<<x<<endl;
    return x;

}

std::string FsSurfaceReader::readTo2eol(ifstream &myfile)
{
    std::string str;

    char prevChar;
    char currentChar;
    myfile.read(&prevChar, 1);
    myfile.read(&currentChar, 1);

    while (!(prevChar == '\n' && currentChar == '\n')) {
        str += prevChar;
        prevChar = currentChar;
        myfile.read(&currentChar, 1);
    }

//    cout<<"user "<< str <<endl;
    return str;
}
