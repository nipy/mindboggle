#include "../GradientComputer.h"

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkMNIObjectReader.h>
#include <vtkSTLReader.h>
#include <vtksys/SystemTools.hxx>


#include <iostream>
#include <fstream>
#include <stdio.h>

template<class TReader> vtkPolyData *ReadObject(const char*fileName)
{
  vtkSmartPointer<TReader> reader =
    vtkSmartPointer<TReader>::New();
  reader->SetFileName(fileName);

  reader->Update();
  reader->GetOutput()->Register(reader);
  return reader->GetOutput();
}

void print_help()
{
    printf(
    "Usage: GradientComputer InputVTKMesh\n"
        );
}

int main(int argc, char** argv)
{
    time_t start= time(NULL);
    vtkPolyData *dataSet;
    std::string extension = ".vtk";
      //vtksys::SystemTools::GetFilenameLastExtension(argv[1]);
    if (extension == ".vtk")
      dataSet = ReadObject<vtkPolyDataReader> (argv[1]);
    else if (extension == ".obj")
      dataSet = ReadObject<vtkMNIObjectReader> (argv[1]);
    else if (extension == ".stl")
      dataSet = ReadObject<vtkSTLReader> (argv[1]);

    GradientComputer* area = new GradientComputer(dataSet);
    area->ComputeGradient();
    area->WriteIntoFile(argv[2]);

    cout<<"Elapsed time (area): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}
