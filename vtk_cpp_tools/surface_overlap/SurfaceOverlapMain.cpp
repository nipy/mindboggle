#include "../Overlap.h"

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
    "Usage: Overlap InputVTKMesh1 InputVTKMesh2 OutputFileName\n"
        );
}

int main(int argc, char** argv)
{
    time_t start= time(NULL);

    vtkPolyData* data1;
    vtkPolyData* data2;

    std::string extension1 = ".vtk";
      //vtksys::SystemTools::GetFilenameLastExtension(argv[1]);
    std::string extension2 = ".vtk";
      //vtksys::SystemTools::GetFilenameLastExtension(argv[2]);

      if (extension1 == ".vtk")
        data1 = ReadObject<vtkPolyDataReader> (argv[1]);
      else if (extension1 == ".obj")
        data1 = ReadObject<vtkMNIObjectReader> (argv[1]);
      else if (extension1 == ".stl")
        data1 = ReadObject<vtkSTLReader> (argv[1]);

      if (extension2 == ".vtk")
        data2 = ReadObject<vtkPolyDataReader> (argv[2]);
      else if (extension2 == ".obj")
        data2 = ReadObject<vtkMNIObjectReader> (argv[2]);
      else if (extension2 == ".stl")
        data2 = ReadObject<vtkSTLReader> (argv[2]);


    Overlap* surfOverlap = new Overlap(data1, data2);
    surfOverlap->ComputeOverlap();
    surfOverlap->WriteIntoFile(argv[3]);


    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}
