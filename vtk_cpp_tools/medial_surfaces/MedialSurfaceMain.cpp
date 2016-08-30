#include "MedialSurfaceComputer.h"

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

int main(int argc, char** argv)
{
    time_t start= time(NULL);

    vtkPolyData* travelData;
    vtkPolyData* fundiData;

    std::string extension1 =
      vtksys::SystemTools::GetFilenameLastExtension(argv[1]);
    std::string extension2 =
      vtksys::SystemTools::GetFilenameLastExtension(argv[2]);

      if (extension1 == ".vtk")
        travelData = ReadObject<vtkPolyDataReader> (argv[1]);
      else if (extension1 == ".obj")
        travelData = ReadObject<vtkMNIObjectReader> (argv[1]);
      else if (extension1 == ".stl")
        travelData = ReadObject<vtkSTLReader> (argv[1]);

      if (extension2 == ".vtk")
        fundiData = ReadObject<vtkPolyDataReader> (argv[2]);
      else if (extension2 == ".obj")
        fundiData = ReadObject<vtkMNIObjectReader> (argv[2]);
      else if (extension2 == ".stl")
        fundiData = ReadObject<vtkSTLReader> (argv[2]);


    MedialSurfaceComputer* msc = new MedialSurfaceComputer( travelData, fundiData );

    msc->WriteIntoFile(argv[3]);

    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}
