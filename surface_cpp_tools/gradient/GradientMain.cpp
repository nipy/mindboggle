#include "../GradientComputer.h"

#include <vtkPolyDataReader.h>

#include <iostream>
#include <fstream>
#include <stdio.h>

void print_help()
{
    printf(
    "Usage: GradientComputer InputVTKMesh\n"
        );
}

int main(int argc, char** argv)
{
    time_t start= time(NULL);

    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName(argv[1]);
//  VTK6 Update: http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Removal_of_Update
//  ???
    reader->Update();

    GradientComputer* area = new GradientComputer(reader->GetOutput());
    area->ComputeGradient();
    area->WriteIntoFile(argv[2]);

    cout<<"Elapsed time (area): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



