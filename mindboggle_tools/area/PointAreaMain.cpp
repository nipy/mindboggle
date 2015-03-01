#include "../PointAreaComputer.h"
<<<<<<< HEAD
#include "../MeshAnalyser.h"
=======
>>>>>>> e4dcb043a2f38490c0fd6c867f5864ca9eb440c4

#include <vtkPolyDataReader.h>

#include <iostream>
#include <fstream>
#include <stdio.h>

void print_help()
{
    printf(
    "Usage: PointAreaComputer InputVTKMesh\n"
        );
}

int main(int argc, char** argv)
{
    time_t start= time(NULL);

    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName(argv[1]);
    reader->Update();

    PointAreaComputer* area = new PointAreaComputer(reader->GetOutput());
    area->ComputeArea();
    area->WriteIntoFile(argv[2]);

    cout<<"Elapsed time (area): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



