#include "../Overlap.h"

#include <vtkPolyDataReader.h>

#include <iostream>
#include <fstream>
#include <stdio.h>

void print_help()
{
    printf(
    "Usage: Overlap InputVTKMesh1 InputVTKMesh2 OutputFileName\n"
        );
}

int main(int argc, char** argv)
{
    time_t start= time(NULL);

    vtkPolyDataReader* reader1 = vtkPolyDataReader::New();
    reader1->SetFileName(argv[1]);
    reader1->Update();

    vtkPolyDataReader* reader2 = vtkPolyDataReader::New();
    reader2->SetFileName(argv[2]);
    reader2->Update();

    Overlap* surfOverlap = new Overlap(reader1->GetOutput(), reader2->GetOutput());
    surfOverlap->ComputeOverlap();
    surfOverlap->WriteIntoFile(argv[3]);




    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



