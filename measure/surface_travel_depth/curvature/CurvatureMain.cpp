#include "../MeshAnalyser.h"

#include <iostream>
#include <fstream>
#include <stdio.h>

int main(int argc, char** argv)
{
    time_t start= time(NULL);

    MeshAnalyser* ma = new MeshAnalyser(argv[1]);
    ma->ComputePrincipalCurvatures();
    ma->WriteIntoFile(argv[2], (char*)"curv");


    cout<<"Elapsed time (meshTest): "<<time(NULL)-start<<" s"<<endl;
    return 0;
}



