#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Modeler.h"

void help();

// masking > 

int main(int argc, char* argv[])
{
    	CModeler* surface=new CModeler;
	
	char* objfile;
	char* valuefile;
	char* pitfile;
	// In input pit file, label number should be smaller than 10000	 
	char* savefile;	// sulcal parcellated file.
	int option1; 	// from min or from max
	float final_value;	// stop value of feature for parcellation
	int i,j;
	
	if(argc<2){
		help();
		return(0);
	}
	
	if (strcmp(argv[1],"help")==0 || (strcmp(argv[1],"-from_max")!=0 && strcmp(argv[1],"-from_min")!=0)){
		help();
		return (0);
	}

	if (strcmp(argv[1],"-from_max")==0)
		option1=0;
	else	// -from_min
		option1=1;
	
	objfile=argv[3];
	if (strcmp(argv[2],"-fs")==0){
		if(!surface->LoadAscFile(objfile)){
			delete surface;
			return(0);
		}
	}
	else if(strcmp(argv[2],"-mni")==0){
		if(!surface->LoadObjFile(objfile)){
			delete surface;
			return(0);
		}
	}
	
	valuefile=argv[4];
	pitfile=argv[5];
	final_value=atof(argv[6]);
	savefile=argv[7];

	if(!surface->LoadDepthFile(valuefile)){
		delete surface;
		return(0);
	}

	if(!surface->LoadValueFile(pitfile)){
		delete surface;
		return(0);
	}

	surface->FindNeighborInfo(surface->m_TriangleInfo, surface->m_VertexInfo, surface->m_VertexNum, surface->m_TriangleNum);
	
	surface->m_label=new float [surface->m_VertexNum];		
	memset( surface->m_label, 0, surface->m_VertexNum*sizeof(float));
	
	surface->SulcalParcel_Pit(surface->m_VertexInfo, surface->m_VertexNum, surface->m_depth, surface->m_value, final_value, option1, surface->m_label);
	surface->SaveMniTexFile(savefile, surface->m_label,surface->m_VertexNum);

	delete surface;
	printf("Finished\n");
	return(0);
}

void help()
{
	printf("\nDeveloped by Kiho Im (kiho.sky@gmail.com)\n\n");
	printf("Usage : mesh_sulcalparcel_usingpit -from_max/-from_min -fs/-mni a.asc/a.obj field.map pit.txt final_value save_region.txt\n");
	printf("watershed progress from_max: max value to final_value (decrease), -from_min: min value to final_value (increase)\n");
	printf("final_value: do watershed until final_value is met\n");
}
