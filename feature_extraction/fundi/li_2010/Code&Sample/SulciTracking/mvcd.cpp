#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mvcd.h"

/*****OPEN FILE*****/

FILE *myopen(char *s, char *io)
{
	FILE *fp;

	if((fp=fopen(s,io))==NULL)
    {
		perror(s);
		
		exit(1);
    }
	
	return(fp);
}
 

/*****COMPUTE MAGNITUDE OF A 2D VECTOR*****/

double Fvector2dmag(Fvector2d vector)
{
	double mag;
  
	mag=sqrt((double) (vector.x*vector.x+vector.y*vector.y));

	return mag;
}

/*****COMPUTE MAGNITUDE OF A 3D VECTOR*****/

double Fvector3dmag(Fvector3d vector)
{
	double mag;

	mag=sqrt((double) (vector.x*vector.x+vector.y*vector.y+vector.z*vector.z));

	return mag;
}

/*****MEMORY ALLOCATION*****/

/*Allocate 1d array of type char*/

char *Calloc1d(int i_size)
{
	char *array;

	array=(char *) calloc(i_size,sizeof(char ));

	return(array);
}
  

/*Allocate 2d array of type char*/

char **Calloc2d(int i_size,int j_size)
{
	char **array;
	int i;

	array=(char **) calloc(i_size,sizeof(char *));

	for(i=0;i<i_size;i++)
		array[i]=(char *) calloc(j_size,sizeof(char ));

	return(array);
}
  
/*Allocate 3d array of type char*/

char ***Calloc3d(int i_size,int j_size,int k_size)
{
	char ***array;
	int i,k;

	array=(char ***) calloc(k_size,sizeof(char **));

	for(k=0;k<k_size;k++)
		array[k]=(char **) calloc(i_size,sizeof(char *));

	for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
		array[k][i]=(char *) calloc(j_size,sizeof(char ));
	
	return(array);
}

/*Allocate 1d array of type unsigned char*/

unsigned char *UCalloc1d(int i_size)
{
	unsigned char *array;

	array=(unsigned char *) calloc(i_size,sizeof(unsigned char ));

	return(array);
}
  

/*Allocate 2d array of type unsigned char*/

unsigned char **UCalloc2d(int i_size,int j_size)
{
	unsigned char **array;
	int i;

	array=(unsigned char **) calloc(i_size,sizeof(unsigned char *));

	for(i=0;i<i_size;i++)
		array[i]=(unsigned char *) calloc(j_size,sizeof(unsigned char ));

	return(array);
}
  
/*Allocate 3d array of type unsigned char*/

unsigned char ***UCalloc3d(int i_size,int j_size,int k_size)
{
	unsigned char ***array;
	int i,k;

	array=(unsigned char ***) calloc(k_size,sizeof(unsigned char **));

	for(k=0;k<k_size;k++)
		array[k]=(unsigned char **) calloc(i_size,sizeof(unsigned char *));

	for(k=0;k<k_size;k++)
	for(i=0;i<i_size;i++)
		array[k][i]=(unsigned char *) calloc(j_size,sizeof(unsigned char ));
	
	return(array);
}

/*Allocate 4d array of type unsigned char*/

unsigned char ****UCalloc4d(int i_size,int j_size,int k_size, int t_size)
{
	unsigned char ****array;
	int i, k, t;

	array=(unsigned char ****) calloc(t_size,sizeof(unsigned char ***));

	for(t=0;t<t_size;t++)
		array[t]=(unsigned char ***) calloc(k_size,sizeof(unsigned char **));

	for(t=0;t<t_size;t++)
	for(k=0;k<k_size;k++)
		array[t][k]=(unsigned char **) calloc(i_size,sizeof(unsigned char *));

	for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
		array[t][k][i]=(unsigned char *) calloc(j_size,sizeof(unsigned char ));
	
	return(array);
}



/*Allocate 1d array of type int*/

int *Ialloc1d(int i_size)
{
	int *array;

	array=(int *) calloc(i_size,sizeof(int ));

	return(array);
}

/*Allocate 2d array of type int*/

int **Ialloc2d(int i_size, int j_size)
{
	int **array;
	int i;

	array=(int **) calloc(i_size,sizeof(int *));

	for(i=0;i<i_size;i++)
		array[i]=(int *) calloc(j_size,sizeof(int ));

	return(array);
}

/*Allocate 3d array of type int*/

int ***Ialloc3d(int i_size, int j_size, int k_size)
{
	int ***array;
	int i,k;

	array=(int ***) calloc(k_size,sizeof(int **));

	for(k=0;k<k_size;k++)
		array[k]=(int **) calloc(i_size,sizeof(int *));

	for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
		array[k][i]=(int *) calloc(j_size,sizeof(int ));
	
	return(array);
}

/*Allocate 4d array of type int*/

int ****Ialloc4d(int i_size, int j_size, int k_size, int t_size)
{
	int ****array;
	int i, k, t;

	array=(int ****) calloc(t_size,sizeof(int ***));

	for(t=0;t<t_size;t++)
		array[t]=(int ***) calloc(k_size,sizeof(int **));

	for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
		array[t][k]=(int **) calloc(i_size,sizeof(int *));

	for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
		array[t][k][i]=(int *) calloc(j_size,sizeof(int ));
	
	return(array);
}



/*Allocate 1d array of type short*/

short *Salloc1d(int i_size)
{
	short *array;

	array=(short *) calloc(i_size,sizeof(short ));

	return(array);
}

/*Allocate 2d array of type short*/

short **Salloc2d(int i_size, int j_size)
{
	short **array;
	int i;

	array=(short **) calloc(i_size,sizeof(short *));

	for(i=0;i<i_size;i++)
		array[i]=(short *) calloc(j_size,sizeof(short ));

	return(array);
}

/*Allocate 3d array of type short*/

short ***Salloc3d(int i_size, int j_size, int k_size)
{
	short ***array;
	int i, k;

	array=(short ***) calloc(k_size,sizeof(short **));

	for(k=0;k<k_size;k++)
		array[k]=(short **) calloc(i_size,sizeof(short *));

	for(k=0;k<k_size;k++)
	for(i=0;i<i_size;i++)
		array[k][i]=(short *) calloc(j_size,sizeof(short ));
	
	return(array);
}

/*Allocate 4d array of type short*/

short ****Salloc4d(int i_size,int j_size,int k_size, int t_size)
{
	short ****array;
	int i, k, t;

	array=(short ****) calloc(t_size,sizeof(short ***));

	for(t=0;t<t_size;t++)
		array[t]=(short ***) calloc(k_size,sizeof(short **));

	for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
		array[t][k]=(short **) calloc(i_size,sizeof(short *));

	for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
		array[t][k][i]=(short *) calloc(j_size,sizeof(short ));
	
	return(array);
}



/*Allocate 1d array of type long*/

long *Lalloc1d(int i_size)
{
	long *array;

	array=(long *) calloc(i_size,sizeof(long ));

	return(array);
}

/*Allocate 2d array of type long*/

long **Lalloc2d(int i_size,int j_size)
{
	long **array;
	int i;

	array=(long **) calloc(i_size,sizeof(long *));

	for(i=0;i<i_size;i++)
		array[i]=(long *) calloc(j_size,sizeof(long ));

	return(array);
}

/*Allocate 3d array of type long*/

long ***Lalloc3d(int i_size,int j_size,int k_size)
{
	long ***array;
	int i, k;

	array=(long ***) calloc(k_size,sizeof(long **));

	for(k=0;k<k_size;k++)
		array[k]=(long **) calloc(i_size,sizeof(long *));

	for(k=0;k<k_size;k++)
	for(i=0;i<i_size;i++)
		array[k][i]=(long *) calloc(j_size,sizeof(long ));
	
	return(array);
}

/*Allocate 4d array of type long*/

long ****Lalloc4d(int i_size,int j_size,int k_size, int t_size)
{
	long ****array;
	int i, k, t;

	array=(long ****) calloc(t_size,sizeof(long ***));

	for(t=0;t<t_size;t++)
		array[t]=(long ***) calloc(k_size,sizeof(long **));

	for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
		array[t][k]=(long **) calloc(i_size,sizeof(long *));

	for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
		array[t][k][i]=(long *) calloc(j_size,sizeof(long ));
	
	return(array);
}


/*Allocate 1d array of type float*/

float *Falloc1d(int i_size)
{
	float *array;

	array=(float *) calloc(i_size,sizeof(float ));

	return(array);
}

/*Allocate 2d array of type float*/

float **Falloc2d(int i_size, int j_size)
{
	float **array;
	int i;

	array=(float **) calloc(i_size,sizeof(float *));

	for(i=0;i<i_size;i++)
		array[i]=(float *) calloc(j_size,sizeof(float ));

	return(array);
}

/*Allocate 3d array of type float*/

float ***Falloc3d(int i_size, int j_size, int k_size)
{
	float ***array;
	int i, k;

	array=(float ***) calloc(k_size,sizeof(float **));

	for(k=0;k<k_size;k++)
		array[k]=(float **) calloc(i_size,sizeof(float *));

	for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
		array[k][i]=(float *) calloc(j_size,sizeof(float ));
	
	return(array);
}

/*Allocate 4d array of type float*/

float ****Falloc4d(int i_size, int j_size, int k_size, int t_size)
{
	float ****array;
	int i, k, t;

	array=(float ****) calloc(t_size,sizeof(float ***));

	for(t=0;t<t_size;t++)
		array[t]=(float ***) calloc(k_size,sizeof(float **));

	for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
		array[t][k]=(float **) calloc(i_size,sizeof(float *));

	for(t=0;t<t_size;t++)
	for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
		array[t][k][i]=(float *) calloc(j_size,sizeof(float ));
	
	return(array);
}



/*Allocate 1d array of type double*/

double *Dalloc1d(int i_size)
{
	double *array;

	array=(double *) calloc(i_size,sizeof(double ));

	return(array);
}

/*Allocate 2d array of type double*/

double **Dalloc2d(int i_size, int j_size)
{
	double **array;
	int i;

	array=(double **) calloc(i_size,sizeof(double *));

	for(i=0;i<i_size;i++)
		array[i]=(double *) calloc(j_size,sizeof(double ));

	return(array);
}

/*Allocate 3d array of type double*/

double ***Dalloc3d(int i_size, int j_size, int k_size)
{
	double ***array;
	int i, k;

	array=(double ***) calloc(k_size,sizeof(double **));

	for(k=0;k<k_size;k++)
		array[k]=(double **) calloc(i_size,sizeof(double *));

	for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
		array[k][i]=(double *) calloc(j_size,sizeof(double ));
	
	return(array);
}

/*Allocate 4d array of type double*/

double ****Dalloc4d(int i_size, int j_size, int k_size, int t_size) /*this one has been changed since June 27, 2002. */
{
	double ****array;
	int i,k, t;

	array=(double ****) calloc(t_size,sizeof(double ***));

	for(t=0;t<t_size;t++)
		array[t]=(double ***) calloc(k_size,sizeof(double **));

	for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
		array[t][k]=(double **) calloc(i_size,sizeof(double *));

	for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
		array[t][k][i]=(double *) calloc(j_size,sizeof(double ));
	
	return(array);
}



/*Allocate 1d array of type UCvector3d*/

UCvector3d *UCvector3dalloc1d(int i_size)
{
	UCvector3d *array;

	array=(UCvector3d *) calloc(i_size,sizeof(UCvector3d ));

	return(array);
}
  

/*Allocate 2d array of type UCvector3d*/

UCvector3d  **UCvector3dalloc2d(int i_size, int j_size)
{
	UCvector3d **array;
	int i;

	array=(UCvector3d **) calloc(i_size,sizeof(UCvector3d *));

	for(i=0;i<i_size;i++)
		array[i]=(UCvector3d *) calloc(j_size,sizeof(UCvector3d ));

	return(array);
}
  
/*Allocate 3d array of type UCvector3d*/

UCvector3d ***UCvector3dalloc3d(int i_size,int j_size,int k_size)
{
	UCvector3d ***array;
	int i,k;

	array=(UCvector3d ***) calloc(k_size,sizeof(UCvector3d **));

	for(k=0;k<k_size;k++)
		array[k]=(UCvector3d **) calloc(i_size,sizeof(UCvector3d *));

	for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
		array[k][i]=(UCvector3d *) calloc(j_size,sizeof(UCvector3d ));
	
	return(array);
}

/*Allocate 4d array of type UCvector3d*/

UCvector3d ****UCvector3dalloc4d(int i_size,int j_size,int k_size, int t_size)
{
	UCvector3d ****array;
	int i, k, t;

	array=(UCvector3d ****) calloc(t_size,sizeof(UCvector3d ***));

	for(t=0;t<t_size;t++)
		array[t]=(UCvector3d ***) calloc(k_size,sizeof(UCvector3d **));

	for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
		array[t][k]=(UCvector3d **) calloc(i_size,sizeof(UCvector3d *));

	for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
		array[t][k][i]=(UCvector3d *) calloc(j_size,sizeof(UCvector3d ));
	
	return(array);
}



/*Allocate 1d array of type UCvector2d*/

UCvector2d *UCvector2dalloc1d(int i_size)
{
  UCvector2d *array;

  array=(UCvector2d *) calloc(i_size,sizeof(UCvector2d ));

  return(array);
}

/*Allocate 2d array of type UCvector2d*/
  
UCvector2d  **UCvector2dalloc2d(int i_size,int j_size)
{
  UCvector2d **array;
  int i;

  array=(UCvector2d **) calloc(i_size,sizeof(UCvector2d *));

  for(i=0;i<i_size;i++)
    array[i]=(UCvector2d *) calloc(j_size,sizeof(UCvector2d ));

  return(array);
}

/*Allocate 3d array of type UCvector2d*/
  
UCvector2d ***UCvector2dalloc3d(int i_size,int j_size,int k_size)
{
  UCvector2d ***array;
  int i, k;

  array=(UCvector2d ***) calloc(k_size,sizeof(UCvector2d **));

  for(k=0;k<k_size;k++)
    array[k]=(UCvector2d **) calloc(i_size,sizeof(UCvector2d *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(UCvector2d *) calloc(j_size,sizeof(UCvector2d ));
	
  return(array);
}


/*Allocate 1d array of type Ivector3d*/

Ivector3d *Ivector3dalloc1d(int i_size)
{
  Ivector3d *array;

  array=(Ivector3d *) calloc(i_size,sizeof(Ivector3d ));

  return(array);
}
  

/*Allocate 2d array of type Ivector3d*/

Ivector3d  **Ivector3dalloc2d(int i_size,int j_size)
{
  Ivector3d **array;
  int i;

  array=(Ivector3d **) calloc(i_size,sizeof(Ivector3d *));

  for(i=0;i<i_size;i++)
    array[i]=(Ivector3d *) calloc(j_size,sizeof(Ivector3d ));

  return(array);
}
  
/*Allocate 3d array of type Ivector3d*/

Ivector3d ***Ivector3dalloc3d(int i_size,int j_size,int k_size)
{
  Ivector3d ***array;
  int i, k;

  array=(Ivector3d ***) calloc(k_size,sizeof(Ivector3d **));

  for(k=0;k<k_size;k++)
    array[k]=(Ivector3d **) calloc(i_size,sizeof(Ivector3d *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(Ivector3d *) calloc(j_size,sizeof(Ivector3d ));
	
  return(array);
}

/*Allocate 4d array of type Ivector3d*/

Ivector3d ****Ivector3dalloc4d(int i_size,int j_size,int k_size, int t_size)
{
  Ivector3d ****array;
  int i, k, t;

  array=(Ivector3d ****) calloc(t_size,sizeof(Ivector3d ***));

  for(t=0;t<t_size;t++)
    array[t]=(Ivector3d ***) calloc(k_size,sizeof(Ivector3d **));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      array[t][k]=(Ivector3d **) calloc(i_size,sizeof(Ivector3d *));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	array[t][k][i]=(Ivector3d *) calloc(j_size,sizeof(Ivector3d ));
	
  return(array);
}


/*Allocate 1d array of type Ivector2d*/

Ivector2d *Ivector2dalloc1d(int i_size)
{
  Ivector2d *array;

  array=(Ivector2d *) calloc(i_size,sizeof(Ivector2d ));

  return(array);
}

/*Allocate 2d array of type Ivector2d*/
  
Ivector2d  **Ivector2dalloc2d(int i_size,int j_size)
{
  Ivector2d **array;
  int i;

  array=(Ivector2d **) calloc(i_size,sizeof(Ivector2d *));

  for(i=0;i<i_size;i++)
    array[i]=(Ivector2d *) calloc(j_size,sizeof(Ivector2d ));

  return(array);
}

/*Allocate 3d array of type Ivector2d*/
  
Ivector2d ***Ivector2dalloc3d(int i_size,int j_size,int k_size)
{
  Ivector2d ***array;
  int i,k;

  array=(Ivector2d ***) calloc(k_size,sizeof(Ivector2d **));

  for(k=0;k<k_size;k++)
    array[k]=(Ivector2d **) calloc(i_size,sizeof(Ivector2d *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(Ivector2d *) calloc(j_size,sizeof(Ivector2d ));
	
  return(array);
}


/*Allocate 1d array of type Svector3d*/

Svector3d *Svector3dalloc1d(int i_size)
{
  Svector3d *array;

  array=(Svector3d *) calloc(i_size,sizeof(Svector3d ));

  return(array);
}
  

/*Allocate 2d array of type Svector3d*/

Svector3d  **Svector3dalloc2d(int i_size,int j_size)
{
  Svector3d **array;
  int i;

  array=(Svector3d **) calloc(i_size,sizeof(Svector3d *));

  for(i=0;i<i_size;i++)
    array[i]=(Svector3d *) calloc(j_size,sizeof(Svector3d ));

  return(array);
}
  
/*Allocate 3d array of type Svector3d*/

Svector3d ***Svector3dalloc3d(int i_size,int j_size,int k_size)
{
  Svector3d ***array;
  int i,k;

  array=(Svector3d ***) calloc(k_size,sizeof(Svector3d **));

  for(k=0;k<k_size;k++)
    array[k]=(Svector3d **) calloc(i_size,sizeof(Svector3d *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(Svector3d *) calloc(j_size,sizeof(Svector3d ));
	
  return(array);
}

/*Allocate 4d array of type Svector3d*/

Svector3d ****Svector3dalloc4d(int i_size,int j_size,int k_size, int t_size)
{
  Svector3d ****array;
  int i,k, t;

  array=(Svector3d ****) calloc(t_size,sizeof(Svector3d ***));

  for(t=0;t<t_size;t++)
    array[t]=(Svector3d ***) calloc(k_size,sizeof(Svector3d **));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      array[t][k]=(Svector3d **) calloc(i_size,sizeof(Svector3d *));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	array[t][k][i]=(Svector3d *) calloc(j_size,sizeof(Svector3d ));
	
  return(array);
}



/*Allocate 1d array of type Svector2d*/

Svector2d *Svector2dalloc1d(int i_size)
{
  Svector2d *array;

  array=(Svector2d *) calloc(i_size,sizeof(Svector2d ));

  return(array);
}

/*Allocate 2d array of type Svector2d*/
  
Svector2d  **Svector2dalloc2d(int i_size,int j_size)
{
  Svector2d **array;
  int i;

  array=(Svector2d **) calloc(i_size,sizeof(Svector2d *));

  for(i=0;i<i_size;i++)
    array[i]=(Svector2d *) calloc(j_size,sizeof(Svector2d ));

  return(array);
}

/*Allocate 3d array of type Svector2d*/
  
Svector2d ***Svector2dalloc3d(int i_size,int j_size,int k_size)
{
  Svector2d ***array;
  int i,k;

  array=(Svector2d ***) calloc(k_size,sizeof(Svector2d **));

  for(k=0;k<k_size;k++)
    array[k]=(Svector2d **) calloc(i_size,sizeof(Svector2d *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(Svector2d *) calloc(j_size,sizeof(Svector2d ));
	
  return(array);
}


/*Allocate 1d array of type Lvector3d*/

Lvector3d *Lvector3dalloc1d(int i_size)
{
  Lvector3d *array;

  array=(Lvector3d *) calloc(i_size,sizeof(Lvector3d ));

  return(array);
}
  

/*Allocate 2d array of type Lvector3d*/

Lvector3d  **Lvector3dalloc2d(int i_size,int j_size)
{
  Lvector3d **array;
  int i;

  array=(Lvector3d **) calloc(i_size,sizeof(Lvector3d *));

  for(i=0;i<i_size;i++)
    array[i]=(Lvector3d *) calloc(j_size,sizeof(Lvector3d ));

  return(array);
}
  
/*Allocate 3d array of type Lvector3d*/

Lvector3d ***Lvector3dalloc3d(int i_size,int j_size,int k_size)
{
  Lvector3d ***array;
  int i,k;

  array=(Lvector3d ***) calloc(k_size,sizeof(Lvector3d **));

  for(k=0;k<k_size;k++)
    array[k]=(Lvector3d **) calloc(i_size,sizeof(Lvector3d *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(Lvector3d *) calloc(j_size,sizeof(Lvector3d ));
	
  return(array);
}

/*Allocate 4d array of type Lvector3d*/

Lvector3d ****Lvector3dalloc4d(int i_size,int j_size,int k_size, int t_size)
{
  Lvector3d ****array;
  int i,k, t;

  array=(Lvector3d ****) calloc(t_size,sizeof(Lvector3d ***));

  for(t=0;t<t_size;t++)
    array[t]=(Lvector3d ***) calloc(k_size,sizeof(Lvector3d **));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      array[t][k]=(Lvector3d **) calloc(i_size,sizeof(Lvector3d *));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	array[t][k][i]=(Lvector3d *) calloc(j_size,sizeof(Lvector3d ));
	
  return(array);
}



/*Allocate 1d array of type Lvector2d*/

Lvector2d *Lvector2dalloc1d(int i_size)
{
  Lvector2d *array;

  array=(Lvector2d *) calloc(i_size,sizeof(Lvector2d ));

  return(array);
}

/*Allocate 2d array of type Lvector2d*/
  
Lvector2d  **Lvector2dalloc2d(int i_size,int j_size)
{
  Lvector2d **array;

  int i;

  array=(Lvector2d **) calloc(i_size,sizeof(Lvector2d *));

  for(i=0;i<i_size;i++)
    array[i]=(Lvector2d *) calloc(j_size,sizeof(Lvector2d ));

  return(array);
}

/*Allocate 3d array of type Lvector2d*/
  
Lvector2d ***Lvector2dalloc3d(int i_size,int j_size,int k_size)
{
  Lvector2d ***array;
  int i,k;

  array=(Lvector2d ***) calloc(k_size,sizeof(Lvector2d **));

  for(k=0;k<k_size;k++)
    array[k]=(Lvector2d **) calloc(i_size,sizeof(Lvector2d *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(Lvector2d *) calloc(j_size,sizeof(Lvector2d ));
	
  return(array);
}


/*Allocate 1d array of type Fvector4d*/

Fvector4d *Fvector4dalloc1d(int i_size)
{
  Fvector4d *array;

  array=(Fvector4d *) calloc(i_size,sizeof(Fvector4d ));

  return(array);
}


/*Allocate 1d array of type Fvector3d*/

Fvector3d *Fvector3dalloc1d(int i_size)
{
  Fvector3d *array;

  array=(Fvector3d *) calloc(i_size,sizeof(Fvector3d ));

  return(array);
}
  

/*Allocate 2d array of type Fvector3d*/

Fvector3d  **Fvector3dalloc2d(int i_size,int j_size)
{
  Fvector3d **array;
  int i;

  array=(Fvector3d **) calloc(i_size,sizeof(Fvector3d *));

  for(i=0;i<i_size;i++)
    array[i]=(Fvector3d *) calloc(j_size,sizeof(Fvector3d ));

  return(array);
}
  
/*Allocate 3d array of type Fvector3d*/

Fvector3d ***Fvector3dalloc3d(int i_size,int j_size,int k_size)
{
  Fvector3d ***array;
  int i,k;

  array=(Fvector3d ***) calloc(k_size,sizeof(Fvector3d **));

  for(k=0;k<k_size;k++)
    array[k]=(Fvector3d **) calloc(i_size,sizeof(Fvector3d *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(Fvector3d *) calloc(j_size,sizeof(Fvector3d ));
	
  return(array);
}

/*Allocate 4d array of type Fvector3d*/

Fvector3d ****Fvector3dalloc4d(int i_size,int j_size,int k_size, int t_size)
{
  Fvector3d ****array;
  int i, k, t;

  array=(Fvector3d ****) calloc(t_size,sizeof(Fvector3d ***));

  for(t=0;t<t_size;t++)
    array[t]=(Fvector3d ***) calloc(k_size,sizeof(Fvector3d **));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      array[t][k]=(Fvector3d **) calloc(i_size,sizeof(Fvector3d *));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	array[t][k][i]=(Fvector3d *) calloc(j_size,sizeof(Fvector3d ));
	
  return(array);
}


/*Allocate 1d array of type Fvector2d*/

Fvector2d *Fvector2dalloc1d(int i_size)
{
  Fvector2d *array;

  array=(Fvector2d *) calloc(i_size,sizeof(Fvector2d ));

  return(array);
}

/*Allocate 2d array of type Fvector2d*/
  
Fvector2d  **Fvector2dalloc2d(int i_size,int j_size)
{
  Fvector2d **array;
  int i;

  array=(Fvector2d **) calloc(i_size,sizeof(Fvector2d *));

  for(i=0;i<i_size;i++)
    array[i]=(Fvector2d *) calloc(j_size,sizeof(Fvector2d ));

  return(array);
}

/*Allocate 3d array of type Fvector2d*/
  
Fvector2d ***Fvector2dalloc3d(int i_size,int j_size,int k_size)
{
  Fvector2d ***array;
  int i,k;

  array=(Fvector2d ***) calloc(k_size,sizeof(Fvector2d **));

  for(k=0;k<k_size;k++)
    array[k]=(Fvector2d **) calloc(i_size,sizeof(Fvector2d *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(Fvector2d *) calloc(j_size,sizeof(Fvector2d ));
	
  return(array);
}

/*Allocate 1d array of type Dvector3d*/

Dvector3d *Dvector3dalloc1d(int i_size)
{
  Dvector3d *array;

  array=(Dvector3d *) calloc(i_size,sizeof(Dvector3d ));

  return(array);
}
  

/*Allocate 2d array of type Dvector3d*/

Dvector3d  **Dvector3dalloc2d(int i_size,int j_size)
{
  Dvector3d **array;
  int i;

  array=(Dvector3d **) calloc(i_size,sizeof(Dvector3d *));

  for(i=0;i<i_size;i++)
    array[i]=(Dvector3d *) calloc(j_size,sizeof(Dvector3d ));

  return(array);
}
  
/*Allocate 3d array of type Dvector3d*/

Dvector3d ***Dvector3dalloc3d(int i_size,int j_size,int k_size)
{
  Dvector3d ***array;
  int i, k;

  array=(Dvector3d ***) calloc(k_size,sizeof(Dvector3d **));

  for(k=0;k<k_size;k++)
    array[k]=(Dvector3d **) calloc(i_size,sizeof(Dvector3d *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(Dvector3d *) calloc(j_size,sizeof(Dvector3d ));
	
  return(array);
}

/*Allocate 4d array of type Dvector3d*/

Dvector3d ****Dvector3dalloc4d(int i_size,int j_size,int k_size, int t_size)
{
  Dvector3d ****array;
  int i, k, t;

  array=(Dvector3d ****) calloc(t_size,sizeof(Dvector3d ***));

  for(t=0;t<t_size;t++)
    array[t]=(Dvector3d ***) calloc(k_size,sizeof(Dvector3d **));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      array[t][k]=(Dvector3d **) calloc(i_size,sizeof(Dvector3d *));

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	array[t][k][i]=(Dvector3d *) calloc(j_size,sizeof(Dvector3d ));
	
  return(array);
}



/*Allocate 1d array of type Dvector2d*/

Dvector2d *Dvector2dalloc1d(int i_size)
{
  Dvector2d *array;

  array=(Dvector2d *) calloc(i_size,sizeof(Dvector2d ));

  return(array);
}

/*Allocate 2d array of type Dvector2d*/
  
Dvector2d  **Dvector2dalloc2d(int i_size,int j_size)
{
  Dvector2d **array;
  int i;

  array=(Dvector2d **) calloc(i_size,sizeof(Dvector2d *));

  for(i=0;i<i_size;i++)
    array[i]=(Dvector2d *) calloc(j_size,sizeof(Dvector2d ));

  return(array);
}

/*Allocate 3d array of type Dvector2d*/
  
Dvector2d ***Dvector2dalloc3d(int i_size,int j_size,int k_size)
{
  Dvector2d ***array;
  int i,k;

  array=(Dvector2d ***) calloc(k_size,sizeof(Dvector2d **));

  for(k=0;k<k_size;k++)
    array[k]=(Dvector2d **) calloc(i_size,sizeof(Dvector2d *));

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      array[k][i]=(Dvector2d *) calloc(j_size,sizeof(Dvector2d ));
	
  return(array);
}

/*Allocate 1d array of type Fsphere*/

Fsphere *Fspherealloc1d(int i_size)
{
  Fsphere *array;

  array=(Fsphere *) calloc(i_size,sizeof(Fsphere ));

  return(array);
}


/*****MEMORY FREEING*****/

/*Free 2d array of type char*/

void Cfree2d(char **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type char*/

void Cfree3d(char ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 2d array of type unsigned char*/

void UCfree2d(unsigned char **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type unsigned char*/

void UCfree3d(unsigned char ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 4d array of type unsigned char*/

void UCfree4d(unsigned char ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}


/*Free 2d array of type int*/

void Ifree2d(int **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}

/*Free 3d array of type int*/

void Ifree3d(int ***array,int k_size,int i_size)
{
  int k,i;


  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);
	
  free(array);
}

/*Free 4d array of type int*/

void Ifree4d(int ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}


/*Free 2d array of type short*/

void Sfree2d(short **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type short*/

void Sfree3d(short ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 4d array of type unsigned short*/

void Sfree4d(short ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}


/*Free 2d array of type long*/

void Lfree2d(long **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type long*/

void Lfree3d(long ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 4d array of type long*/

void Lfree4d(long ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}


/*Free 2d array of type float*/

void Ffree2d(float **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type float*/

void Ffree3d(float ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 4d array of type float*/

void Ffree4d(float ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}

  
/*Free 2d array of type double*/

void Dfree2d(double **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type double*/

void Dfree3d(double ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 4d array of type double*/

void Dfree4d(double ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}


/*Free 2d array of type UCvector3d*/

void UCvector3dfree2d(UCvector3d **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type UCvector3d*/

void UCvector3dfree3d(UCvector3d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 4d array of type UCvector3d*/

void UCvector3dfree4d(UCvector3d ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}


/*Free 2d array of type UCvector2d*/

void UCvector2dfree2d(UCvector2d **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type UCvector2d*/

void UCvector2dfree3d(UCvector2d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 2d array of type Ivector2d*/

void Ivector2dfree2d(Ivector2d **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type Ivector2d*/

void Ivector2dfree3d(Ivector2d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}
  

/*Allocate 1d array of type Ivector4d*/

Ivector4d *Ivector4dalloc1d(int i_size)
{
  Ivector4d *array;

  array=(Ivector4d *) calloc(i_size,sizeof(Ivector4d ));

  return(array);
}


/*Free 2d array of type Ivector3d*/

void Ivector3dfree2d(Ivector3d **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type Ivector3d*/

void Ivector3dfree3d(Ivector3d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 4d array of type Ivector3d*/

void Ivector3dfree4d(Ivector3d ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}


/*Free 2d array of type Svector3d*/

void Svector3dfree2d(Svector3d **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type Svector3d*/

void Svector3dfree3d(Svector3d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 4d array of type Svector3d*/

void Svector3dfree4d(Svector3d ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}


/*Free 2d array of type Svector2d*/

void Svector2dfree2d(Svector2d **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type Svector2d*/

void Svector2dfree3d(Svector2d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}
  
/*Free 2d array of type Lvector3d*/

void Lvector3dfree2d(Lvector3d **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type Lvector3d*/

void Lvector3dfree3d(Lvector3d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 4d array of type Lvector3d*/

void Lvector3dfree4d(Lvector3d ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}


/*Free 2d array of type Lvector2d*/

void Lvector2dfree2d(Lvector2d **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type Lvector2d*/

void Lvector2dfree3d(Lvector2d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 2d array of type Fvector3d*/

void Fvector3dfree2d(Fvector3d **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type Fvector3d*/

void Fvector3dfree3d(Fvector3d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 4d array of type Fvector3d*/

void Fvector3dfree4d(Fvector3d ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}


/*Free 2d array of type Fvector2d*/

void Fvector2dfree2d(Fvector2d **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type Fvector2d*/

void Fvector2dfree3d(Fvector2d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 2d array of type Dvector3d*/

void Dvector3dfree2d(Dvector3d **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}
  
/*Free 3d array of type Dvector3d*/

void Dvector3dfree3d(Dvector3d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/*Free 4d array of type Dvector3d*/

void Dvector3dfree4d(Dvector3d ****array,int t_size,int k_size,int i_size)
{
  int t,k,i;

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      for(i=0;i<i_size;i++)
	free(array[t][k][i]);

  for(t=0;t<t_size;t++)
    for(k=0;k<k_size;k++)
      free(array[t][k]);

  for(t=0;t<t_size;t++)
    free(array[t]);

  free(array);
}


/*Free 2d array of type Dvector2d*/

void Dvector2dfree2d(Dvector2d **array,int i_size)
{
  int i;

  for(i=0;i<i_size;i++)
    free(array[i]);

  free(array);
}

/*Free 3d array of type Dvector2d*/

void Dvector2dfree3d(Dvector2d ***array,int k_size,int i_size)
{
  int k,i;

  for(k=0;k<k_size;k++)
    for(i=0;i<i_size;i++)
      free(array[k][i]);

  for(k=0;k<k_size;k++)
    free(array[k]);

  free(array);
}

/******VARIOUS ROUTINES******/

/******
Redistributes a 1d array of type Fvector2d
Notes: M = output # of points
       N = input  # of points
       MAKE SURE that p has been allocated for at least M points of type Fvector3d
******/

void Fvector2dredistribute1d(Fvector2d *p,int M,int N)
{
  int i,j,done,max_size;
  double length,qstep,temp,dist,*len;
  Fvector2d *q;

  if(N>M)
    max_size=N;
  else
    max_size=M;

  len=Dalloc1d(max_size);
  q=Fvector2dalloc1d(max_size);

  for(i=0;i<M;i++)
    len[i]=0.0;

  len[0]=0;
  length=0;
  for(j=0;j<(N-1);j++)
    {
      length+=sqrt((double) ((p[(j+1)].x-p[j].x)*(p[(j+1)].x-p[j].x)+(p[(j+1)].y-p[j].y)*(p[(j+1)].y-p[j].y)));
      len[(j+1)]=length;
    }
  qstep=(length-.000001)/((double) (M-1));
  
  q[0].x=p[0].x;
  q[0].y=p[0].y;

  dist=0;
  
  for(i=1;i<M;i++)
    {
      j=0;
      done=0;
      temp=0;
      while(done==0)
	{
	  j++;
	  temp=len[j];
	  if(temp>=((double) i)*qstep)
	    {
	      temp=(((double) i)*qstep-len[j-1])/(len[j]-len[j-1]);
	      q[i].x=p[j-1].x+temp*(p[j].x-p[j-1].x);
	      q[i].y=p[j-1].y+temp*(p[j].y-p[j-1].y);
	      done=1;
	    }
	}
    }

  for(i=0;i<M;i++)
    {
      p[i].x=q[i].x;
      p[i].y=q[i].y;
    }

  free(len);
  free(q);
}

/*Redistributes a 1d array of type Fvector3d*/

void Fvector3dredistribute1d(Fvector3d *p,int M,int N)
{
  int i,j,done,max_size;
  double length,qstep,temp,dist,*len;
  Fvector3d *q;

  if(N>1)
    {
      if(N>M)
	max_size=N;
      else
	max_size=M;
      /*
      printf("N=%d M=%d max_size = %d",N,M,max_size);
      fflush(stdout);
      */
      len=Dalloc1d(max_size);
      q=Fvector3dalloc1d(max_size);
      
      for(i=0;i<M;i++)
	len[i]=0.0;
      
      len[0]=0;
      length=0;

      for(j=0;j<(N-1);j++)
	{
	  length+=sqrt((double) ((p[(j+1)].x-p[j].x)*(p[(j+1)].x-p[j].x)+(p[(j+1)].y-p[j].y)*(p[(j+1)].y-p[j].y)+(p[(j+1)].z-p[j].z)*(p[(j+1)].z-p[j].z)));
	  len[(j+1)]=length;
	}
      if(length<=0)
	printf("Warning:length <= 0\n");

      qstep=(length-.000001)/((double) (M-1));
      
      q[0].x=p[0].x;
      q[0].y=p[0].y;
      q[0].z=p[0].z;
      
      dist=0;
      
      for(i=1;i<M;i++)
	{
	  j=0;
	  done=0;
	  temp=0;
	  while(done==0)
	    {
	      j++;
	      temp=len[j];
	      if(temp>=((double) i)*qstep)
		{
		  temp=(((double) i)*qstep-len[j-1])/(len[j]-len[j-1]);
		  q[i].x=p[j-1].x+temp*(p[j].x-p[j-1].x);
		  q[i].y=p[j-1].y+temp*(p[j].y-p[j-1].y);
		  q[i].z=p[j-1].z+temp*(p[j].z-p[j-1].z);
		  done=1;
		}
	    }
	}
      
      for(i=0;i<M;i++)
	{
	  p[i].x=q[i].x;
	  p[i].y=q[i].y;
	  p[i].z=q[i].z;
	}
      
      free(len);
      free(q);
    }
  else
    {
      for(i=1;i<M;i++)
	{
	  p[i].x=p[i-1].x+.0001;
	  p[i].y=p[i-1].y+.0001;
	  p[i].z=p[i-1].z+.0001;
	}
    }
}

void Fvector3dinterpCurve1d(Fvector3d *v,Fvector3d *n,float *c,int M,int N)
{
  int i,j,done,max_size;
  double length,newstep,temp,dist,*len,mag;
  Fvector3d *vnew,*nnew;
  float *cnew;

  if(N>1)
    {
      if(N>M)
	max_size=N;
      else
	max_size=M;
      /*
      printf("N=%d M=%d max_size = %d",N,M,max_size);
      fflush(stdout);
      */
      len=Dalloc1d(max_size);
      vnew=Fvector3dalloc1d(max_size);
      nnew=Fvector3dalloc1d(max_size);
      cnew=Falloc1d(max_size);

      for(i=0;i<max_size;i++)
	len[i]=0.0;
      
      len[0]=0;
      length=0;

      for(j=0;j<(N-1);j++)
	{
	  length+=sqrt((double) ((v[(j+1)].x-v[j].x)*(v[(j+1)].x-v[j].x)+(v[(j+1)].y-v[j].y)*(v[(j+1)].y-v[j].y)+(v[(j+1)].z-v[j].z)*(v[(j+1)].z-v[j].z)));
	  len[(j+1)]=length;
	}
      if(length==0)
	printf("Warning:length = 0\n");

      newstep=(length-.000001)/((double) (M-1));
      
      vnew[0].x=v[0].x;
      vnew[0].y=v[0].y;
      vnew[0].z=v[0].z;
      nnew[0].x=n[0].x;
      nnew[0].y=n[0].y;
      nnew[0].z=n[0].z;
      cnew[0]=c[0];
      
      dist=0;
      
      for(i=1;i<M;i++)
	{
	  j=0;
	  done=0;
	  temp=0;
	  while(done==0)
	    {
	      j++;
	      temp=len[j];
	      if(temp>=((double) i)*newstep)
		{
		  temp=(((double) i)*newstep-len[j-1])/(len[j]-len[j-1]);
		  vnew[i].x=v[j-1].x+temp*(v[j].x-v[j-1].x);
		  vnew[i].y=v[j-1].y+temp*(v[j].y-v[j-1].y);
		  vnew[i].z=v[j-1].z+temp*(v[j].z-v[j-1].z);
		  nnew[i].x=n[j-1].x+temp*(n[j].x-n[j-1].x);
		  nnew[i].y=n[j-1].y+temp*(n[j].y-n[j-1].y);
		  nnew[i].z=n[j-1].z+temp*(n[j].z-n[j-1].z);
		  cnew[i]=c[j-1]+temp*(c[j]-c[j-1]);
		  done=1;
		}
	    }
	}

      for(i=0;i<M;i++)
	{
	  mag=Fvector3dmag(nnew[i]);
	  nnew[i].x/=mag;
	  nnew[i].y/=mag;
	  nnew[i].z/=mag;
	}
      
      for(i=0;i<M;i++)
	{
	  v[i].x=vnew[i].x;
	  v[i].y=vnew[i].y;
	  v[i].z=vnew[i].z;
	  n[i].x=nnew[i].x;
	  n[i].y=nnew[i].y;
	  n[i].z=nnew[i].z;
	  c[i]=cnew[i];
	}
      
      free(len);
      free(vnew);
      free(nnew);
      free(cnew);
    }
  else
    {
      for(i=1;i<M;i++)
	{
	  v[i].x=v[i-1].x+.0001;
	  v[i].y=v[i-1].y+.0001;
	  v[i].z=v[i-1].z+.0001;
	}
    }
}

/*Calculate Normals*/

void Fvector3dnorms2d(int NNN,int MMM,Fvector3d **p,Fvector3d **norm)
{
  int i,j;
  float mag;

  for(i=1;i<(NNN-1);i++)
    for(j=1;j<(MMM-1);j++)
      {
	norm[i][j].x=(p[i+1][j].y-p[i][j].y)*(p[i][j+1].z-p[i][j].z)-(p[i+1][j].z-p[i][j].z)*(p[i][j+1].y-p[i][j].y);
	norm[i][j].y=-(p[i+1][j].x-p[i][j].x)*(p[i][j+1].z-p[i][j].z)+(p[i+1][j].z-p[i][j].z)*(p[i][j+1].x-p[i][j].x);
	norm[i][j].z=(p[i+1][j].x-p[i][j].x)*(p[i][j+1].y-p[i][j].y)-(p[i+1][j].y-p[i][j].y)*(p[i][j+1].x-p[i][j].x);
	mag=Fvector3dmag(norm[i][j]);
	norm[i][j].x/=mag;
	norm[i][j].y/=mag;
	norm[i][j].z/=mag;
      }	

}


/*****Cartesian to spherical coordinate transformation*****/

void cartesianToSpherical(Fsphere *dataout,Fvector3d *datain,int N)
{
	int i;
	/*
	for(i=0;i<N;i++)
    {
		dataout[i].r=Fvector3dmag(datain[i]);
		dataout[i].t=((float) (acos(((double) (datain[i].z))/((double) (dataout[i].r)))*radtodeg));
		dataout[i].p=((float) (atan2(((double) (datain[i].y)),((double) (datain[i].x)))*radtodeg));
    }
	*/
}  
