/**************************************
*** Header file for the mvcd library***
***************************************/

#define degtorad M_PI/180
#define radtodeg 180/M_PI

typedef struct 
{
  unsigned char x;
  unsigned char y;
} UCvector2d;

typedef struct
{
  unsigned char x;
  unsigned char y;
  unsigned char z;
} UCvector3d;

typedef struct 
{
  int x;
  int y;
} Ivector2d;

typedef struct
{
  int x;
  int y;
  int z;
} Ivector3d;

typedef struct
{
  int x;
  int y;
  int z;
  int t;
} Ivector4d;

typedef struct 
{
  short x;
  short y;
} Svector2d;

typedef struct
{
  short x;
  short y;
  short z;
} Svector3d;

typedef struct 
{
  long x;
  long y;
} Lvector2d;

typedef struct
{
  long x;
  long y;
  long z;
} Lvector3d;

typedef struct 
{
  float x;
  float y;
} Fvector2d;

typedef struct
{
  float x;
  float y;
  float z;
} Fvector3d;

typedef struct
{
  float x;
  float y;
  float z;
  float t;
} Fvector4d;

typedef struct
{
  double x;
  double y;
} Dvector2d;

typedef struct
{
  double x;
  double y;
  double z;
} Dvector3d;

typedef struct
{
  float r;
  float t;
  float p;
} Fsphere;

/*****Memory Allocation*****/

char *Calloc1d(int);
char **Calloc2d(int,int);
char ***Calloc3d(int,int,int);
unsigned char *UCalloc1d(int);
unsigned char **UCalloc2d(int,int);
unsigned char ***UCalloc3d(int,int,int);
unsigned char ****UCalloc4d(int,int,int,int);
int *Ialloc1d(int);
int **Ialloc2d(int,int);
int ***Ialloc3d(int,int,int);
int ****Ialloc4d(int,int,int,int);
short *Salloc1d(int);
short **Salloc2d(int,int);
short ***Salloc3d(int,int,int);
short ****Salloc4d(int,int,int,int);
long *Lalloc1d(int);
long **Lalloc2d(int,int);
long ***Lalloc3d(int,int,int);
long ****Lalloc4d(int,int,int,int);
float *Falloc1d(int);
float **Falloc2d(int,int);
float ***Falloc3d(int,int,int);
float ****Falloc4d(int,int,int,int);
double *Dalloc1d(int);
double **Dalloc2d(int,int);
double ***Dalloc3d(int,int,int);
double ****Dalloc4d(int,int,int,int);
UCvector3d *UCvector3dalloc1d(int);
UCvector3d **UCvector3dalloc2d(int,int);
UCvector3d ***UCvector3dalloc3d(int,int,int);
UCvector3d ****UCvector3dalloc4d(int,int,int,int);
UCvector2d *UCvector2dalloc1d(int);
UCvector2d **UCvector2dalloc2d(int,int);
UCvector2d ***UCvector2dalloc3d(int,int,int);
Ivector3d *Ivector3dalloc1d(int);
Ivector3d **Ivector3dalloc2d(int,int);
Ivector3d ***Ivector3dalloc3d(int,int,int);
Ivector3d ****Ivector3dalloc4d(int,int,int,int);
Ivector2d *Ivector2dalloc1d(int);
Ivector2d **Ivector2dalloc2d(int,int);
Ivector2d ***Ivector2dalloc3d(int,int,int);
Svector3d *Svector3dalloc1d(int);
Svector3d **Svector3dalloc2d(int,int);
Svector3d ***Svector3dalloc3d(int,int,int);
Svector3d ****Svector3dalloc4d(int,int,int,int);
Svector2d *Svector2dalloc1d(int);
Svector2d **Svector2dalloc2d(int,int);
Svector2d ***Svector2dalloc3d(int,int,int);
Lvector3d *Lvector3dalloc1d(int);
Lvector3d **Lvector3dalloc2d(int,int);
Lvector3d ***Lvector3dalloc3d(int,int,int);
Lvector3d ****Lvector3dalloc4d(int,int,int,int);
Lvector2d *Lvector2dalloc1d(int);
Lvector2d **Lvector2dalloc2d(int,int);
Lvector2d ***Lvector2dalloc3d(int,int,int);
Fvector3d *Fvector3dalloc1d(int);
Fvector3d **Fvector3dalloc2d(int,int);
Fvector3d ***Fvector3dalloc3d(int,int,int);
Fvector3d ****Fvector3dalloc4d(int,int,int,int);
Fvector2d *Fvector2dalloc1d(int);
Fvector2d **Fvector2dalloc2d(int,int);
Fvector2d ***Fvector2dalloc3d(int,int,int);
Dvector3d *Dvector3dalloc1d(int);
Dvector3d **Dvector3dalloc2d(int,int);
Dvector3d ***Dvector3dalloc3d(int,int,int);
Dvector3d ****Dvector3dalloc4d(int,int,int,int);
Dvector2d *Dvector2dalloc1d(int);
Dvector2d **Dvector2dalloc2d(int,int);
Dvector2d ***Dvector2dalloc3d(int,int,int);
Fsphere *Fspherealloc1d(int);
Ivector4d *Ivector4dalloc1d(int i_size);
Fvector4d *Fvector4dalloc1d(int i_size);


void Cfree2d(char **,int);
void Cfree3d(char ***,int,int);
void UCfree2d(unsigned char **,int);
void UCfree3d(unsigned char ***,int,int);
void UCfree4d(unsigned char ****,int,int,int);
void Ifree2d(int **,int);
void Ifree3d(int ***,int,int);
void Ifree4d(int ****,int,int,int);
void Sfree2d(short **,int);
void Sfree3d(short ***,int,int);
void Sfree4d(short ****,int,int,int);
void Lfree2d(long **,int);
void Lfree3d(long ***,int,int);
void Lfree4d(long ****,int,int,int);
void Ffree2d(float **,int);
void Ffree3d(float ***,int,int);
void Ffree4d(float ****,int,int,int);
void Dfree2d(double **,int);
void Dfree3d(double ***,int,int);
void Dfree4d(double ****,int,int,int);
void UCvector3dfree2d(UCvector3d **,int);
void UCvector3dfree3d(UCvector3d ***,int,int);
void UCvector3dfree4d(UCvector3d ****,int,int,int);
void UCvector2dfree2d(UCvector2d **,int);
void UCvector2dfree3d(UCvector2d ***,int,int);
void Ivector3dfree2d(Ivector3d **,int);
void Ivector3dfree3d(Ivector3d ***,int,int);
void Ivector3dfree4d(Ivector3d ****,int,int,int);
void Ivector2dfree2d(Ivector2d **,int);
void Ivector2dfree3d(Ivector2d ***,int,int);
void Svector3dfree2d(Svector3d **,int);
void Svector3dfree3d(Svector3d ***,int,int);
void Svector3dfree4d(Svector3d ****,int,int,int);
void Svector2dfree2d(Svector2d **,int);
void Svector2dfree3d(Svector2d ***,int,int);
void Lvector3dfree2d(Lvector3d **,int);
void Lvector3dfree3d(Lvector3d ***,int,int);
void Lvector3dfree4d(Lvector3d ****,int,int,int);
void Lvector2dfree2d(Lvector2d **,int);
void Lvector2dfree3d(Lvector2d ***,int,int);
void Fvector3dfree2d(Fvector3d **,int);
void Fvector3dfree3d(Fvector3d ***,int,int);
void Fvector3dfree4d(Fvector3d ****,int,int,int);
void Fvector2dfree2d(Fvector2d **,int);
void Fvector2dfree3d(Fvector2d ***,int,int);
void Dvector3dfree2d(Dvector3d **,int);
void Dvector3dfree3d(Dvector3d ***,int,int);
void Dvector3dfree4d(Dvector3d ****,int,int,int);
void Dvector2dfree2d(Dvector2d **,int);
void Dvector2dfree3d(Dvector2d ***,int,int);

/*****Cartesian to spherical coordinate transformation*****/

void cartesianToSpherical(Fsphere *,Fvector3d *,int);
/*****Error reporting file open function****/

FILE *myopen(char *,char *);

/*****Redistribute Points*****/

void Fvector2dredistribute1d(Fvector2d *,int,int);
void Fvector3dredistribute1d(Fvector3d *,int,int);

void Fvector3dinterpCurve1d(Fvector3d *,Fvector3d *,float *,int,int);


/*****Compute magnitude of 2D and 3D vectors*****/

double Fvector2dmag(Fvector2d);
double Fvector3dmag(Fvector3d);

/*****Compute Normals*****/

void Fvector3dnorms2d(int,int,Fvector3d **,Fvector3d **);

/*****
Note: Fvector3dcurvature2d() computes curvature and is in the cres library
*****/


