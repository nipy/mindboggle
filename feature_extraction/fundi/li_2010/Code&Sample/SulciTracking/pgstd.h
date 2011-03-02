/************************************************/
/*						*/
/*	File: pgstd.h	  		        */
/*      Purpose: Header file for pgstd library  */
/*						*/
/************************************************/

/*
 * Some old RCS info.
 * Revision 1.10  1998/07/24  17:39:52  pham
 * Started 4-D array allocation functions, but not finished.
**/
#ifndef PGSTD_TOOLS
#define PGSTD_TOOLS

/* PG Constants */
#define PG_ZERO    0
#define PG_BOUNDARY 1
#define PG_EXTRAPOLATE 2

#define PG_PI      3.141592654

#define PG_TRUE    1
#define PG_FALSE   0

#define PG_OK      0
#define PG_ERROR  -1

/* PG Macros */
#define PG_MAX(a,b)  (a>b? a: b)
#define PG_MIN(a,b)  (a<b? a: b)
#define PG_SQR(a)    ((a)*(a))
#define PG_TOGGLE(a) a= (a==0)? 1: 0

/* Pseudo type definitions */
#define PGuchar PGbyte

/* PG Type definitions */
typedef int PGerror;

typedef unsigned char  PGbyte;
typedef          char  PGchar;
typedef unsigned short PGushort;
typedef          short PGshort; 
typedef unsigned int   PGuint;
typedef          int   PGint;
typedef float          PGfloat;
typedef double         PGdouble;
typedef void           PGvoid;

typedef struct { PGbyte r;
                 PGbyte g;
                 PGbyte b;
               } PGrgb;

typedef struct { PGbyte x;
                 PGbyte y;
               } PGbyte2d;

typedef struct { PGint x;
                 PGint y;
               } PGint2d;

typedef struct { PGfloat x;
                 PGfloat y;
               } PGfloat2d;

typedef struct { PGdouble x;
                 PGdouble y;
	       } PGdouble2d;

typedef struct { PGbyte x;
                 PGbyte y;
                 PGbyte z;
               } PGbyte3d;

typedef struct { PGint x;
                 PGint y;
                 PGint z;
               } PGint3d;

typedef struct { PGfloat x;
                 PGfloat y;
                 PGfloat z;
               } PGfloat3d;

typedef struct { PGdouble x;
                 PGdouble y;
                 PGdouble z;
	       } PGdouble3d;

typedef struct { PGdouble3d x;
                 PGdouble3d y;
                 PGdouble3d z;
               } PGdtensor3d;

/* Standard error routine */
PGvoid pgError();

/* Misc routines */
void *pgMalloc(int size);
void *pgRealloc(void *ptr, int size);

/* Allocation and freeing routines from */
/* Vector routines - dynamically allocates one dimensional arrays */
PGbyte      *pgBvector();
PGvoid       pgFreeBvector();
PGushort    *pgUsvector();
PGvoid       pgFreeUsvector();
PGshort     *pgSvector();
PGvoid       pgFreeSvector();
PGuint      *pgUivector();
PGvoid       pgFreeUivector();
PGint       *pgIvector();
PGvoid       pgFreeIvector();
PGfloat     *pgFvector();
PGvoid       pgFreeFvector();
PGdouble    *pgDvector();
PGvoid       pgFreeDvector();
PGbyte2d    *pgB2dvector();       
PGvoid       pgFreeB2dvector();   
PGint2d     *pgI2dvector();       
PGvoid       pgFreeI2dvector();   
PGfloat2d   *pgF2dvector();       
PGvoid       pgFreeF2dvector();   
PGdouble2d  *pgD2dvector();
PGvoid       pgFreeD2dvector();
PGbyte3d    *pgB3dvector();        
PGvoid       pgFreeB3dvector();   
PGint3d     *pgI3dvector();     
PGvoid       pgFreeI3dvector(); 
PGfloat3d   *pgF3dvector();
PGvoid       pgFreeF3dvector();
PGdouble3d  *pgD3dvector();
PGvoid       pgFreeD3dvector();

/* Matrix Routines - dynamically allocates two dimensional arrays*/
PGbyte      **pgBmatrix();
PGvoid        pgFreeBmatrix();
PGint       **pgImatrix();
PGvoid        pgFreeImatrix();
PGushort    **pgUsmatrix();
PGvoid        pgFreeUsmatrix();
PGshort     **pgSmatrix();
PGvoid        pgFreeSmatrix();
PGuint      **pgUimatrix();
PGvoid        pgFreeUimatrix();
PGfloat     **pgFmatrix();
PGvoid        pgFreeFmatrix();
PGdouble    **pgDmatrix();
PGvoid        pgFreeDmatrix();
PGbyte2d    **pgB2dmatrix();      
PGvoid        pgFreeB2dmatrix();  
PGint2d     **pgI2dmatrix();      
PGvoid        pgFreeI2dmatrix();   
PGfloat2d   **pgF2dmatrix();
PGvoid        pgFreeF2dmatrix();
PGdouble2d  **pgD2dmatrix();
PGvoid        pgFreeD2dmatrix();
PGbyte3d    **pgB3dmatrix();    
PGvoid        pgFreeB3dmatrix(); 
PGint3d     **pgI3dmatrix();     
PGvoid        pgFreeI3dmatrix(); 
PGfloat3d   **pgF3dmatrix();     
PGvoid        pgFreeF3dmatrix(); 
PGdouble3d  **pgD3dmatrix();
PGvoid        pgFreeD3dmatrix();

/* Cube Routines - dynamically allocate three dimensional arrays*/
PGbyte      ***pgBcube();
PGvoid         pgFreeBcube();
PGint       ***pgIcube();
PGvoid         pgFreeIcube();
PGushort    ***pgUscube();
PGvoid         pgFreeUscube();
PGshort     ***pgScube();
PGvoid         pgFreeScube();
PGfloat     ***pgFcube();
PGvoid         pgFreeFcube();
PGdouble    ***pgDcube();
PGvoid         pgFreeDcube();
PGbyte2d    ***pgB2dcube(); 
PGvoid         pgFreeB2dcube(); 
PGint2d     ***pgI2dcube();     
PGvoid         pgFreeI2dcube(); 
PGfloat2d   ***pgF2dcube();     
PGvoid         pgFreeF2dcube(); 
PGdouble2d  ***pgD2dcube();
PGvoid         pgFreeD2dcube();
PGbyte3d    ***pgB3dcube();     
PGvoid         pgFreeB3dcube(); 
PGint3d     ***pgI3dcube();     
PGvoid         pgFreeI3dcube(); 
PGfloat3d   ***pgF3dcube();
PGvoid         pgFreeF3dcube();
PGdouble3d  ***pgD3dcube();
PGvoid         pgFreeD3dcube();

/*  four Routines - dynamically allocate four dimensional arrays*/
PGbyte     ****pgBfour();
PGvoid         pgFreeBfour();
PGint      ****pgIfour();
PGvoid         pgFreeIfour();
PGushort   ****pgUsfour();
PGvoid         pgFreeUsfour();
PGshort    ****pgSfour();
PGvoid         pgFreeSfour();
PGfloat    ****pgFfour();
PGvoid         pgFreeFfour();
PGdouble   ****pgDfour();
PGvoid         pgFreeDfour();

/* File I/O Routines */
int pgReadByte();
int pgGetByte();
int pgReadUshort();
int pgGetUshort();
int pgReadShort();
int pgGetShort();
int pgReadUint();
int pgGetUint();
int pgReadInt();
int pgGetInt();
int pgReadFloat();
int pgGetFloat();
int pgReadDouble();
int pgGetDouble();

int pgWriteByte();
int pgPrintByte();
int pgWriteUshort();
int pgPrintUshort();
int pgWriteShort();
int pgPrintShort();
int pgWriteUint();
int pgPrintUint();
int pgWriteInt();
int pgPrintInt();
int pgWriteFloat();
int pgPrintFloat();
int pgWriteDouble();
int pgPrintDouble();

void pgGetNextLine();
int pgFileSize();
int pgFileExists();

/* DX interface routines */
int pgWriteCubeDXhead(char*, int, int, int, PGfloat3d*, PGfloat3d*, 
		      char*, int, int, char*);


/* PGM I/O routines */
PGbyte **pgPGMRead();
int pgPGMWrite();

/* math and image processing routines */

void pgBcubeGradient(PGfloat3d ***g, PGbyte  ***I, int ZN, int YN, int XN);
void pgBcubeCopy(PGbyte ***dst, PGbyte ***src, int ZN, int YN, int XN);
void pgBcubeToFcube(PGfloat ***dst, PGbyte ***src, int ZN, int YN, int XN);
void pgFcubeGradient(PGfloat3d ***g, PGfloat ***I, int ZN, int YN, int XN);
void pgFcubeCopy(PGfloat ***dst, PGfloat ***src, int ZN, int YN, int XN);
void pgFcubeAddFcube(PGfloat ***sum, PGfloat ***I1, PGfloat ***I2, 
		     int ZN, int YN, int XN);
void pgFcubeSelfAddFcube(PGfloat ***I, PGfloat ***add, int ZN, int YN, int XN);
void pgFcubeSubFcube(PGfloat ***sum, PGfloat ***I1, PGfloat ***I2, 
		     int ZN, int YN, int XN);
void pgFcubeSelfSubFcube(PGfloat ***I, PGfloat ***sub, int ZN, int YN, int XN);
void pgFcubeSelfDivScalarF(PGfloat ***I, PGfloat value, 
			   int ZN, int YN, int XN);
void pgF3dcubeMag(PGfloat ***mag, PGfloat3d ***v, int ZN, int YN, int XN);
void pgF3dcubeCopy(PGfloat3d ***dst, PGfloat3d ***src, int ZN, int YN, int XN);

/* Endian adaptive functions */
int pgEndianTest();
void pgEndianConvert(void* data, int type, int nItem);
#endif