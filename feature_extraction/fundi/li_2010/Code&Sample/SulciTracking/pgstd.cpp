
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pgstd.h"

/* Prints error message */
PGvoid pgError(char error_text[])
{
	//PGvoid exit();

	fprintf(stderr,"PG Tools run-time error:\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"Now exiting to system.\n");
	exit(0);
}



/********************************************************/
/* Vector allocation and freeing functions              */
/********************************************************/

/* Allocates a vector of unsigned chars with range [nl..nh] */
PGbyte *pgBvector(int nl,
		  int nh)
{
        PGbyte *v;

        v=(PGbyte *)malloc((unsigned) (nh-nl+1)*sizeof(PGbyte));
        if (!v) pgError("allocation failure in pgBvector()");
        return v-nl;
}

/* Frees a vector allocated by pgBvector */
PGvoid pgFreeBvector(PGbyte *v,
		     int nl,
		     int nh)
{
        free((char*) (v+nl));
}

/* Allocates a vector of unsigned shorts with range [nl..nh] */
PGushort *pgUsvector(int nl,
		     int nh)
{
        PGushort *v;

        v=(PGushort *)malloc((unsigned) (nh-nl+1)*sizeof(PGushort));
        if (!v) pgError("allocation failure in pgUsvector()");
        return v-nl;
}



/* Frees a vector allocated by pgUsvector */
PGvoid pgFreeUsvector(PGushort *v,
		       int nl,
		       int nh)
{
        free((char*) (v+nl));
}
/* Allocates a vector of unsigned shorts with range [nl..nh] */
PGshort *pgSvector(int nl,
		     int nh)
{
        PGshort *v;

        v=(PGshort *)malloc((unsigned) (nh-nl+1)*sizeof(PGshort));
        if (!v) pgError("allocation failure in pgSvector()");
        return v-nl;
}



/* Frees a vector allocated by pgUsvector */
PGvoid pgFreeSvector(PGshort *v,
		       int nl,
		       int nh)
{
        free((char*) (v+nl));
}


/* Allocates a vector of unsigned ints with range [nl..nh] */
PGuint *pgUivector(int nl,
		   int nh)
{
        PGuint *v;

        v=(PGuint *)malloc((unsigned) (nh-nl+1)*sizeof(PGuint));
        if (!v) pgError("allocation failure in pgUivector()");
        return v-nl;
}

/* Frees a vector allocated by pgUivector */
PGvoid pgFreeUivector(PGuint *v,
		       int nl,
		       int nh)
{
        free((char*) (v+nl));
}

/* Allocates a vector of ints with range [nl..nh] */
PGint *pgIvector(int nl, int nh)
{
	PGint *v;

	v=(PGint *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) pgError("allocation failure in pgIvector()");
	return v-nl;
}



/* Frees a vector allocated by pgIvector */
PGvoid pgFreeIvector(int *v,int nl,int nh)
{
	free((char*) (v+nl));
}





/* Allocates a vector of floats with range [nl..nh] */
PGfloat *pgFvector(int nl, int nh)
{
	PGfloat *v;

	v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
	if (!v) pgError("allocation failure in pgFvector()");
	return v-nl;
}


/* Frees a vector allocated by fvector */
PGvoid pgFreeFvector(PGfloat *v, int nl,int nh)
{
	free((char*) (v+nl));
}






/* Allocates a vector of doubles with range [nl..nh] */
PGdouble *pgDvector(int nl, int nh)
{
	PGdouble *v;

	v=(PGdouble *)malloc((unsigned) (nh-nl+1)*sizeof(PGdouble));
	if (!v) pgError("allocation failure in pgDvector()");
	return v-nl;
}

/* Frees a vector allocated by dvector */
PGvoid pgFreeDvector(PGdouble *v, int nl, int nh)
{
	free((char*) (v+nl));
}

/* Allocates a vector of byte, 2-D vectors with range [nl..nh] */
PGbyte2d *pgB2dvector(int nl, int nh)
{
	PGbyte2d *v;

	v=(PGbyte2d *)malloc((unsigned) (nh-nl+1)*sizeof(PGbyte2d));
	if (!v) pgError("allocation failure in pgB2dvector()");
	return v-nl;
}


/* Frees a vector allocated by pgB2dvector */
PGvoid pgFreeB2dvector(PGbyte2d *v,int nl, int nh)
{
	free((char*) (v+nl));
}

/* Allocates a vector of int, 2-D vectors with range [nl..nh] */
PGint2d *pgI2dvector(int nl, int nh)
{
	PGint2d *v;

	v=(PGint2d *)malloc((unsigned) (nh-nl+1)*sizeof(PGint2d));
	if (!v) pgError("allocation failure in pgI2dvector()");
	return v-nl;
}


/* Frees a vector allocated by pgI2dvector */
PGvoid pgFreeI2dvector(PGint2d *v, int nl,int nh)
{
	free((char*) (v+nl));
}

/* Allocates a vector of float, 2-D vectors with range [nl..nh] */
PGfloat2d *pgF2dvector(int nl, int nh)
{
	PGfloat2d *v;

	v=(PGfloat2d *)malloc((unsigned) (nh-nl+1)*sizeof(PGfloat2d));
	if (!v) pgError("allocation failure in pgF2dvector()");
	return v-nl;
}


/* Frees a vector allocated by pgF2dvector */
PGvoid pgFreeF2dvector(PGfloat2d *v, int nl, int nh)
{
	free((char*) (v+nl));
}

/* Allocates a vector of double, 2-D vectors with range [nl..nh] */
PGdouble2d *pgD2dvector(int nl, int nh)
{
	PGdouble2d *v;

	v=(PGdouble2d *)malloc((unsigned) (nh-nl+1)*sizeof(PGdouble2d));
	if (!v) pgError("allocation failure in pgD2dvector()");
	return v-nl;
}


/* Frees a vector allocated by d2dvector */
PGvoid pgFreeD2dvector(PGdouble2d *v,int nl,int nh)
{
	free((char*) (v+nl));
}




/* Allocates a vector of byte, 3-D vectors with range [nl..nh] */
PGbyte3d *pgB3dvector(int nl,int nh)
{
	PGbyte3d *v;

	v=(PGbyte3d *)malloc((unsigned) (nh-nl+1)*sizeof(PGbyte3d));
	if (!v) pgError("allocation failure in pgB3dvector()");
	return v-nl;
}



/* Frees a vector allocated by pgB3dvector */
PGvoid pgFreeB3dvector(PGbyte3d *v, int nl,int nh)
{
	free((char*) (v+nl));
}


/* Allocates a vector of int, 3-D vectors with range [nl..nh] */
PGint3d *pgI3dvector(int nl,int nh)
{
	PGint3d *v;

	v=(PGint3d *)malloc((unsigned) (nh-nl+1)*sizeof(PGint3d));
	if (!v) pgError("allocation failure in pgI3dvector()");
	return v-nl;
}



/* Frees a vector allocated by pgI3dvector */
PGvoid pgFreeI3dvector(PGint3d *v, int nl,int nh)
{
	free((char*) (v+nl));
}




/* Allocates a vector of float, 3-D vectors with range [nl..nh] */
PGfloat3d *pgF3dvector(int nl,int nh)
{
	PGfloat3d *v;

	v=(PGfloat3d *)malloc((unsigned) (nh-nl+1)*sizeof(PGfloat3d));
	if (!v) pgError("allocation failure in pgF3dvector()");
	return v-nl;
}


/* Frees a vector allocated by pgF3dvector */
PGvoid pgFreeF3dvector(PGfloat3d *v,int nl,int nh)
{
	free((char*) (v+nl));
}




/* Allocates a vector of double, 3-D vectors with range [nl..nh] */
PGdouble3d *pgD3dvector(int nl,int nh)
{
	PGdouble3d *v;

	v=(PGdouble3d *)malloc((unsigned) (nh-nl+1)*sizeof(PGdouble3d));
	if (!v) pgError("allocation failure in pgD3dvector()");
	return v-nl;
}



/* Frees a vector allocated by pgD3dvector */
PGvoid pgFreeD3dvector(PGdouble3d *v,int nl,int nh)
{
	free((char*) (v+nl));
}




/********************************************************/
/* Matrix allocation and freeing functions              */
/********************************************************/

/* Allocates a matrix of unsigned chars with range [nrl..nrh][ncl..nch] */
PGbyte **pgBmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGbyte **m;

  bufsize = (nrh-nrl+1)*sizeof(PGbyte*)
            + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGbyte);

  m=(PGbyte **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgBmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGbyte*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGbyte *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGbyte)));
    m[j] -= ncl;
  }

  return m;
}




/* Frees a matrix allocated by bmatrix */
PGvoid pgFreeBmatrix(PGbyte **m,int nrl,int nrh,int ncl,int nch)
{
        free((char*) (m+nrl));
}

/* Allocates a matrix of unsigned shorts with range [nrl..nrh][ncl..nch] */
PGushort **pgUsmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGushort bufsize,bufptr;
  PGushort **m;

  bufsize = (nrh-nrl+1)*sizeof(PGushort*)
            + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGushort);

  m=(PGushort **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgUsmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGushort*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGushort *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGushort)));
    m[j] -= ncl;
  }

  return m;
}



/* frees a matrix allocated by usmatrix */
PGvoid pgFreeUsmatrix(PGushort **m,int nrl,int nrh,int ncl,int nch)
{
        free((char*) (m+nrl));
}

/* Allocates a matrix of signed shorts with range [nrl..nrh][ncl..nch] */
PGshort **pgSmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGshort bufsize,bufptr;
  PGshort **m;

  bufsize = (nrh-nrl+1)*sizeof(PGshort*)
            + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGshort);

  m=(PGshort **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgSmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGshort*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGshort *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGshort)));
    m[j] -= ncl;
  }

  return m;
}



/* frees a matrix allocated by Smatrix */
PGvoid pgFreeSmatrix(PGshort **m,int nrl,int nrh,int ncl,int nch)
{
        free((char*) (m+nrl));
}


/* Allocates a matrix of unsigned ints with range [nrl..nrh][ncl..nch] */
PGuint **pgUimatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGuint **m;

  bufsize = (nrh-nrl+1)*sizeof(PGuint*)
            + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGuint);

  m=(PGuint **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgUimatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGuint*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGuint *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGuint)));
    m[j] -= ncl;
  }

  return m;
}



/* Frees a matrix allocated by uimatrix */
PGvoid pgFreeUimatrix(PGuint **m,int nrl,int nrh,int ncl,int nch)
{
        free((char*) (m+nrl));
}



/* Allocates a matrix of doubles with range [nrl..nrh][ncl..nch] */
PGdouble **pgDmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGdouble **m;

  bufsize = (nrh-nrl+1)*sizeof(PGdouble*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGdouble);

  m=(PGdouble **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgDmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGdouble*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGdouble *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGdouble)));
    m[j] -= ncl;
  }

  return m;
}



/* Frees a matrix allocated by dmatrix */
PGvoid pgFreeDmatrix(PGdouble **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}



/* Allocates a matrix of floats with range [nrl..nrh][ncl..nch] */
PGfloat **pgFmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGfloat **m;

  bufsize = (nrh-nrl+1)*sizeof(PGfloat*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGfloat);

  m=(PGfloat **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgFmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGfloat*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGfloat *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGfloat)));
    m[j] -= ncl;
  }

  return m;
}


/* Frees a matrix allocated by fmatrix */
PGvoid pgFreeFmatrix(PGfloat **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}




/* Allocates a matrix of ints with range [nrl..nrh][ncl..nch] */
PGint **pgImatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGint **m;

  bufsize = (nrh-nrl+1)*sizeof(PGint*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGint);

  m=(PGint **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgImatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGint*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGint *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGint)));
    m[j] -= ncl;
  }

  return m;
}



/* Frees a matrix allocated by pgImatrix */
PGvoid pgFreeImatrix(PGint **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}


/* Allocates a matrix of byte, 2-D vectors with range [nrl..nrh][ncl..nch] */
PGbyte2d **pgB2dmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGbyte2d **m;

  bufsize = (nrh-nrl+1)*sizeof(PGbyte2d*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGbyte2d);

  m=(PGbyte2d **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgB2dmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGbyte2d*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGbyte2d *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGbyte2d)));
    m[j] -= ncl;
  }

  return m;
}

/* Frees a matrix allocated by pgB2dmatrix */
PGvoid pgFreeB2dmatrix(PGbyte2d **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}


/* Allocates a matrix of int, 2-D vectors with range [nrl..nrh][ncl..nch] */
PGint2d **pgI2dmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGint2d **m;

  bufsize = (nrh-nrl+1)*sizeof(PGint2d*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGint2d);

  m=(PGint2d **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgI2dmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGint2d*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGint2d *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGint2d)));
    m[j] -= ncl;
  }

  return m;
}

/* Frees a matrix allocated by pgI2dmatrix */
PGvoid pgFreeI2dmatrix(PGint2d **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}


/* Allocates a matrix of float, 2-D vectors with range [nrl..nrh][ncl..nch] */
PGfloat2d **pgF2dmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGfloat2d **m;

  bufsize = (nrh-nrl+1)*sizeof(PGfloat2d*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGfloat2d);

  m=(PGfloat2d **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgF2dmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGfloat2d*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGfloat2d *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGfloat2d)));
    m[j] -= ncl;
  }

  return m;
}

/* Frees a matrix allocated by pgF2dmatrix */
PGvoid pgFreeF2dmatrix(PGfloat2d **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}




/* Allocates a matrix of double, 2-D vectors with range [nrl..nrh][ncl..nch] */
PGdouble2d **pgD2dmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGdouble2d **m;

  bufsize = (nrh-nrl+1)*sizeof(PGdouble2d*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGdouble2d);

  m=(PGdouble2d **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgD2dmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGdouble2d*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGdouble2d *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGdouble2d)));
    m[j] -= ncl;
  }

  return m;
}

/* Frees a matrix allocated by d2dmatrix */
PGvoid pgFreeD2dmatrix(PGdouble2d **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}

/* Allocates a matrix of byte, 3-D vectors with range [nrl..nrh][ncl..nch] */
PGbyte3d **pgB3dmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGbyte3d **m;

  bufsize = (nrh-nrl+1)*sizeof(PGbyte3d*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGbyte3d);

  m=(PGbyte3d **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgB3dmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGbyte3d*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGbyte3d *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGbyte3d)));
    m[j] -= ncl;
  }

  return m;
}


/* Frees a matrix allocated by pgB3dmatrix */
PGvoid pgFreeB3dmatrix(PGbyte3d **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}


/* Allocates a matrix of int, 3-D vectors with range [nrl..nrh][ncl..nch] */
PGint3d **pgI3dmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGint3d **m;

  bufsize = (nrh-nrl+1)*sizeof(PGint3d*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGint3d);

  m=(PGint3d **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgI3dmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGint3d*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGint3d *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGint3d)));
    m[j] -= ncl;
  }

  return m;
}


/* Frees a matrix allocated by pgI3dmatrix */
PGvoid pgFreeI3dmatrix(PGint3d **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}


/* Allocates a matrix of float, 3-D vectors with range [nrl..nrh][ncl..nch] */
PGfloat3d **pgF3dmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGfloat3d **m;

  bufsize = (nrh-nrl+1)*sizeof(PGfloat3d*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGfloat3d);

  m=(PGfloat3d **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgF3dmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGfloat3d*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGfloat3d *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGfloat3d)));
    m[j] -= ncl;
  }

  return m;
}


/* Frees a matrix allocated by pgF3dmatrix */
PGvoid pgFreeF3dmatrix(PGfloat3d **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}



/* Allocates a matrix of double, 3-D vectors with range [nrl..nrh][ncl..nch] */
PGdouble3d **pgD3dmatrix(int nrl,int nrh,int ncl,int nch)
{
  int j;
  PGuint bufsize,bufptr;
  PGdouble3d **m;

  bufsize = (nrh-nrl+1)*sizeof(PGdouble3d*)
	    + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGdouble3d);

  m=(PGdouble3d **) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgD3dmatrix()");
  m -= nrl;

  bufptr = ((unsigned) (m+nrl)) + (nrh-nrl+1)*sizeof(PGdouble3d*);
  for(j=nrl;j<=nrh;j++) {
    m[j] = ((PGdouble3d *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGdouble3d)));
    m[j] -= ncl;
  }

  return m;
}


/* Frees a matrix allocated by d3dmatrix */
PGvoid pgFreeD3dmatrix(PGdouble3d **m,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+nrl));
}



/********************************************************/
/* Cube allocation and freeing functions                */
/********************************************************/

/* Allocates a cube of unsigned shorts with 
   range [npl..nph][nrl..nrh][ncl..nch] */
PGushort ***pgUscube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGushort ***m;

  bufsize = (nph-npl+1)*sizeof(PGushort**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGushort*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGushort);

  m=(PGushort ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgUscube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGushort**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGushort **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGushort*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGushort*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGushort *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGushort)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGushort)));
      m[i][j] -= ncl;
    }

  return m;
}

/* Frees a cube allocated by uscube */
PGvoid pgFreeUscube(PGushort ***m,int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}

/* Allocates a cube of signed shorts with 
   range [npl..nph][nrl..nrh][ncl..nch] */
PGshort ***pgScube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGshort ***m;

  bufsize = (nph-npl+1)*sizeof(PGshort**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGshort*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGshort);

  m=(PGshort ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgScube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGshort**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGshort **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGshort*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGshort*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGshort *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGshort)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGshort)));
      m[i][j] -= ncl;
    }

  return m;
}

/* Frees a cube allocated by scube */
PGvoid pgFreeScube(PGshort ***m,int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}

/* Allocates a cube of doubles with range [npl..nph][nrl..nrh][ncl..nch] */
PGdouble ***pgDcube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGdouble ***m;

  bufsize = (nph-npl+1)*sizeof(PGdouble**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGdouble*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGdouble);

  m=(PGdouble ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgDcube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGdouble**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGdouble **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGdouble*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGdouble*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGdouble *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGdouble)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGdouble)));
      m[i][j] -= ncl;
    }

  return m;
}

/* Frees a cube allocated by dcube */
PGvoid pgFreeDcube(PGdouble ***m,int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}




/* Allocates a cube of integers with range [npl..nph][nrl..nrh][ncl..nch] */
PGint ***pgIcube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGint ***m;

  bufsize = (nph-npl+1)*sizeof(PGint**)
            + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGint*)
            + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGint);

  m=(PGint ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgIcube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGint**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGint **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGint*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGint*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGint *) (bufptr
                + (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGint)
                + (j-nrl)*(nch-ncl+1)*sizeof(PGint)));
      m[i][j] -= ncl;
    }

  return m;
}


/* Frees a cube allocated by icube */
PGvoid pgFreeIcube(PGint ***m,int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}



/* Allocates a cube of PGbyte2d, 2D vectors with 
   range [npl..nph][nrl..nrh][ncl..nch] */
PGbyte2d ***pgB2dcube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGbyte2d ***m;

  bufsize = (nph-npl+1)*sizeof(PGbyte2d**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGbyte2d*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGbyte2d);

  m=(PGbyte2d ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgB2dcube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGbyte2d**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGbyte2d **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGbyte2d*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGbyte2d*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGbyte2d *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGbyte2d)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGbyte2d)));
      m[i][j] -= ncl;
    }

  return m;
}

/* Frees a cube allocated by pgB2dcube */
PGvoid pgFreeB2dcube(PGbyte2d ***m,int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}

/* Allocates a cube of PGint2d, 2D vectors with 
   range [npl..nph][nrl..nrh][ncl..nch] */
PGint2d ***pgI2dcube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGint2d ***m;

  bufsize = (nph-npl+1)*sizeof(PGint2d**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGint2d*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGint2d);

  m=(PGint2d ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgI2dcube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGint2d**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGint2d **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGint2d*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGint2d*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGint2d *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGint2d)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGint2d)));
      m[i][j] -= ncl;
    }

  return m;
}

/* Frees a cube allocated by pgI2dcube */
PGvoid pgFreeI2dcube(PGint2d ***m,int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}


/* Allocates a cube of PGfloat2d, 2D vectors with 
   range [npl..nph][nrl..nrh][ncl..nch] */
PGfloat2d ***pgF2dcube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGfloat2d ***m;

  bufsize = (nph-npl+1)*sizeof(PGfloat2d**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGfloat2d*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGfloat2d);

  m=(PGfloat2d ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgF2dcube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGfloat2d**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGfloat2d **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGfloat2d*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGfloat2d*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGfloat2d *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGfloat2d)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGfloat2d)));
      m[i][j] -= ncl;
    }

  return m;
}

/* Frees a cube allocated by pgF2dcube */
PGvoid pgFreeF2dcube(PGfloat2d ***m,int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}


/* Allocates a cube of double, 2D vectors with 
   range [npl..nph][nrl..nrh][ncl..nch] */
PGdouble2d ***pgD2dcube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGdouble2d ***m;

  bufsize = (nph-npl+1)*sizeof(PGdouble2d**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGdouble2d*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGdouble2d);

  m=(PGdouble2d ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgD2dcube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGdouble2d**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGdouble2d **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGdouble2d*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGdouble2d*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGdouble2d *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGdouble2d)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGdouble2d)));
      m[i][j] -= ncl;
    }

  return m;
}

/* Frees a cube allocated by pgD2dcube */
PGvoid pgFreeD2dcube(PGdouble2d ***m,int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}


/* Allocates a cube of byte, 3D vectors with 
   range [npl..nph][nrl..nrh][ncl..nch] */
PGbyte3d ***pgB3dcube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGbyte3d ***m;

  bufsize = (nph-npl+1)*sizeof(PGbyte3d**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGbyte3d*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGbyte3d);

  m=(PGbyte3d ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgB3dcube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGbyte3d**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGbyte3d **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGbyte3d*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGbyte3d*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGbyte3d *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGbyte3d)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGbyte3d)));
      m[i][j] -= ncl;
    }

  return m;
}

/* Frees a cube allocated by pgB3dcube */
PGvoid pgFreeB3dcube(PGbyte3d ***m,int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}

/* Allocates a cube of int, 3D vectors with 
   range [npl..nph][nrl..nrh][ncl..nch] */
PGint3d ***pgI3dcube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGint3d ***m;

  bufsize = (nph-npl+1)*sizeof(PGint3d**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGint3d*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGint3d);

  m=(PGint3d ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgI3dcube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGint3d**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGint3d **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGint3d*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGint3d*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGint3d *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGint3d)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGint3d)));
      m[i][j] -= ncl;
    }

  return m;
}

/* Frees a cube allocated by pgI3dcube */
PGvoid pgFreeI3dcube(PGint3d ***m,int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}



/* Allocates a cube of float, 3-D vectors with 
   range [npl..nph][nrl..nrh][ncl..nch] */
PGfloat3d ***pgF3dcube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGfloat3d ***m;

  bufsize = (nph-npl+1)*sizeof(PGfloat3d**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGfloat3d*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGfloat3d);

  m=(PGfloat3d ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgF3dcube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGfloat3d**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGfloat3d **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGfloat3d*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGfloat3d*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGfloat3d *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGfloat3d)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGfloat3d)));
      m[i][j] -= ncl;
    }

  return m;
}




/* Frees a cube allocated by f3dcube */
PGvoid pgFreeF3dcube(PGfloat3d ***m,int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}



/* Allocates a cube of double, 3D vectors with 
   range [npl..nph][nrl..nrh][ncl..nch] */
PGdouble3d ***pgD3dcube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGdouble3d ***m;

  bufsize = (nph-npl+1)*sizeof(PGdouble3d**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGdouble3d*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGdouble3d);

  m=(PGdouble3d ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgD3dcube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGdouble3d**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGdouble3d **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGdouble3d*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGdouble3d*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGdouble3d *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGdouble3d)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGdouble3d)));
      m[i][j] -= ncl;
    }

  return m;
}

/* Frees a cube allocated by d3dcube */
PGvoid pgFreeD3dcube(PGdouble3d ***m,int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}



/* Allocates a cube of floats with range [npl..nph][nrl..nrh][ncl..nch] */
PGfloat ***pgFcube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGfloat ***m;

  bufsize = (nph-npl+1)*sizeof(PGfloat**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGfloat*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGfloat);

  m=(PGfloat ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgFcube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGfloat**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGfloat **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGfloat*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGfloat*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGfloat *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGfloat)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGfloat)));
      m[i][j] -= ncl;
    }

  return m;
}




/* Frees a cube allocated by fcube */
PGvoid pgFreeFcube(PGfloat ***m, int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}




/* Allocates a cube of unsigned chars with 
   range [npl..nph][nrl..nrh][ncl..nch] */
PGbyte ***pgBcube(int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
  int i,j;
  PGuint bufsize,bufptr;
  PGbyte ***m;

  bufsize = (nph-npl+1)*sizeof(PGbyte**)
	    + (nph-npl+1)*(nrh-nrl+1)*sizeof(PGbyte*)
	    + (nph-npl+1)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGbyte);

  m=(PGbyte ***) malloc(bufsize);
  if (!m) pgError("allocation failure 1 in pgBcube()");
  m -= npl;

  bufptr = ((unsigned) (m+npl)) + (nph-npl+1)*sizeof(PGbyte**);
  for(i=npl;i<=nph;i++) {
    m[i] = ((PGbyte **) (bufptr+(i-npl)*(nrh-nrl+1)*sizeof(PGbyte*)));
    m[i] -= nrl;
  }

  bufptr += (unsigned) (nph-npl+1)*(nrh-nrl+1)*sizeof(PGbyte*);
  for(i=npl;i<=nph;i++)
    for(j=nrl;j<=nrh;j++) {
      m[i][j] = ((PGbyte *) (bufptr
		+ (i-npl)*(nrh-nrl+1)*(nch-ncl+1)*sizeof(PGbyte)
		+ (j-nrl)*(nch-ncl+1)*sizeof(PGbyte)));
      m[i][j] -= ncl;
    }

  return m;
}




/* Frees a cube allocated by bcube */
PGvoid pgFreeBcube(PGbyte ***m, int npl,int nph,int nrl,int nrh,int ncl,int nch)
{
	free((char*) (m+npl));
}













