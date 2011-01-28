/*
 * Exact Medial Axis with Euclidean Distance in 3D
 *
 * Copyright (C) Eric Remy and Edouard Thiel - Oct 2005
 *
 * http://www.lif-sud.univ-mrs.fr/~thiel/IVC2004/LutEucli3D.c
 *
 * This program is free software under the terms of the 
 * GNU Lesser General Public License (LGPL) version 2.1.
 *
 * To compile:  cc -O2 LutEucli3D.c -o LutEucli3D
 * Usage:       ./LutEucli3D
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 * Types and macros. Modify macros to enlarge masks and look-up tables.
*/

typedef struct {
    int x, y, z, r;
} Weighting;

#define MAXWEIGHTING 2000

typedef struct {
    Weighting vg[MAXWEIGHTING];
    int ng;
} MaskG;

typedef unsigned long ImagePoint;
typedef ImagePoint *Image;

#define MAXLUTENTRY 100000

typedef ImagePoint LookUpTable[MAXWEIGHTING][MAXLUTENTRY];


/*
 * Create a new volume L*L*L and initialize the voxels to 0.
 *
 * Input:  L the side length.
 * Output: the pointer to the volume, or NULL if memory allocation failed.
*/

Image NewImage (int L)
{
    int n = L*L*L;
    Image res;

    printf ("NewImage: allocating %dk ... ", n*sizeof(ImagePoint)/1024);
    res = calloc (n, sizeof(ImagePoint));
    if (res == NULL)
         printf ("malloc error\n");
    else printf ("ok\n");
    return res;
}


/*
 * Add a weighting in the mask M.
 *
 * Input:  M the generator of the mask, (x,y,z,r) the weighting.
 * Output: the number of the weighting in the mask.
*/

int AddWeighting (MaskG *M, int x, int y, int z, int r)
{
    int i = M->ng;

    if (i >= MAXWEIGHTING) {
        printf ("AddWeighting: weighting number limited to MAXWEIGHTING = %d\n",
            MAXWEIGHTING);
        return 0;
    }
    if (! (0 <= z && z <= y && y <= x && 0 < x && 0 < r)) {
        printf ("AddWeighting: (x = %d, y = %d, z = %d, r = %d)\n",
            x, y, z, r);
        printf ("  does not respect 0 <= z <= y <= x, 0 < x, 0 < r\n");
        return 0;
    }
    M->vg[i].x = x;
    M->vg[i].y = y;
    M->vg[i].z = z;
    M->vg[i].r = r;
    M->ng++;
    return i;
}


/*
 * Fast Cone Distance Transform Algorithm.
 *
 * Input:  L the side length.
 * Output: CTg the L^3 distance image to the origin.
*/

void CompCTg (int L, Image CTg)
{
    int x, y, z, L2 = L*L;

    CTg[0] = 0;

    for (x = 1; x < L; x++)
    for (y = 0; y <= x; y++)
    for (z = 0; z <= y; z++)
    {
        CTg[ z*L2 + y*L + x ] = x*x + y*y + z*z;
    }
}


/*
 * Fast Distance Transform in G(Z^3).
 *
 * Input:  L the side length,
 *         CTg the distance cone to the origin,
 *         R the radius of the ball in CTg.
 * Output: DTg the distance image to the background of the ball.
 *
 * Algorithm derived from Hirata's SEDT, with improvements (mark unpropagated 
 * distances with -1 in cones, stop loops ASAP., store intersections).
 * See: T. Hirata, "A unified linear-time algorithm fo computing distance maps",
 * in Information Processing Letters 58:129-133, 1996.
*/


double D_Intersec (int u, int gu, int v, int gv)
{
    return (u+v + ((double)gu-gv)/(u-v)) / 2;
}

void CompDTg_Hirata (int L, Image CTg, Image DTg, int R)
{
    int x, xM, y, z, L2=L*L, p, dp, k, propag;
    int Si[L], Sv[L], Sn;  /* index, values, number */
    double Sr[L];          /* intersections */
 
    /* compute bound xM[ and verify that xM < L */
    for (xM = 0; xM < L; xM++) if (CTg[0*L2+0*L+xM] > R) break;
    if (xM >= L) printf ("WARNING xM is not < L\n");
    
    /* First scan: x++, y++, z-- */
    for (x = 0; x <= xM; x++)
    for (y = 0; y <= x; y++)
    {
        k = 0; propag = 0;
        for (z = y; z >= 0; z--)
        {
            p = z*L2+y*L+x;
            if (CTg[p] > R)   /* outside the ball : background */
                propag = 1;
            else if (propag)  /* inside the ball, mark to dist. k*k from bg */
              { k++; DTg[p] = k*k; }
            else DTg[p] = -1; /* inside the ball, no distance propagated */
        }
    }
    
    /* Intermediate scan: x++, z++, y++ */
    for (x = 0; x <= xM; x++)
    for (z = 0; z <= x; z++)
    {
        Sn = 0;  /* empty stacks */
        
        /* Compute stacks indices Si[Sn] and values Sv[Sn] */
        for (y = z; y <= x; y++)
        {
            p = z*L2+y*L+x; dp = DTg[p];
            if (dp < 0) continue;  /* Non propagated value */
            
            /* To speedup algorithm, stop at the second consecutive 0 */
            if (dp == 0 && y > z && DTg[p-L*1] == 0) break;
            
            while (Sn >= 2 && D_Intersec (Si[Sn-1], Sv[Sn-1], y, dp) < Sr[Sn-1])
                Sn--;  /* pop */
                  
            /* push */    
            Si[Sn] = y; Sv[Sn] = dp; 
            if (Sn >= 1) Sr[Sn] = D_Intersec (Si[Sn-1], Sv[Sn-1], Si[Sn], Sv[Sn]);
            Sn++;
        }
        
        if (Sn == 0) continue;  /* Empty stack */
        
        /* Compute new DTg values using stacks */
        for (y = x; y >= z; y--)
        {
            p = z*L2+y*L+x;
            if (DTg[p] == 0) continue;
            
            while (Sn >= 2 && y < Sr[Sn-1]) Sn--;  /* pop */
                  
            DTg[p] = (y-Si[Sn-1])*(y-Si[Sn-1]) + Sv[Sn-1];
        } 
    }
    
    /* Final scan: y++, z++, x++ */
    for (y = 0; y <= xM; y++)
    for (z = 0; z <= y; z++)
    {
        Sn = 0;  /* empty stacks */
        
        /* Compute stacks indices Si[Sn] and values Sv[Sn] */
        for (x = y; x <= xM; x++)
        {
            p = z*L2+y*L+x; dp = DTg[p];
            if (dp < 0) continue;  /* Non propagated value */
            
            /* To speedup algorithm, stop at the second consecutive 0 */
            if (dp == 0 && x > y && DTg[p-1] == 0) break;
            
            while (Sn >= 2 && D_Intersec (Si[Sn-1], Sv[Sn-1], x, dp) < Sr[Sn-1])
                Sn--;  /* pop */
                  
            /* push */    
            Si[Sn] = x; Sv[Sn] = dp; 
            if (Sn >= 1) Sr[Sn] = D_Intersec (Si[Sn-1], Sv[Sn-1], Si[Sn], Sv[Sn]);
            Sn++;
        }
        
        if (Sn == 0) continue;  /* Empty stack */
        
        /* Compute new DTg values using stacks */
        for (x = xM; x >= y; x--)
        {
            p = z*L2+y*L+x;
            if (DTg[p] == 0) continue;
            
            while (Sn >= 2 && x < Sr[Sn-1]) Sn--;  /* pop */
                  
            DTg[p] = (x-Si[Sn-1])*(x-Si[Sn-1]) + Sv[Sn-1];
        } 
    }
}


/* Miscellaneous */

int CompGCD (int a, int b)
{
    if (a*b == 0) return a+b;
    if (a <= b) return CompGCD (a, b % a);
    return CompGCD (a % b, b);
}

int IsVisible (int x, int y, int z)
{
    return CompGCD(x, CompGCD(y,z)) == 1;
}


/*
 * Lut column Computation Algorithm.
 *
 * Input:  CTg the distance cone to origin, L the side length, vg the
 *         weighting and ivg its number, Rmax the greatest radius
 *         to be verified in Lut.
 * Output: The Lut column Lut[ivg] is filled with the correct values.
*/

void CompLutCol (Image CTg, int L, Weighting *vg, int ivg,
                    int Rmax, LookUpTable Lut)
{
    int x, y, z, r, r1, r2, ra, rb, L2 = L*L;

    /* Initializes Lut[ivg] to 0 */
    for (r = 0; r <= Rmax; r++)
        Lut[ivg][r] = 0;

    for (x = 0; x < L - vg->x; x++)
    for (y = 0; y <= x; y++)
    for (z = 0; z <= y; z++)
    {
        /* Radius of the ball where p1 is located */
        r1 = CTg[ z*L2 + y*L + x ] +1;

        /* Same for p2 */
        r2 = CTg[ (z+vg->z)*L2 + (y+vg->y)*L + x+vg->x ] +1;

        if (r1 <= Rmax && r2 > Lut[ivg][r1]) Lut[ivg][r1] = r2;
    }

    rb = 0;
    for (ra = 0; ra <= Rmax; ra++)
    {
        if (Lut[ivg][ra] > rb)
             rb = Lut[ivg][ra];
        else Lut[ivg][ra] = rb;
    }
}


/*
 * Fast extraction of MA points from Bd inter G(Z^3).
 *
 * Input:  x,y,z the point to test, MgL the generator of the Lut mask,
 *         Lut the look-up table, DTg the distance transform of the section of
 *         the ball, L the side length.
 * Output: returns 1 if point x,y,z is detected as MA in the DTg.
*/

int IsMAg (int x, int y, int z, MaskG *MgL, LookUpTable Lut, Image DTg, int L)
{
    int xx, yy, zz, val, i, L2 = L*L;

    val = DTg[ z*L2 + y*L + x ];

    for (i = 0; i < MgL->ng; i++)
    {
        xx = x - MgL->vg[i].x;
        yy = y - MgL->vg[i].y;
        zz = z - MgL->vg[i].z;

        if (0 <= zz && zz <= yy && yy <= xx)
        {
            if ( DTg[ zz*L2 + yy*L + xx ] >= Lut[i][val] )
                return 0;
        }
    }

    return 1;
}


/*
 * Full Lut Computation Algorithm with determination of MgLut.
 *
 * Input:  CTg and DTg two images, L the side length, MgL the generator of the
 *         Lut mask, Lut the look-up table, Rknown the last verified radius, 
 *         Rtarget the maximum radius to be verified.
 * Output: MgL and Lut are filled with the correct values.
 *
 * To compute MgL and Lut from beginning to the maximum radius allowed with L:
 *    MgL.ng = 0;
 *    Rknown = 0;
 *    Rtarget = GreatestRadius (L);
 *    CompLutMask (CTg, DTg, L, &MgL, Lut, Rknown, Rtarget);
 *    Rknown = Rtarget;
*/

void CompLutMask (Image CTg, Image DTg, int L,
                  MaskG *MgL, LookUpTable Lut,
                  int Rknown, int Rtarget)
{
    int x, y, z, R, i, val, L2 = L*L, Possible[MAXLUTENTRY];

    CompCTg (L, CTg);

    /* Init DTg to 0 */ 
    DTg[0] = 0;
    for (x = 1; x < L; x++)
    for (y = 0; y <= x; y++)
    for (z = 0; z <= y; z++)
    {
        DTg[ z*L2 + y*L + x ] = 0;
    }
    
    /* Mark possible values in Possible[] */
    for (i = 1; i < MAXLUTENTRY; i++)
        Possible[i] = 0;

    for (x = 1; x < L; x++)
    for (y = 0; y <= x; y++)
    for (z = 0; z <= y; z++)
    {
        val = CTg[ z*L2 + y*L + x ];
        if (val < MAXLUTENTRY) Possible[val] = 1;
    }

    /* Compute Lut columns for current Lut mask */
    for (i = 0; i < MgL->ng; i++)
        CompLutCol (CTg, L, MgL->vg + i, i, Rtarget, Lut);

    for (R = Rknown+1; R <= Rtarget; R++)
    if (Possible[R])
    {
        printf ("R = %d / %d\r", R, Rtarget); fflush (stdout);

        /* Here we can avoid to init DTg because the ball grows with R. */
        
        CompDTg_Hirata (L, CTg, DTg, R);

        for (x = 1; x < L; x++)
        {
            if (DTg[0*L2 + 0*L + x] == 0) break;  /* Outside the ball */
            
            for (y = 0; y <= x; y++)
            {
                if (DTg[0*L2 + y*L + x] == 0) break; /* Outside the ball */
                
                for (z = 0; z <= y; z++)
                {
                    if (DTg[z*L2 + y*L + x] == 0) break; /* Outside the ball */
                
                    if ( IsMAg (x, y, z, MgL, Lut, DTg, L) )
                    {
                        /* Add a new weighting to MgL */
                        i = AddWeighting (MgL, x, y, z, R);
        
                        printf ("i=%3d  (x,y,z)= ( %2d, %2d, %2d)  added for R= %7d   %s\n",
                            i+1, x, y, z, R,
                            IsVisible(x,y,z)?"visible":"** NON VISIBLE");
        
                        /* New column in Lut */
                        CompLutCol (CTg, L, MgL->vg + i, i, Rtarget, Lut);
        
                        if (IsMAg (x, y, z, MgL, Lut, DTg, L))
                        {
                            printf ("\nCompLutMask: ERROR for R = %d\n", R);
                            return;
                        }
                    }
                }
            }
        }
    }
}


/*
 * Compute the greatest verifiable radius of balls.
 *
 * Input:  L the side length.
 * Output: returns the greatest verifiable radius in the image.
*/

int GreatestRadius (int L)
{
    int res = (L-1)*(L-1)-1;
    
    if (res >= MAXLUTENTRY) {
        printf ("GreatestRadius: maximum radius limited to MAXLUTENTRY = %d\n",
            MAXLUTENTRY);
        res = MAXLUTENTRY-1;
    }
    return res;
}


/*
 * Print a mask M in a file.
 *
 * Input:  f a file, M the generator of a mask.
 * Output: to file f.
*/

void PrintMask (FILE *f, MaskG *M)
{
    int i;

    fprintf (f, "# Computed Lut Mask:\n");
    fprintf (f, "# i  (  x,  y,  z)        R\n");
    for (i = 0; i < M->ng; i++)
        fprintf (f, "%3d  ( %2d, %2d, %2d)  %7d\n",
            i+1, M->vg[i].x, M->vg[i].y, M->vg[i].z, M->vg[i].r);
    fprintf (f, "\n");
}


/*
 * Print the LUT in a file.
 *
 * Input:  f a file, MgL the lut mask, Lut the look-up table, CTg the distance
 *         cone to origin, L the side length, Rknown the maximum verified radius.
 * Output: print the lut colums to file f.
*/

void PrintLut (FILE *f, MaskG *MgL, LookUpTable Lut, Image CTg,
               int L, int Rknown)
{
    int x, y, z, i, j, val, Possible[MAXLUTENTRY], L2 = L*L;

    /* Mark possible values in Possible[] */
    for (i = 1; i < MAXLUTENTRY; i++)
        Possible[i] = 0;

    for (x = 1; x < L; x++)
    for (y = 0; y <= x; y++)
    for (z = 0; z <= y; z++)
    {
        val = CTg[ z*L2 + y*L + x ];
        if (val < MAXLUTENTRY) Possible[val] = 1;
    }

    fprintf (f, "# Look-Up Table: (L = %d, R <= %d)\n", L, Rknown);

    for (j = 1; j < Rknown; j++)
    if (Possible[j])
    {
        fprintf (f, "Lut[][%5d ] = ", j);

        for (i = 0; i < MgL->ng; i++)
            fprintf (f, "%7ld ", Lut[i][j]);

        fprintf (f, "\n");
    }
    fprintf (f, "\n");
}


void AffiUsage ()
{
    fprintf (stderr, "USAGE: LutEucli3D [-m] [-lut] L\n");
    fprintf (stderr, "  -m    save mask in file lut3D-euclidean.txt\n");
    fprintf (stderr, "  -lut  save LUT columns in file lut3D-euclidean.txt; might be HUGE!!\n");
    fprintf (stderr, "  L     image side length in pixels, > 0 \n");
    fprintf (stderr, "\n");
    fprintf (stderr, "Example: ./LutEucli3D 30\n");
}


/*
 * Main program
*/

int main (int argc, char **argv)
{

    MaskG MgL;
    LookUpTable *Lut;
    Image CTg, DTg;
    int L = 0, Rknown, Rtarget, opt_m = 0, opt_lut = 0;
    char name[256];
    FILE *f;

    printf ("LutEucli3D - (C) Eric Remy and Edouard Thiel - Oct 2005\n\n");
    
    /* Parse arguments */
    if (argc == 1) { AffiUsage(); exit (1); }
    while (argc > 1) {
        if (strcmp (argv[1], "-m") == 0) 
          { opt_m = 1; argc--; argv++; }
        else if (strcmp (argv[1], "-lut") == 0) 
          { opt_lut = 1; argc--; argv++; }
        else { L = atoi(argv[1]); break; }
    }
    if (argc != 2 || L <= 0) { AffiUsage(); exit (1); }
    
    printf ("Lut: malloc of %d Kb\n", sizeof(LookUpTable)/1024);
    Lut = malloc (sizeof(LookUpTable));
    if (Lut == NULL) {
        printf ("malloc error\n");
    }

    CTg = NewImage (L);
    DTg = NewImage (L);
    
    /*
     * Computation of Lut and MgL
    */
    printf ("\nComputed Lut Mask:\n");
    MgL.ng = 0;
    Rknown = 0;
    Rtarget = GreatestRadius (L);
    CompLutMask (CTg, DTg, L, &MgL, *Lut, Rknown, Rtarget);
    Rknown = Rtarget;

    /*
     * Write results into a file
    */
    if (opt_m || opt_lut) {
        sprintf (name, "lut3D-euclidean.txt");
        printf ("\n\nSaving results in %s ... ", name);
        f = fopen (name, "w");
        if (f == NULL) printf ("error while opening file\n\n");
        else {
            if (opt_m) PrintMask (f, &MgL);
            if (opt_lut) PrintLut (f, &MgL, *Lut, CTg, L, Rknown); 
            printf ("ok\n\n");
            fclose (f);
        }
    }
    
    free (CTg); free (DTg);
    exit (0);
}
