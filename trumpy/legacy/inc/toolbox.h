#ifndef _TOOLBOX_H_
#define _TOOLBOX_H_

#include <stdio.h>
#include "constants.h"
#include "control.h"

#include "convcoord.h"

#include "convtime.h"

#include "stools.h"

/*
 * Math functions that occur frequently, but can be done in a few lines.
 */
/* Replaced with the greater value of X or Y. */
#define max(X,Y)       ((X>Y)?X:Y)

/* Replaced with the lesser value of X or Y. */
#define min(X,Y)       ((X>Y)?Y:X)

/* copies the contents of A to B */
/* #define cpyvec(A,B)    {B[0]=A[0];B[1]=A[1];B[2]=A[2];} */
#define cpyvec(A,B)    {int _i;for (_i=0;_i<3;_i++) B[_i]=A[_i];}

/* Sets vector C equal to A + B */
/* #define addvec(A,B,C)  {C[0]=A[0]+B[0];C[1]=A[1]+B[1];C[2]=A[2]+B[2];} */
#define addvec(A,B,C)  {int _i;for (_i=0;_i<3;_i++) C[_i]=A[_i]+B[_i];}

/* Sets vector C equal to A - B */
/* #define addvec(A,B,C)  {C[0]=A[0]-B[0];C[1]=A[1]-B[1];C[2]=A[2]-B[2];} */
#define subvec(A,B,C)  {int _i;for (_i=0;_i<3;_i++) C[_i]=A[_i]-B[_i];}

/* Replaced with the dot product of vectors A and B */
#define dotprod(A,B)   (A[0]*B[0]+A[1]*B[1]+A[2]*B[2])

/* Replaced with the magnitude of Vector A */
#define magvec(A)      sqrt(dotprod(A,A))

/* Sets vector C equal to the cross product A x B */
#define crsprod(A,B,C) {C[0]=A[1]*B[2]-A[2]*B[1];\
                        C[1]=A[2]*B[0]-A[0]*B[2];\
                        C[2]=A[0]*B[1]-A[1]*B[0];}

/* Replaced with X^2 */
#define sqr(X)         ((X)*(X))

/* Replaced with X^4 */
#define pow4(X)        ((X)*(X)*(X)*(X))

/* Age formula */
#define ageS(X,XMAX)   (3./(1.+2.*(XMAX)/(X)))

/* for whatever reason, this function is not included in math.h */
#define round(X)       ((double)((int)((X)+0.50)))

/* these have been retired (conflict with pthread(?) - SS 10.4.11 */
/* #define clone(V,A,N)  {int _i;for (_i=0;_i<N;_i++) A[_i]=V;} */
/* #define zeros(A,N)    clone(0.0,A,N) */
/* #define ones(A,N)     clone(1.0,A,N) */

/* in "list-items.c" */
#ifdef __cplusplus
extern "C"{
#endif
int addListItem (int val, int *n, int max, int a[]);

/* in "indexx.c" */
void indexx(int n, double *arr, int *indx);

/* in "jday_to_GMST.c" */
double jday_to_GMST (double jday);

/* in "linefit.c" */
void linefit (int num, double x[], double y[], double ey[],
              double *m, double *em, double *b, double *eb, double *chi2);

/* in "radec_to_zenazm.c" */
void radec_to_zenazm (double jday, double lat, double lon,
                      double ra, double dec, double *zen, double *azm);
void radec_to_zenazm_ecc (double jday, double lat, double lon, double alt,
			  double ra, double dec, double *zen, double *azm);

/* in "zenazm_to_radec.c" */
void zenazm_to_radec (double jday, double lat, double lon,
                      double zen, double azm, double *ra, double *dec);
void zenazm_to_radec_ecc (double jday, double lat, double lon, double alt,
			  double zen, double azm, double *ra, double *dec);
/* in "range.c" */
double range (double X, double Y);

/* in "significance.c" */
double significance(int mean, int var, int *wf, int ns);

/* in "vectormath.c" */
double dotProduct(double a[3], double b[3]);
void crossProduct(double a[3], double b[3], double c[3]);
void unitVector(double r[3], double n[3]);
double magnitude(double a[3]);
void rotx(double vin[3], double p, double vout[3]);
void roty(double vin[3], double p, double vout[3]);
void rotz(double vin[3], double p, double vout[3]);
void applyRotation(double R[3][3], double xo[3], double xp[3]);
void sdpRotationMatrix(double n[3], double M[3][3]);
void matrixInverse(double a[3][3], double b[3][3]);
void matrixTranspose(double M[3][3], double Mt[3][3]);
void matrixMultiply(double a[3][3], double b[3][3], double c[3][3]);
void zmatrix(double z[3], double m[3][3]);

/* in qsimp.c */
double qsimp(double (*func)(double), double a, double b, double eps);

/* in trapzd.c */
double trapzd(double (*func)(double), double a, double b, int n);

/* in nrutil.c */
void nrerror(char error_text[]);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long newrl, long newcl);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(float *v, long nl);
void free_ivector(int *v, long nl);
void free_cvector(unsigned char *v, long nl);
void free_lvector(unsigned long *v, long nl);
void free_dvector(double *v, long nl);
void free_matrix(float **m, long nrl, long ncl);
void free_dmatrix(double **m, long nrl, long ncl);
void free_imatrix(int **m, long nrl, long ncl);
void free_submatrix(float **b, long nrl);
void free_convert_matrix(float **b, long nrl);
void free_f3tensor(float ***t, long nrl, long ncl, long ndl);

/* in 'binsearch.c' */
int binsearch(double x, int j, double *xi);

/* in 'readline.c' */
int readline(FILE *fp, char *line, char args[][MAX_STRLEN]);

/* in 'fadccl.c' */
double fadccl(int mean, int vari, int nbin, short *wf);

/* in 'disttoalt.c' */
double getDistanceToAltitude_ECC(double rr_ec[3], double uv_ec[3], double z, double rz_ec[3]);
double getDistanceToAltitude(double rr[3], double uv[3], double z, double rz[3]);

/* in 'heapsort.c' */
void heapsort(int n, double x[]);
void heapsortp(int n, double x[], double *p[]);

#ifdef __cplusplus
}
#endif

#endif

