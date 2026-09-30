#include <math.h>
#include "toolbox.h"

inline double dotProduct(double a[3], double b[3]) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline void crossProduct(double a[3], double b[3], double c[3]) {
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = a[2]*b[0]-a[0]*b[2];
  c[2] = a[0]*b[1]-a[1]*b[0];
}

void unitVector(double r[3], double n[3]) {
  double l = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];

  /* If 'r' is a null vector, return a unit vector pointing along z-axis */
  if ( l <= 0.00 ) {
    n[0] = 0.00;
    n[1] = 0.00;
    n[2] = 1.00;
  }
  else {
    l = sqrt( l );

    n[0] = r[0] / l;
    n[1] = r[1] / l;
    n[2] = r[2] / l;
  }
}

inline double magnitude(double a[3]) {
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

void rotx(double vin[3], double p, double vout[3]) {
  double cp, sp, v[3];

  v[0] = vin[0];
  v[1] = vin[1];
  v[2] = vin[2];

  cp = cos(p);
  sp = sin(p);

  vout[0] = v[0];
  vout[1] = cp*v[1] + sp*v[2];
  vout[2] = -sp*v[1] + cp*v[2];
}

void roty(double vin[3], double p, double vout[3]) {
  double cp, sp, v[3];

  v[0] = vin[0];
  v[1] = vin[1];
  v[2] = vin[2];

  cp = cos(p);
  sp = sin(p);

  vout[0] = cp*v[0] - sp*v[2];
  vout[1] = v[1];
  vout[2] = sp*v[0] + cp*v[2];
}

void rotz(double vin[3], double p, double vout[3]) {
  double cp, sp, v[3];

  v[0] = vin[0];
  v[1] = vin[1];
  v[2] = vin[2];

  cp = cos(p);
  sp = sin(p);

  vout[0] = cp*v[0] + sp*v[1];
  vout[1] = -sp*v[0] + cp*v[1];
  vout[2] = v[2];
}

inline void applyRotation(double R[3][3], double xo[3], double xp[3]) {
  double x=xo[0];
  double y=xo[1];
  double z=xo[2];

  xp[0] = R[0][0]*x + R[0][1]*y + R[0][2]*z;
  xp[1] = R[1][0]*x + R[1][1]*y + R[1][2]*z;
  xp[2] = R[2][0]*x + R[2][1]*y + R[2][2]*z;
}

void sdpRotationMatrix(double n[3], double M[3][3]) {
  if ( n[2]*n[2] < 1.00 ) {
    double den = sqrt( 1.00 - n[2]*n[2] );

    /* bottom row is just Shower-Detector Plane Normal vector */
    M[2][0] = n[0];
    M[2][1] = n[1];
    M[2][2] = n[2];

    /*
     * x' axis is a vector pointing along the ground, perpendicular to both 
     *   the z & z' axes.  That is, x' = z' x z / | z' x z |.
     */
    M[0][0] = n[1] / den;
    M[0][1] = -n[0] / den;
    M[0][2] = 0.00;
    
    /* y' axis is simply z' x x'. */
    M[1][0] = M[2][1]*M[0][2]-M[2][2]*M[0][1];
    M[1][1] = M[2][2]*M[0][0]-M[2][0]*M[0][2];
    M[1][2] = M[2][0]*M[0][1]-M[2][1]*M[0][0];
  }
  else {
    /* return the identity matrix */
    M[0][0] = M[1][1] = M[2][2] = 1.00;
    M[0][1] = M[0][2] = M[1][2] = 0.00;
    M[1][0] = M[2][0] = M[2][1] = 0.00;
  }
}

void matrixMultiply(double a[3][3], double b[3][3], double c[3][3]) {
  c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
  c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
  c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];

  c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
  c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
  c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];

  c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
  c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
  c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
}

void matrixInverse(double a[3][3], double b[3][3]) {
  int i, j;
  double x[3][3];
  double det;

  for ( i=0; i<3; i++ )
    for ( j=0; j<3; j++ )
      x[i][j] = a[i][j];

  det = x[0][0]*x[1][1]*x[2][2] - x[0][0]*x[1][2]*x[2][1] + 
    x[0][1]*x[1][2]*x[2][0] - x[0][1]*x[2][2]*x[1][0] + 
    x[0][2]*x[1][0]*x[2][1] - x[0][2]*x[2][0]*x[1][1];

  b[0][0] = ( x[1][1]*x[2][2]-x[1][2]*x[2][1] ) / det;
  b[0][1] = ( x[2][1]*x[0][2]-x[2][2]*x[0][1] ) / det;
  b[0][2] = ( x[0][1]*x[1][2]-x[0][2]*x[1][1] ) / det;

  b[1][0] = ( x[1][2]*x[2][0]-x[1][0]*x[2][2] ) / det;
  b[1][1] = ( x[2][2]*x[0][0]-x[2][0]*x[0][2] ) / det;
  b[1][2] = ( x[0][2]*x[1][0]-x[0][0]*x[1][2] ) / det;

  b[2][0] = ( x[1][0]*x[2][1]-x[1][1]*x[2][0] ) / det;
  b[2][1] = ( x[2][0]*x[0][1]-x[2][1]*x[0][0] ) / det;
  b[2][2] = ( x[0][0]*x[1][1]-x[0][1]*x[1][0] ) / det;
}

/*
 *  Generates a rotation matrix that tranforms a coordinate system with z-axis 
 *    pointing along (0,0,1) to one with z'=(z[0],z[1],z[2]), and x' pointing 
 *    in the direction of z x z' (z cross z'), i.e. along the x-y plane and 
 *    perpendicular to both z and z'.  This function is useful for finding 
 *    the inter-site transformations, site-to-SDP transformations, and 
 *    site-to-telescope transformations.  Since this matrix is always unitary, 
 *    inverse matrix is just itself with swapped row & column indices.  If 
 *    z' points straight up, the identity matrix (I) is returned.  If 
 *    z' points straight down (0,0,-1), then -I (minus identity) matrix is 
 *    returned.
 */
void zmatrix(double z[3], double m[3][3]) {
  double l = sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);

  /* ensure that bottom row has unit length */
  m[2][0] = z[0] / l;
  m[2][1] = z[1] / l;
  m[2][2] = z[2] / l;

  if ( fabs(m[2][2]) < 1.00 ) {  /* avoid floating point errors */
    m[0][0] = -m[2][1] / sqrt(1.00 - m[2][2]*m[2][2]);
    m[0][1] =  m[2][0] / sqrt(1.00 - m[2][2]*m[2][2]);
    m[0][2] =  0.00;

    m[1][0] = m[2][1]*m[0][2] - m[2][2]*m[0][1];
    m[1][1] = m[2][2]*m[0][0] - m[2][0]*m[0][2];
    m[1][2] = m[2][0]*m[0][1] - m[2][1]*m[0][0];
  }
  else {  /* return the identity matrix */
    m[0][0] = 1.00;
    m[0][1] = m[0][2] = 0.00;

    m[1][1] = 1.00;
    m[1][0] = m[1][2] = 0.00;

    m[2][2] = 1.00;
    m[2][0] = m[2][1] = 0.00;
  }  
}

inline void matrixTranspose(double M[3][3], double Mt[3][3]) {
  int i, j;
  double Mp[3][3];

  /* copy the input matrix, in case M & Mt point to the same matrix */
  for ( i=0; i<3; i++ )
    for ( j=0; j<3; j++ )
      Mp[i][j] = M[i][j];

  /* fill the contents of Mt with the input matrix, with indices reversed */
  for ( i=0; i<3; i++ )
    for ( j=0; j<3; j++ )
      Mt[i][j] = Mp[j][i];

}
