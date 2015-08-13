#ifndef __MD_DEFS_H__
#define __MD_DEFS_H__

#define LINE 0
#define HELIX 1
#define SAW 2

#define OSEEN 1
#define ROTNE 2
#define NOHI 3

#define I(X,Y) ((X==Y) ? 1.0 : 0.0)
#define tri(X,Y,N) X*N - X*(X-1)/2 + Y - X 
#define squ(i,j,N) i*N + j
#define rank4(i,j,n,m,N) (N*DIM*(DIM*i + n) + (DIM*j + m)) 
#define FOR_ALL_K for(int k=0;k<DIM;k++)


// Constants
#define q0 1.12246204831
#define q0_2 1.25992104989
#define dt 0.00004
#define zeta 1
#define kT 1.0
#define TWOPI 6.28318530718

#define DIM 3

#define VSAdd( r1, r2, s, r3 ) \
  FOR_ALL_K r1[k] = r2[k] + s*r3[k]
#define VAdd( r1, r2, r3 ) \
  VSAdd( r1, r2, 1, r3 )
#define VSub( r1, r2, r3 ) \
  VSAdd( r1, r2, -1, r3 )
#if DIM == 3
  #define VDot( r1, r2 ) \
    (r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2])
#elif DIM == 2
  #define VDot( r1, r2 ) \
    (r1[0]*r2[0] + r1[1]*r2[1])
#endif

#define VSqr( r ) VDot( r, r )
#define VMod( r ) sqrt( VSqr(r) )
#define VZero( r ) \
  FOR_ALL_K r[k] = 0.0
#define VPrint( r ) \
  FOR_ALL_K printf("%lf\t", r[k]); \
  printf("\n") 

;


typedef struct{

  unsigned int numAtoms, NOISE, FORCE, HYDRO; // Population, thermal on/off, type of force, hydro on/off
  unsigned int firstAtomID; // 
  double spring_constant, radius; // used if NOISE is hookean


} Polymers;

/*
typedef struct{

  double *Q, *V, *A; // Position, velocity, ac.
  double mass, size, zeta; 
  int POLID; // Identifies which polymer this one is in for bonding

} Atom;


typedef struct{

  Atom* Afirst;
  int numAtoms;
  int SPID;

} Polymer;
*/
#endif
