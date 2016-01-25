#ifndef __MD_DEFS_H__
#define __MD_DEFS_H__

// Types of initialising shape
#define LINE 0
#define NORMAL 0
#define HELIX 1
#define SAW 2
#define TRACER 3


// HI tensors ... WARNING: Only use ROTNE
#define NOHI 0
#define OSEEN 1
#define ROTNE2 2      // 2D Rotne-Prager Tensor
#define ROTNE3 3      // 3D RP Tensor
#define ROTNE ROTNE3  // Default will be interpreted as 3D

#define I(X,Y) ((X==Y) ? 1.0 : 0.0)

// Triangular, N*N, and DIM*N*DIM*N matrix indexing
#define tri(X,Y,N) (X*N - X*(X-1)/2 + Y - X) // Not used, but could be useful in future.
#define squ(i,j,N) (i*N + j)
#define rank4(i,j,n,m,N) (N*DIM*(DIM*i + n) + (DIM*j + m))

// I use this a lot
#define FOR_ALL_K for(int k=0;k<DIM;k++)

// Constants
#define MD_q0 1.12246204831
#define MD_q0_2 1.25992104989
#define MD_fene_rr0 2.25 // This is r0^2 where the r0 in the fene potential is 1.5
#define MD_dt 0.000001
#define MD_zeta 1
#define MD_kT 1.0 // This currently doesn't enter into the code anywhere... 
#define MD_TWOPI 6.28318530718

// Dimensions... so far untested in 2D
#define DIM 3



// Position data is stored in r[3*N*DIM], and a list of Polymers described how the N particles are partitioned
// into molecules with different properties
typedef struct{

  unsigned int  numAtoms, 
                firstAtomID, 
                perscription;     // NORMAL (=0), HELIX (=1), can add more later
  double        spring_constant,  // Currently unused
                radius,           // Currently unused
                h_w,              // helix frequency
                h_v,              // helix velocity
                h_R,              // helix radius
                h_l;              // helix pitch [cos(h_l z), sin(h_l z), z]
} Polymers;

typedef struct{

  int total_time, hydro;
  double temperature;

} Params;

typedef struct{

  double mean_Rg, mean_r_com;

} Results;


#endif
