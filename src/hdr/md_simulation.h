#ifndef __MD_SIMULATION_H__
#define __MD_SIMULATION_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "md_defs.h"
#include "md_init.h"
#include "md_forces.h"
#include "md_io.h"

void simulation( Polymers* Ply, Params parameters, double* r_init, int numPolymers, int simnum  );
void simulation2( Polymers* Ply, Params parameters, double* r_init, int numPolymers, int simnum, double* sqDisp, double* Rg );
void simulation3( Polymers* Ply, Params parameters, double* r_init, int numPolymers, int simnum );
double update( Polymers* Ply, int numPolymers, double* r, double* f, double* dw, double* r_ij, double* D, double* B, double* Df, double* Bdw, int N_tot, int HYDRO, double t );
void flowfield( int numAtoms, double* r, double* f, double* vx0, double* vy0, double* vz0, int M, int L  );
void flowfield3( int numAtoms, double* r, double* f, double* vx0, double* vy0, double* vz0, int N, int M, int L  );
void meanforce( double* r, double* f, int numAtoms, FILE* fp );
void calculate_CoM( Polymers* Ply, double* r, double* r_com );
double calculate_Rg( Polymers* Ply, double* r, double* r_com );
double calculate_Rg2( Polymers* Ply, double* r, double* r_com );
void shift( double* r, double* r_com, int N );
double calculate_SqDisplacement( double* r_com );

#endif
