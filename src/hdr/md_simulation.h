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

void simulation( Polymers* Ply, Params parameters, double* r_init, int numPolymers, int simnum, int poly1, int polyn  );
void update( Polymers* Ply, int numPolymers, double* r, double* f, double* dw, double* r_ij, double* D, double* B, double* Df, double* Bdw, int N_tot, int HYDRO, double t, int poly1, int polyn );
void flowfield( int numAtoms, double* r, double* f );
void meanforce( double* r, double* f, int numAtoms );
void calculate_CoM( Polymers Ply, double* r, double* r_com );
double calculate_Rg( Polymers Ply, double* r, double* r_com );


#endif
