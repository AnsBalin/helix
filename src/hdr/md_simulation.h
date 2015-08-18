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

void simulation( Polymers* Ply, Params parameters, double* r_init, int numPolymers, int simnum );
void update( Polymers* Ply, int numPolymers, double* r, double* f, double* dw, double* r_ij, double* D, double* B, int N_tot, int HYDRO );
void calculate_CoM( Polymers Ply, double* r, double* r_com );
double calculate_Rg( Polymers Ply, double* r, double* r_com );


#endif
