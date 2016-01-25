#ifndef __MD_INIT_H__
#define __MD_INIT_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "md_defs.h"

void allocAll( double** r,  double** f, double** dw, double** D, double** B, double** Df, double** Bdw, double** rr, const int N );
void placeAtom( int atom, double* r, double* r_new );
void placePolymer( Polymers* Ply, int SHAPE, double* r, double* r0, double* dr, double param1 );
int overlap( Polymers* Ply, int atom, double* r_test, double* r );

#endif
