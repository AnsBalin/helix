#ifndef __MD_ALGEBRA_H__
#define __MD_ALGEBRA_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "md_defs.h"

void mat_multiply_sym_A_x( double* mat, int N, double* x, double* y );
void mat_multiply_A_x( double* mat, int N, double* x, double* y );
int mat_threshhold( double* mat, int N, double thresh, int** overlap_inds );
void mat_outer_product( double* x, double* y, double* xy, int N );



#endif
