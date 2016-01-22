#ifndef __MD_FORCES_H__
#define __MD_FORCES_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "/opt/NAG/clmi623dgl/include/nag.h"
//#include "/opt/NAG/clmi623dgl/include/nag_stdlib.h"
//#include "/opt/NAG/clmi623dgl/include/nagf07.h"
#include "nag.h"
#include "nag_stdlib.h"
#include "nagf07.h"
//#include <nagx04.h>
#include "md_defs.h"
#include "md_algebra.h"

//void computeForces( Polymers* PlyList, int numPolymers, double temperature );
void computeForces( Polymers* PlyList, int numPolymers, double* r, double* f, double* r_ij, int N, double t );
void repel_old( int i, int j, double r_ij, double* f, double* r );
void repel( Polymers* PlyList, int numPolymers, double* r, double* f, double* r_ij, int N  );
void prescribedForces( Polymers* PlyList, int numPolymers, double* r, double* f, double* r_ij, int N, double t );
void attract( Polymers* PlyList, int numPolymers, double* r, double* f, double* r_ij, int N );
void bending( Polymers* PlyList, int numPolymers, double* r, double* f, double* r_ij, int N );
void calc_dw( Polymers* PlyList, int numPolymers, double* dw, int N_tot );
double gaussian( double sigma );
double softsphere( double r );
int compute_r_ij( double* r, double* r_ij, int numAtoms );
double hook( double r );
double fene( double rr );
void calcD( double* r, double* r_ij, double* D, int N, int TENSOR );
int calcB( double* D, int N, double* B, int poly1, int polyn);

#endif