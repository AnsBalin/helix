#ifndef __MD_IO_H__
#define __MD_IO_H__

#include <stdio.h>
#include <stdlib.h>
#include "md_defs.h"

void printAtomR( int atom, double* r );
void saveRtofile( double* r, int N_tot, FILE* fp );
void saveXYZtofile( Polymers* PlyList, int numPolymers, double* r, int N_tot, FILE* fp );
void saveStatstofile( double t, double* r_com, double Rg, FILE* fp );
#endif

