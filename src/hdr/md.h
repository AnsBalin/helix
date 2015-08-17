#ifndef __MD_H__
#define __MD_H__

#include "md_defs.h"
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>

void printAtomQ( Atom* Atm );
void printAtomV( Atom* Atm );
void printAtomA( Atom* Atm );
void saveQtofile( Atom* AtmList, int GPOP, FILE* fp );
void allocAll( Atom** AtmList, Polymer** PlyList, int N, Spec* SpList );
void initAll( Polymer* PlyList, Spec* SpList, int N );
void zeroPolymer( Polymer* Ply );
void zeroAtom( Atom* Atm );
void placeAtom( Atom* Atm, double* Qnew );
void placePolymer( Polymer* Ply, int SHAPE, double* Q0, double* dr, double param1 );
int overlap( Polymer* Ply, int num, double* r );


#endif