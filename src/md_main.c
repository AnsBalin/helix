#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "hdr/md_init.h"
#include "hdr/md_defs.h"
#include "hdr/md_simulation.h"
#include "hdr/md_algebra.h"

void helix_singlePolymer();
void helix_manyPolymers();
void N_colloids( int N, double L );

int main(){

	N_colloids( 100, 10 );

	return 0;
}

void helix_singlePolymer(){

	int numPolymers = 2;
	Polymers* Ply = malloc( numPolymers*sizeof(Polymers) );

	
	/* Create polymer */
	int polymerN = 10;
	Ply[0].numAtoms = polymerN;
	Ply[0].firstAtomID = 0;
	Ply[0].perscription = NORMAL;

	double* r_init_poly = malloc( DIM*polymerN*sizeof(double) );
	/* Dz =  */
	double r0[3] = {10.0, 0.0, 10.0};
	double Dr[3] = {1.0, 0.0, 0.0};
	placePolymer( Ply, SAW, r_init_poly, r0, Dr, 0.0 );


	/* Create helix */
	int helixN = 150;
	Ply[1].numAtoms = helixN;
	Ply[1].firstAtomID = polymerN;
	Ply[1].perscription = HELIX;
	Ply[1].h_w = 4.0;
	Ply[1].h_v = 0.0;
	Ply[1].h_R = 3.0;
	Ply[1].h_l = MD_TWOPI/20.0;

	double* r_init_helix = malloc( DIM*helixN*sizeof(double) );
	placePolymer( Ply+1, HELIX, r_init_helix, r0, Dr, 0.0 );

	double* r_init = malloc( DIM*(helixN + polymerN)*sizeof(double) );
	
	int i = 0;
	for (int n = 0; n < DIM*polymerN; ++n) 	{
		r_init[i] 	= r_init_poly[n];
		i++;

	}
	for (int n = 0; n < DIM*helixN; ++n) 	{
		r_init[i] 	= r_init_helix[n];
		i++;

	}


	Params parameters = { .total_time = 100000, .hydro = 1, .temperature = 1};

	simulation( Ply, parameters, r_init, numPolymers, 1);

	free(Ply);
	free(r_init);

}

void N_colloids( int N, double L ) {


	int numPolymers = 1;
	Polymers* colloids = malloc( numPolymers*sizeof(Polymers) );

	colloids[0].numAtoms = N;
	colloids[0].firstAtomID = 0;
	colloids[0].perscription = NORMAL;

	double* r_init = malloc( DIM*N*sizeof(double) );

	double r0[3] = {10.0, 0.0, 10.0}; //does nothing but needed for placePolymer
	double Dr[3] = {1.0, 0.0, 0.0}; // DOes nothing in this case
	placePolymer( colloids, COL, r_init, r0, Dr, L );

	Params parameters = { .total_time = 10000, .hydro = 0, .temperature = 1, .size=L };
	printf("hi...\n");
	simulation( colloids, parameters, r_init, numPolymers, 1);

	free(colloids);
	free(r_init);

} 