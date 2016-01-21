#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "hdr/md_init.h"
#include "hdr/md_defs.h"
#include "hdr/md_simulation.h"
#include "hdr/md_algebra.h"

void md_init_SinglePolymer();
void simple_test(int in_N, double in_R, double in_l, double in_w, double in_v );
void single_polymer_test( int in_simnum );
void CholeskyTest();
void matrixMultiplyTest();

int main( int argc, char *argv[]  ){

	//int in_simnum;
	int in_simnum = 1;
	//double in_R=10.0, in_l=1.0, in_w=0.0, in_v=0.0;
	//sscanf(argv[2],"%d",&in_N);
	//sscanf(argv[3],"%lf",&in_R);
	//sscanf(argv[4],"%lf",&in_w);
	//sscanf(argv[5],"%lf",&in_l);	
	//sscanf(argv[1],"%lf",&in_v);
	//printf("#\t%d\t%lf\t%lf\t%lf\t./git%lf\n", in_N, in_R, in_w, in_l, in_v);
	//simple_test( in_N, in_R, in_l, in_w, in_v );
	sscanf( argv[1], "%d", &in_simnum );
	single_polymer_test( in_simnum );

	return 0;
}

// void md_init_SinglePolymer()
// {

// 	double *r, *dr, *f, *dw, *D, *B, *rr;
// 	const int N_tot = 2;
// 	int numPolymers = 1;

// 	srand(time(NULL));
// 	allocAll( &r, &f, &dw, &D, &B, &rr, N_tot );

// 	Polymers* Ply = malloc( numPolymers*sizeof(Polymers) );
// 	int id = 0;
// 	Ply->numAtoms = N_tot/numPolymers;
// 	Ply->FORCE = 1;
// 	Ply->firstAtomID = id;
// 	id += Ply->numAtoms;
// 	Ply->spring_constant = 1.0;
// 	Ply->radius = 1.0;

// 	double r0[3] = {0.0, 0.0, 0.0};
// 	double Dr[3] = {1.0, 0.0, 0.0};
// 	placePolymer( Ply, SAW, r, r0, Dr, 0.0 );

// 	for (int i = 0; i < N_tot; ++i)
// 	{
// 		printf("%lf\t%lf\t%lf\n", r[DIM*i], r[DIM*i+1], r[DIM*i+2] ); 
// 	}

// 	free(r);
// 	free(dr);
// 	free(f);
// 	free(D);
// 	free(B);
// 	free(rr);
// 	free(Ply);

//}

void single_polymer_test( int simnum ){

	int numPolymers = 2;
	Polymers* Ply = malloc( numPolymers*sizeof(Polymers) );

	/* Create helix */
	// Vectors not needed for helix but needed for polymer
	double r0[3] = {15.0, 0.0, 0.0};
	double Dr[3] = {0.0, 0.0, 0.0};

	int helixN = 120;
	Ply[0].numAtoms = helixN;
	Ply[0].firstAtomID = 0;
	Ply[0].perscription = HELIX;
	Ply[0].h_w = 20.0;
	Ply[0].h_v = 0.0;
	Ply[0].h_R = 4.0;
	Ply[0].h_l = MD_TWOPI/15;

	double* r_init_helix = malloc( DIM*helixN*sizeof(double) );
	placePolymer( Ply, HELIX, r_init_helix, r0, Dr, 0.0 );


	int polymerN = 25;
	Ply[1].numAtoms = polymerN;
	Ply[1].firstAtomID = helixN;
	Ply[1].perscription = NORMAL;
	

	double* r_init_poly = malloc( DIM*polymerN*sizeof(double) );
	placePolymer( Ply+1, SAW, r_init_poly, r0, Dr, 0.0 );

	double* r_init = malloc( DIM*(helixN + polymerN)*sizeof(double) );
	
	int ii = 0;
	for (int i = 0; i < DIM*helixN; ++i)
	{
		r_init[ii] = r_init_helix[i];
		++ii;
	}
	for (int i = 0; i < DIM*polymerN; ++i)
	{
		r_init[ii] = r_init_poly[i];
		printf("%f\n", r_init[i]);
		++ii;
	}
	
	Params parameters = { .total_time = 1000000, .hydro = 1 , .temperature = 1};

	simulation( Ply, parameters, r_init, numPolymers, simnum);

	free(Ply);
	free(r_init_helix);
	free(r_init_poly);
	free(r_init);

}

void simple_test(int in_N, double in_R, double in_l, double in_w, double in_v ){

	int numPolymers = 1;
	Polymers* Ply = malloc( numPolymers*sizeof(Polymers) );

	
	/* Dz =  */
	double r0[3] = {15.0, 15.0, 20.0};
	double Dr[3] = {1.0, 0.0, 0.0};

	/* Create helix */
	int helixN = in_N;
	Ply[0].numAtoms = helixN;
	Ply[0].firstAtomID = 0;
	Ply[0].perscription = HELIX;
	Ply[0].h_w = in_w;
	Ply[0].h_v = in_v;
	Ply[0].h_R = in_R;
	Ply[0].h_l = MD_TWOPI/in_l;

	double* r_init_helix = malloc( DIM*helixN*sizeof(double) );
	placePolymer( Ply, HELIX, r_init_helix, r0, Dr, 0.0 );

	
	/* Create plane of tracers */
	/*int tracersN = 6;
	Ply[1].numAtoms = tracersN*tracersN;
	Ply[1].firstAtomID = helixN;
	Ply[1].perscription = TRACER;

	double* r_init_tracer = malloc( DIM*tracersN*tracersN*sizeof(double) );
	placePolymer( Ply+1, TRACER, r_init_tracer, r0, Dr, 0.0 );

	
	

	
	for (int i = 0; i < DIM*(tracersN*tracersN); ++i)
	{
		r_init[ii] = r_init_tracer[i];
		printf("%f\n", r_init[i]);
		++ii;
	}
	*/

	/* Create polymer */
	int polymerN = 0;
	//Ply[1].numAtoms = polymerN;
	//Ply[1].firstAtomID = helixN;
	//Ply[1].perscription = NORMAL;

	//double* r_init_poly = malloc( DIM*polymerN*sizeof(double) );
	//placePolymer( Ply+1, SAW, r_init_poly, r0, Dr, 0.0 );

	double* r_init = malloc( DIM*(helixN + polymerN)*sizeof(double) );
	int ii = 0;
	for (int i = 0; i < DIM*helixN; ++i)
	{
		r_init[ii] = r_init_helix[i];
		++ii;
	}
	/*for (int i = 0; i < DIM*polymerN; ++i)
	{
		r_init[ii] = r_init_poly[i];
		printf("%f\n", r_init[i]);
		++ii;
	}*/
	Params parameters = { .total_time = 10000, .hydro = 1 , .temperature = 1};

	simulation( Ply, parameters, r_init, numPolymers, 1);

	free(Ply);
	free(r_init_helix);
	//free(r_init_poly);
	free(r_init);

}


void CholeskyTest(){

	double *U, *UT, *D, *B;
	int N = 4;
	U  = malloc( N*N*sizeof(double) );
	UT = malloc( N*N*sizeof(double) );
	D  = malloc( N*N*sizeof(double) );
	B  = malloc( N*N*sizeof(double) );

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			U[squ(i,j,N)]  = 0.0;
			UT[squ(j,i,N)] = 0.0;
			printf("\t\t");
		}
		for (int j = i; j < 4; ++j )
		{
			U[squ(i,j,N)] = i + (double)j/((double)i+0.1) + 0.1;
			UT[squ(j,i,N)] = U[squ(i,j,N)];
			printf("%lf\t",U[squ(i,j,N)] );
		}
		printf("\n");
	}
	printf("\n");
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			for (int m = 0; m < N; ++m)
			{
				D[squ(i,j,N)] += UT[squ(i,m,N)] * U[squ(m,j,N)];
			}
			printf("%lf\t", D[squ(i,j,N)]);
		}
		printf("\n");
	}
	printf("\n");
	calcB(D,N,B);


	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			printf("%lf\t", B[squ(i,j,N)]);
		}
		printf("\n");
	}

	free(U);
	free(UT);
	free(D);
	free(B);

}

void matrixMultiplyTest(){


	int N = 3;

	double* M = malloc(N*N*sizeof(double));
	double* x = malloc(N*sizeof(double));
	double* y = malloc(N*sizeof(double));

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			M[squ(i,j,N)] = i-j+1;
		}
		x[i] = i;
		y[i] = 0.0;

	}

	mat_multiply_A_x( M, N, x, y );

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{

			printf("%f\t", M[squ(i,j,N)]);
		}
		printf("\t%f\t\t%f\n", x[i], y[i]);
		
	}

}
