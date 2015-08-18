#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "hdr/md_init.h"
#include "hdr/md_defs.h"
#include "hdr/md_simulation.h"
#include "hdr/md_algebra.h"

void md_init_SinglePolymer();
void simple_test();
void CholeskyTest();
void matrixMultiplyTest();

int main(){

	simple_test();

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

void simple_test(){

	int N_arr[6] = {160, 80, 40, 20, 10, 5};
	int N_tot;
	Results results; 
	for (int nn = 0; nn < 1 ; ++nn)
	{

		N_tot = N_arr[ nn ];
		N_tot = 50;
		for (int sim_number = 0; sim_number < 1; ++sim_number)
		{
			
			printf("%d, %d\n", N_tot, sim_number);
			int numPolymers = 1;
			srand(time(NULL)+sim_number);

			/* Create 1 polymer */
			Polymers* Ply = malloc( numPolymers*sizeof(Polymers) );
			int id = 0;
			Ply->numAtoms = N_tot/numPolymers;
			Ply->perscription = 0;
			Ply->firstAtomID = id;
			id += Ply->numAtoms;
			Ply->spring_constant = 1.0;
			Ply->radius = 1.0;
			
			/* Plase it at the origin, with unit spacing according to SAW algorithm */
			double* r_init = malloc( DIM*N_tot*sizeof(double) );
			double r0[3] = {0.0, 0.0, 0.0};
			double Dr[3] = {1.0, 0.0, 0.0};
			placePolymer( Ply, SAW, r_init, r0, Dr, 0.0 );

			/* Simulation parameters (NOTE temperature does nothing) */
			Params parameters = { .total_time = 100000, .hydro =1, .temperature=1 };
		
			simulation( Ply, parameters, r_init, 1, sim_number );
			printf("simulation() finished\n");

			free(Ply);
			free(r_init);
		}
	}

	
	

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