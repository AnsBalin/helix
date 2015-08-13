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

void md_init_SinglePolymer()
{

	double *r, *dr, *f, *dw, *D, *B, *rr;
	const int N_tot = 2;
	int numPolymers = 1;

	srand(time(NULL));
	allocAll( &r, &f, &dw, &D, &B, &rr, N_tot );

	Polymers* Ply = malloc( numPolymers*sizeof(Polymers) );
	int id = 0;
	Ply->numAtoms = N_tot/numPolymers;
	Ply->NOISE = 1;
	Ply->FORCE = 1;
	Ply->HYDRO = 1;
	Ply->firstAtomID = id;
	id += Ply->numAtoms;
	Ply->spring_constant = 1.0;
	Ply->radius = 1.0;

	double r0[3] = {0.0, 0.0, 0.0};
	double Dr[3] = {1.0, 0.0, 0.0};
	placePolymer( Ply, SAW, r, r0, Dr, 0.0 );

	for (int i = 0; i < N_tot; ++i)
	{
		printf("%lf\t%lf\t%lf\n", r[DIM*i], r[DIM*i+1], r[DIM*i+2] ); 
	}

	free(r);
	free(dr);
	free(f);
	free(D);
	free(B);
	free(rr);
	free(Ply);

}

void simple_test(){

//	int N_arr[7] = {10, 20, 40, 80, 160, 320, 640};

	int N_arr[6] = {160, 80, 40, 20, 10, 5};
	int N_tot;
	Results results; 
	FILE *fp;
	fp = fopen("dat/Rg_r_com.dat", "w");
	for (int nn = 0; nn < 6; ++nn)
	{

		double mean_Rg = 0.0;
		double mean_r_com = 0.0;
		N_tot = N_arr[ nn ];
		//N_tot = 320;
		for (int sim_number = 0; sim_number < 15; ++sim_number)
		{
			
			printf("%d, %d\n", N_tot, sim_number);
			int numPolymers = 1;
			srand(time(NULL)+sim_number);

			Polymers* Ply = malloc( numPolymers*sizeof(Polymers) );
			int id = 0;
			Ply->numAtoms = N_tot/numPolymers;
			Ply->NOISE = 1;
			Ply->FORCE = 1;
			Ply->HYDRO = 0;
			Ply->firstAtomID = id;
			id += Ply->numAtoms;
			Ply->spring_constant = 1.0;
			Ply->radius = 1.0;
			
			double* r_init = malloc( DIM*N_tot*sizeof(double) );
			double r0[3] = {0.0, 0.0, 0.0};
			double Dr[3] = {1.0, 0.0, 0.0};
			placePolymer( Ply, SAW, r_init, r0, Dr, 0.0 );

			for (int i = 0; i < N_tot; ++i)
			{
				//printf("%lf\t%lf\t%lf\n", r_init[DIM*i], r_init[DIM*i+1], r_init[DIM*i+2] ); 
			}


			Params parameters = { .total_time = 100000, .HYDRO =1, .temperature=1 };
			
			//printf("Starting simulation...\n");
			results = simulation( Ply, parameters, r_init, 1, sim_number );
			printf("simulation() finished\n");
			//fflush(stdout);
			mean_Rg += 0.001*results.mean_Rg;
			mean_r_com += 0.001*results.mean_r_com;
			fprintf(fp,"%d\t%lf\t%lf\n", N_tot, results.mean_Rg, results.mean_r_com);

			free(Ply);
			free(r_init);
		}
	}

	fclose(fp);
	

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