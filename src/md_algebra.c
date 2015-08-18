//#include "md_h"
#include "hdr/md_algebra.h"
#include "hdr/md_defs.h"
//#include "omp.h"

void mat_multiply_sym_A_x( double* mat, int N, double* x, double* y ){
	
	/* Takes NxN matrix A, and assumes it is symmetric */
	// This is marginally faster than mat_multiply_A_x() for N ~ O(1000)
	double* A = mat;
	double a;

	for (int i = 0; i < N; ++i)
	{
		for (int j = i; j < N; ++j)
		{
			a = A[ squ( i, j, N ) ];
			y[i] += a*x[j];
			y[j] += a*x[i];
		}
	}

}

void mat_multiply_A_x( double* mat, int N, double* x, double* y ){

	double* A = mat;
	double a;

	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			a = A[ squ( i, j, N ) ];
			y[i] += a*x[j];
			
		}
	}

}


int mat_threshhold( double* mat, int N, double thresh, int** overlap_inds_result ){

	/* This is from a bygone era and can largely be ignored */

	// Returns array containing indices i,j, s.t. i>j
	// for A_ij < thresh. expands array in amortised time
	// Only evaluate pair repulsions for these particles i,j
	double* A = mat;
	unsigned int overlap_inds_size = 512; // Number of pairs of indices to store
	int* overlap_inds =   (int*)  malloc( 2*overlap_inds_size * sizeof(int) );
	int* overlap_inds_tmp =(int*) malloc( 2*overlap_inds_size * sizeof(int) );
	//printf("overlap_inds: %p\noverlap_inds_tmp: %p\n*overlap_inds_result: %p\n", (void*)overlap_inds, (void*)overlap_inds_tmp, (void*)*(overlap_inds_result));
	unsigned int counter = 0;
	double* ij_arr;
	for (int i = 0; i < N; ++i)
	{
		for (int j = i+1; j < N; ++j)
		{
			if( A[ squ(i,j,N) ] <= thresh ){
				//printf("sdasdasdasda\n");
				overlap_inds[2*counter  ] = i;
				overlap_inds[2*counter+1] = j;				
				overlap_inds_tmp[2*counter  ] = i;
				overlap_inds_tmp[2*counter+1] = j;
				counter++;
				
				// If array limit is reached, double size of array
				if(counter==overlap_inds_size){

					overlap_inds_size = overlap_inds_size*2;
					
					free(overlap_inds);
					overlap_inds = malloc(2*overlap_inds_size * sizeof(int));
					
					for (int i = 0; i < overlap_inds_size; ++i)
					{
						overlap_inds[i] = overlap_inds_tmp[i];
					}
					free(overlap_inds_tmp);
					overlap_inds_tmp = malloc(2*overlap_inds_size * sizeof(int));
					for (int i = 0; i < overlap_inds_size; ++i)
					{
						overlap_inds_tmp[i] = overlap_inds[i];
					}

				}
			}
		}
	}
	
	*overlap_inds_result = malloc(2*counter*sizeof(int));
	for (int i = 0; i < 2*counter; ++i)
	{
		(*overlap_inds_result)[i] = overlap_inds[i];
		
	}
	//*overlap_inds_result = overlap_inds;
	//printf("overlap_inds: %p\noverlap_inds_tmp: %p\n*overlap_inds_result: %p\n", (void*)overlap_inds, (void*)overlap_inds_tmp, (void*)*(overlap_inds_result));

	free(overlap_inds);		
	free(overlap_inds_tmp);
	//printf("hello...\n");
	return counter;

}

void mat_outer_product( double* x, double* y, double* xy, int N ){

	
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			xy[ squ(i,j,N) ] += x[i]*y[j];
		}
	}

}