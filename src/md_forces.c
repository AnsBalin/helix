#include "hdr/md_forces.h"
#include "hdr/md_defs.h"

void computeForces( Polymers* PlyList, int numPolymers, double* r, double* f, double* r_ij, int N ){
	
	
	
	int r_ij_fail;
	r_ij_fail = compute_r_ij( r, r_ij, N );

	if( r_ij_fail ){
		printf("r_ij fail\n");
	}


 	for (int i = 0; i < N; ++i)
	{
		for (int j = i+1; j < N; ++j)
		{		

			if( r_ij[ squ(i,j,N) ] <= q0 ){
				repel( i, j, r_ij[squ(i,j,N)], f, r );
			}
		}
	}
 	
 	attract( PlyList, numPolymers, r, f, r_ij, N );

}

void repel( Polymers* PlyList, int numPolymers, double* r, double* f, double* r_ij, int N  ){

	double* r_diff = malloc( DIM*sizeof(double));
	double force;
	

	force = softsphere( r_ij_renorm );
	//force = (force>hook(2*q0)) ? hook(2*q0) : force;
	
	// r_diff points from rj to ri
	FOR_ALL_K r_diff[k] = r[ DIM*i + k ] - r[ DIM*j + k ];
	FOR_ALL_K f[ DIM*i + k ] +=  force*r_diff[k];
	FOR_ALL_K f[ DIM*j + k ] += -force*r_diff[k];

	for (int i = 0; i < N; ++i)
	{
		for (int j = i+1; j < N; ++j)
		{		

			if( r_ij[ squ(i,j,N) ] <= q0 ){
				force = softsphere( r_ij );
				FOR_ALL_K r_diff[k] = r[ DIM*i + k ] - r[ DIM*j + k ];
				FOR_ALL_K f[ DIM*i + k ] +=  force*r_diff[k];
				FOR_ALL_K f[ DIM*j + k ] += -force*r_diff[k];

			}
		}
	}


}

void attract( Polymers* PlyList, int numPolymers, double* r, double* f, double* r_ij, int N ){

	double* r_diff = malloc( DIM*sizeof(double));
	double force;
	int p1, p2;
	int numAtoms;

	for (int p = 0; p < numPolymers; ++p)
 	{
 		numAtoms = (PlyList+p)->numAtoms;
 		if( (PlyList+p)->FORCE ){
 			for (int i = 0; i < numAtoms-1; ++i)
 			{

 				p1 = (PlyList+p)->firstAtomID + i;
 				p2 = p1 + 1;

 				force = fene( r_ij[ N*p1 + p2 ] );
 				//if(r_ij[N*p1 + p2]>1.5) printf("%f\n", force);

 				FOR_ALL_K r_diff[k]  =  r[  DIM*p1 + k ] - r[  DIM*p2 + k];
				FOR_ALL_K f[ DIM*p1 + k ] +=  force*r_diff[k];
				FOR_ALL_K f[ DIM*p2 + k ] += -force*r_diff[k];
 			}
 		}
 	}
}

/*void addNoise( Polymers* PlyList, double* r, double* f, double* r_ij, int N ){

	double* r_diff = malloc( DIM*sizeof(double));
	double force;
	int p1;

	for (int p = 0; p < numPolymers; ++p)
 	{
 		numAtoms = (PlyList+p)->numAtoms;
 		if( (PlyList+p)->NOISE ){
 			for (int i = 0; i < numAtoms-1; ++i)
 			{
 				p1 = (PlyList+p)->firstAtomID + i;
 				FOR_ALL_K f[ DIM*p1 + k ] += gaussian( ((PlyList+p)->noise_factor) );

 			}
 		}
 	}


}*/


void calc_dw( double* dw, int N_tot ){

	for (int n = 0; n < 3*N_tot; ++n)
	{
		dw[n] = gaussian( 1.0 );
	}

}


double gaussian( double sigma ){

  // Applies Box-Muller transform to generate gaussian white noise
  // Adapted from wikipedia page (March 2015)
  

  static double z0, z1;
  static int generate=0;
  //double epsilon = DBL_MIN;
  generate++;

  if( generate % 2 == 0 )
    return z1 * sigma;

  double u1, u2;
  do{

    u1 = rand() * (1.0 / RAND_MAX);
    u2 = rand() * (1.0 / RAND_MAX);
  
  } while( u1 == 0.0 );

  z0 = sqrt( -2.0 * log(u1) ) * cos( TWOPI * u2 ); 
  z1 = sqrt( -2.0 * log(u1) ) * sin( TWOPI * u2 );

  return z0 * sigma;

}

double softsphere( double r ){

	double rri,rri3;
  
	rri = 1./(r*r);
	rri3 = rri*rri*rri;
  	return 48. * rri3 * (rri3 - 0.5) * rri;

}

int compute_r_ij( double* r, double* r_ij, int numAtoms ){
	
	double *ri, *rj, *r_diff;
	double rix, riy, riz, rjx, rjy, rjz;
	r_diff = malloc( DIM*sizeof(double));
	for (int i = 0; i < numAtoms; ++i)
	{
		//ri = r + DIM*i;
		rix = r[DIM*i];
		riy = r[DIM*i+1];
		riz = r[DIM*i+2];

		for (int j = 0; j < numAtoms; ++j)
		{
			rjx = r[DIM*j];
			rjy = r[DIM*j+1];
			rjz = r[DIM*j+2];

			//rj = r + DIM*j;
			//VSub( r_diff, ri, rj );
			//r_ij[ squ(i,j,numAtoms) ] = sqrt(VSqr( r_diff )); 
			r_ij[ squ(i,j,numAtoms) ] = sqrt( (rix-rjx)*(rix-rjx) + (riy-rjy)*(riy-rjy) + (riz-rjz)*(riz-rjz) );
			
			/*if(r_ij[squ(i,j,numAtoms)] > 2.0){
				free(r_diff);
				return squ(i,j,numAtoms) + 1;

			}*/

		}

	}
	free(r_diff);

	return 0;


}

double hook( double r ){
  
  /* Hookean spring F = -k*r */
  return -r;
   
}

double fene( double rr ){

	const static double rr_0 = 1.5*1.5;
	double f_max = 30*rr_0*sqrt(0.6*rr_0)/(0.6*rr_0-rr_0);
  	return rr > 0.6*rr_0 ? f_max : 30*rr_0*sqrt(rr)/(rr-rr_0);

}

void calcD( double* r_arr, double* r_ij, double* D, int N, int TENSOR ){

	double a,r,rr, r1, r2,C1,C2,C3;
	a = q0/2.0;

	switch ( TENSOR ){

		case NOHI:

			for (int i = 0; i < DIM*N; ++i)
			{
				D[squ(i,i,DIM*N)] = (4.0 / (3.0*a));
			}
			break; 

		/* I got rid of case OSEEN altogether */ 
		case ROTNE:

			for (int i = 0; i < N; ++i)
			{
				for (int n = 0; n < DIM; ++n)
				{
					for (int j = 0; j < N; ++j)
					{

						
						r = r_ij[ squ(i,j,N) ];
						rr =  r*r;

						C1 = 1 + 2.0*a*a/(3.0*rr);
						C2 = 1 - 2.0*a*a/rr;
						C3 = 1 - 9.0*r/(32.0*a);

						for (int m = 0; m < DIM; ++m)
						{
							r1 = r_arr[DIM*i+n] - r_arr[DIM*j+n];
							r2 = r_arr[DIM*i+m] - r_arr[DIM*j+m];
							if( r >= 2.0*a ) {

								D[rank4(i,j,n,m,N)] = I(i,j) ? I(n,m) : (1-I(i,j))*( 3.0*a/(4.0*r) )*( C1*I(n,m) + C2*r1*r2/rr );
							}
							else{

								D[rank4(i,j,n,m,N)] = I(i,j) ? I(n,m) : (1-I(i,j))*( C3*I(n,m) + (3.0/(32.0*a))*r1*r2/r );
								
							}

							
						}
					}
				}
			}

			break;

		
		default:
			break;


	}

	
}

int calcB( double* D, int N, double* B ){
	double* a = malloc( N*N*sizeof(double) );
	NagError fail;
	INIT_FAIL(fail);
	Nag_OrderType order = Nag_RowMajor;
	Nag_UploType uplo = Nag_Upper;
	int n = N;
	int pda = n;

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			a[squ(i,j,n)] = 0.0;
		}
		for (int j = i; j < n; ++j)
		{
			a[squ(i,j,n)] = D[squ(i,j,n)];
		}
	}

	f07fdc( order, uplo, n, a, pda, &fail );

	/*
	if( fail.code == NE_POS_DEF ){

		printf("Error! Matrix not positive-definite!\n");

	}else{

		printf("Success!!!!\n");
	}*/

	for (int i = 0; i < n*n; ++i)
	{
		B[i] = a[i];
	}

	free(a);

	return fail.code == NE_POS_DEF;

}