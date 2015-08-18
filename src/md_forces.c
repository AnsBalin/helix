#include "hdr/md_forces.h"
#include "hdr/md_defs.h"

void computeForces( Polymers* PlyList, int numPolymers, double* r, double* f, double* r_ij, int N ){
	
	
	int r_ij_fail;
	r_ij_fail = compute_r_ij( r, r_ij, N );

	if( r_ij_fail ){
		printf("r_ij fail\n");
	}

	/* 	Old repel function exerts mutual force between all monomers irrespective of type
		New repel function will push normal monomers off perscribed ones (with twice the force) */
 	
 	// for (int i = 0; i < N; ++i)
	// {
	// 	for (int j = i+1; j < N; ++j)
	// 	{		
	// 		if( r_ij[ squ(i,j,N) ] <= MD_q0 ){
	// 			repel_old( i, j, r_ij[squ(i,j,N)], f, r );
	// 		}
	// 	}
	// }
 	
 	repel( PlyList, numPolymers, r, f, r_ij, N  );
 	attract( PlyList, numPolymers, r, f, r_ij, N );

 	/* We need to add force f = 6*pi*mu v to monomers perscribed to move at speed v */
 	prescribedForces( PlyList, numPolymers, r, f, r_ij, N );

}

void repel_old( int i, int j, double r_ij, double* f, double* r ){

	/* This function is depricated and does not support simulations containing
		monomers with perscribed motion. */

	double* r_diff = malloc( DIM*sizeof(double));
	double force;
	double r_ij_renorm = r_ij;

	force = softsphere( r_ij_renorm );
	//force = (force>hook(2*MD_q0)) ? hook(2*MD_q0) : force;
	

	// r_diff points from rj to ri
	FOR_ALL_K {
		r_diff[k] 		 =  r[ DIM*i + k ] - r[ DIM*j + k ];
		f[ DIM*i + k ] 	+=  force*r_diff[k];
		f[ DIM*j + k ] 	+= -force*r_diff[k];
	}

	free(r_diff);

}

void repel( Polymers* PlyList, int numPolymers, double* r, double* f, double* r_ij, int N  ){

	double* r_diff = malloc( DIM*sizeof(double));
	double force, separation;
	int numAtomsP, numAtomsQ, perscribedP, perscribedQ, p1, p2;
	
	for (int p = 0; p < numPolymers; ++p){
 		
 		numAtomsP = (PlyList+p)->numAtoms;
 		perscribedP = (PlyList+p)->perscription;
 		
 		for (int q = p; q < numPolymers; ++q){
 			
 			numAtomsQ = (PlyList+q)->numAtoms;
 			perscribedQ = (PlyList+q)->perscription;

 			if ( !perscribedP || !perscribedQ ){
	 			for (int i = 0; i < numAtomsP; ++i){
	 				p1 = (PlyList+p)->firstAtomID + i;
 					for (int j = ((p==q)? i+1 : 0); j < numAtomsQ; ++j){
 						
 						p2 = (PlyList+q)->firstAtomID + j;
 						separation = r_ij[ squ(p1,p2,N) ];
 						if( separation <= MD_q0 ){
 							force = softsphere( separation );
 							

							FOR_ALL_K {
								r_diff[k] = r[ DIM*p1 + k ] - r[ DIM*p2 + k ];
								/* Apologies for this somewhat complicated logic */
								f[ DIM*p1 + k ] +=  (!perscribedP)*(2 - !perscribedQ)*force*r_diff[k];
								f[ DIM*p2 + k ] += -(!perscribedQ)*(2 - !perscribedP)*force*r_diff[k];
							}
						}
 					}
 				}
 			}
 		

 		}
 	}
 	free(r_diff);
}

void attract( Polymers* PlyList, int numPolymers, double* r, double* f, double* r_ij, int N ){

	double* r_diff = malloc( DIM*sizeof(double));
	double force;
	int p1, p2;
	int numAtoms, perscribed;

	for (int p = 0; p < numPolymers; ++p)
 	{
 		numAtoms = (PlyList+p)->numAtoms;
 		perscribed = (PlyList+p)->perscription;
 		if( !perscribed ){
 			for (int i = 0; i < numAtoms-1; ++i)
 			{

 				p1 = (PlyList+p)->firstAtomID + i;
 				p2 = p1 + 1;

 				force = fene( r_ij[ N*p1 + p2 ] );
 				
				for (int k = 0; k < DIM; ++k)
				{
					r_diff[k]  =  r[  DIM*p1 + k ] - r[  DIM*p2 + k];
					f[ DIM*p1 + k ] +=  force*r_diff[k];
					f[ DIM*p2 + k ] += -force*r_diff[k];
				}
 			}
 		}
 	}

 	free(r_diff);
}

void prescribedForces( Polymers* PlyList, int numPolymers, double* r, double* f, double* r_ij, int N ){

// 	for (int p = 0; p < numPolymers; ++p)
//  	{
//  		numAtoms = (PlyList+p)->numAtoms;
//  		if( (PlyList+p)->FORCE ){
//  			for (int i = 0; i < numAtoms-1; ++i)
//  			{

//  				p1 = (PlyList+p)->firstAtomID + i;
	
// 				switch( (PlyList+p)->PERSCRIPTION ){

// 					case HELIX:



// 					case default:



// 				} 				


//  			}
//  		}
//  	}

}


void calc_dw( Polymers* PlyList, int numPolymers, double* dw, int N_tot ){

	/* dw is a Gaussian vector with:
			mean = 0.0
			var  = 1.0 

		For monomers with perscribed motion, dw = 0 */
	
	int numAtoms, perscribed;
	int n=0;

	for (int p = 0; p < numPolymers; ++p)
	{
		
		perscribed = (PlyList+p)->perscription;
		
		if ( !perscribed )
		{
			numAtoms = (PlyList+p)->numAtoms;
			for (int i = 0; i < numAtoms; ++i)
			{
				FOR_ALL_K dw[DIM*n+k] = gaussian( 1.0 );
				n++;
			}
		}
		else{
			numAtoms = (PlyList+p)->numAtoms;
			for (int i = 0; i < numAtoms; ++i)
			{
				FOR_ALL_K dw[DIM*n+k] = 0.0;
				n++;
			}
		}

		
	}

	

}


double gaussian( double sigma ){

  /* 	Applies Box-Muller transform to generate gaussian white noise
   		Adapted from wikipedia page (March 2015) */
  
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

  z0 = sqrt( -2.0 * log(u1) ) * cos( MD_TWOPI * u2 ); 
  z1 = sqrt( -2.0 * log(u1) ) * sin( MD_TWOPI * u2 );

  return z0 * sigma;

}

double softsphere( double r ){

	double rri,rri3;
  
	rri = 1./(r*r);
	rri3 = rri*rri*rri;
  	return 48. * rri3 * (rri3 - 0.5) * rri;

}

int compute_r_ij( double* r, double* r_ij, int numAtoms ){
	
	/*  NOTE Consider triangular for loop so symmetric elements of r_ij need not be computed twice. 
		In fact... if j>i always in the future, lower elements need not be computed at all...? */

	double *ri, *rj, *r_diff;
	double rix, riy, riz, rjx, rjy, rjz;
	r_diff = malloc( DIM*sizeof(double));
	for (int i = 0; i < numAtoms; ++i)
	{	
		rix = r[DIM*i];
		riy = r[DIM*i+1];
		riz = (DIM==3 ? r[DIM*i+2] : 0.0);
		for (int j = 0; j < numAtoms; ++j)
		{
			rjx = r[DIM*j];
			rjy = r[DIM*j+1];
			rjz = (DIM==3 ? r[DIM*j+2] : 0.0);
			r_ij[ squ(i,j,numAtoms) ] = sqrt( (rix-rjx)*(rix-rjx) + (riy-rjy)*(riy-rjy) + (riz-rjz)*(riz-rjz) );

		}

	}
	free(r_diff);

	return 0;

}

double hook( double r ){
  
  /* Hookean spring F = -k*r */
  /*  NOTE Currently no spring constant implemented */
  return -r;
   
}

double fene( double rr ){

	/* 	NOTE I just noticed that f_max is const for all runtime. Consider making rr_0 and f_max both
		global constants to avoid revaluation of f_max on every call. */
	const static double rr_0 = 1.5*1.5;
	double f_max = 30*rr_0*sqrt(0.6*rr_0)/(0.6*rr_0-rr_0);
  	return rr > 0.6*rr_0 ? f_max : 30*rr_0*sqrt(rr)/(rr-rr_0);

}

void calcD( double* r_arr, double* r_ij, double* D, int N, int TENSOR ){

	/*  */
	double a, r, rr, r1, r2, C1, C2, C3;
	
	a = MD_q0/2.0; // a is the 'radius' where MD_q0 is the separation that minimises WCA 

	switch ( TENSOR ){
		case NOHI:

			/* NOTE for NOHI we don't need to evaluate this every timestep. We can initialise D as diag
			   and then add the corrections only in the +HI cases */
			for (int i = 0; i < DIM*N; ++i)
			{
				D[squ(i,i,DIM*N)] = (4.0 / (3.0*a));
			}
			break; 

		/* I got rid of case OSEEN altogether... we will never use it */ 
		case OSEEN:

		case ROTNE2:
			/* TYLER Add ROTNE2 here, see below for 3D implementation */

		case ROTNE3:

			for (int i = 0; i < N; ++i){
				for (int n = 0; n < DIM; ++n){
					for (int j = 0; j < N; ++j){

						r = r_ij[ squ(i,j,N) ];
						rr =  r*r;

						C1 = 1 + 2.0*a*a/(3.0*rr);
						C2 = 1 - 2.0*a*a/rr;
						C3 = 1 - 9.0*r/(32.0*a);

						for (int m = 0; m < DIM; ++m)
						{
							r1 = r_arr[DIM*i+n] - r_arr[DIM*j+n];
							r2 = r_arr[DIM*i+m] - r_arr[DIM*j+m];
							if( r >= 2.0*a ) 	D[rank4(i,j,n,m,N)] = I(i,j) ? I(n,m) : (1-I(i,j))*( 3.0*a/(4.0*r) )*( C1*I(n,m) + C2*r1*r2/rr );
							else 				D[rank4(i,j,n,m,N)] = I(i,j) ? I(n,m) : (1-I(i,j))*( C3*I(n,m) + (3.0/(32.0*a))*r1*r2/r );

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

	/* 	Uses the NAG cholesky decomposition to ensure <B.B> = D 
		Documentation: http://www.nag.co.uk/numeric/cl/manual/pdf/F07/f07fdc.pdf 
		Returns 1 if D is not positive definite, 0 for success */
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

	/* NOTE Consider passing B directly to Cholesky? This might marginally reduce overhead */
	for (int i = 0; i < n*n; ++i)
	{
		B[i] = a[i];
	}
	free(a);
	return fail.code == NE_POS_DEF;
}

