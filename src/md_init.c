#include "hdr/md_init.h"

void allocAll( double** r, double** f, double** dw, double** D, double** B, double** rr, const int N ){

	(*r) = malloc( DIM*N*sizeof(double) );
	(*f) = malloc( DIM*N*sizeof(double) );
	(*dw) = malloc( DIM*N*sizeof(double) );
	(*D) = malloc( DIM*DIM*N*N*sizeof(double) );
  (*B) = malloc( DIM*DIM*N*N*sizeof(double) );
	(*rr) = malloc( N*N*sizeof(double) );


}



void placeAtom( int atom, double* r, double* r_new ){
  
  double* r_tmp;
  r_tmp = r + atom*DIM;
  FOR_ALL_K r_tmp[k] = r_new[k];


}

void placePolymer( Polymers* Ply, int SHAPE, double* r, double* r0, double* dr, double param1 ){

  
  double* r_new = (double*) malloc( DIM*sizeof(double) );
  double theta,phi;

  switch(SHAPE)
  {

    case LINE:
      
      // Place atom i at position r = Q0 + i*dr
      for( int i=0; i<Ply->numAtoms; ++i ){
        
        //printf("%d\n", i);
        FOR_ALL_K r_new[k] = r0[k] + i*dr[k]; //+ 0.1*((double)rand()/(double)RAND_MAX - 0.5);

        placeAtom( Ply->firstAtomID + i, r, r_new );
      
      }

      break;

    case HELIX:

      break;

    case SAW:


      placeAtom( Ply->firstAtomID, r, r0);
      for( int i=1; i<Ply->numAtoms; ++i ){
        do {
          phi = TWOPI * (double)rand()/(double)RAND_MAX;
          theta=0;
          if( DIM==3 ) {
            dr[2] = q0*( 2.0*(double)rand()/(double)RAND_MAX - 1.0);
            theta = asin(dr[2]/q0);
          }
          else printf("DIM must be 2 or 3 for SAW.\n");
          dr[0] = q0*cos(theta)*cos(phi);
          dr[1] = q0*cos(theta)*sin(phi);
          FOR_ALL_K r_new[k] = r[(i-1)*DIM + k] + dr[k];
          
        }while (overlap( Ply, i, r_new, r ) != 0);
        
        placeAtom( Ply->firstAtomID + i, r, r_new );
      }

      break;

    default:
      printf("Did not supply correct SHAPE.\n");
      break;
  }
  free(r_new);
}

int overlap( Polymers* Ply, int atom, double* r_test, double* r ){
  
  // Returns the number of particles a particle at r_test overlaps with 

  int result = 0;
  int i,k;
  double* r_i;
  double r2;

  for( i=0; i<atom; i++  ){

    r2 = 0.0;
    r_i = r + i*DIM;
    FOR_ALL_K r2 += ( r_test[k] - r_i[k] )*( r_test[k] - r_i[k] ); 
    if( sqrt(r2)<1.0*q0 ) result++;

  }
  return result;

}





