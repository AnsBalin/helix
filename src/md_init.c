#include "hdr/md_init.h"

void allocAll( double** r, double** f, double** dw, double** D, double** B, double** Df, double** Bdw, double** rr, const int N ){

	(*r) = malloc( DIM*N*sizeof(double) );
	(*f) = malloc( DIM*N*sizeof(double) );
	(*dw) = malloc( DIM*N*sizeof(double) );
	(*D) = malloc( DIM*DIM*N*N*sizeof(double) );
  (*B) = malloc( DIM*DIM*N*N*sizeof(double) );
  (*Df) = malloc( DIM*N*sizeof(double) );
  (*Bdw) = malloc( DIM*N*sizeof(double) );
	(*rr) = malloc( N*N*sizeof(double) );

}



void placeAtom( int atom, double* r, double* r_new ){
  
  double* r_tmp;
  r_tmp = r + atom*DIM;
  FOR_ALL_K r_tmp[k] = r_new[k];

}

void placePolymer( Polymers* Ply, int SHAPE, double* r, double* r0, double* dr, double param1 ){

  /* NOTE The program can only place 1 polymer in the simulation right now
      Need to fix this next. */
  
  double* r_new = (double*) malloc( DIM*sizeof(double) );
  double r_com[3] = {0.0, 0.0, 0.0};
  double r_shifted[3] = {0.0,0.0,0.0};
  double theta,phi, a, R, l, Dz, z;
  int ii = 0;

  double meanx = 0;
  double meany = 0;
  double meanz = 0;
  switch(SHAPE)
  {

    case LINE:
      
      /* Place atom i at position r = MD_q0 + i*dr */
      for( int i=0; i<Ply->numAtoms; ++i ){
        
        //printf("%d\n", i);
        FOR_ALL_K r_new[k] = r0[k] + MD_q0*i*dr[k]; //+ 0.1*((double)rand()/(double)RAND_MAX - 0.5);

        placeAtom( Ply->firstAtomID + i, r, r_new );
     
      }

      break;

    case HELIX:

      /* For inter-particle spacing of a, the z-spacing must be:
                    dz = a/√( 1 + (R*l)^2 )
      */
      a = MD_q0;
      R = (Ply->h_R);
      l = Ply->h_l;
      Dz = a / sqrt( 1 + R*R*l*l );

      

      z = -(Dz*Ply->numAtoms)/2;

      for (int n = 0; n < Ply->numAtoms; ++n)
      {
        r_new[0] = (R /*/sqrt( (1.0/25.0)*(Ply->h_w)*(Ply->h_w) + 1 )*/)*cos(l*z);
        r_new[1] = (R /*/sqrt( (1.0/25.0)*(Ply->h_w)*(Ply->h_w) + 1 )*/)*sin(l*z);
        r_new[2] = z;

        

        placeAtom( n, r, r_new );
        z += Dz;
      }

      break;

    case SAW:
      
      placeAtom( 0, r, r0);
      for( int i=1; i<Ply->numAtoms; ++i ){
        do {
          phi = MD_TWOPI * (double)rand()/(double)RAND_MAX;
          theta=0;
          if( DIM==3 ) {
            dr[2] = MD_q0*( 2.0*(double)rand()/(double)RAND_MAX - 1.0);
            theta = asin(dr[2]/MD_q0);
          }
          else printf("DIM must be 2 or 3 for SAW.\n");
          dr[0] = MD_q0*cos(theta)*cos(phi);
          dr[1] = MD_q0*cos(theta)*sin(phi);

          meanz += dr[2];
          meany += dr[1];
          meanx += dr[0];

          FOR_ALL_K r_new[k] = r[(i-1)*DIM + k] + dr[k];
          
          

        }while (overlap( Ply, i, r_new, r ) != 0);
        r_com[0] += (1.0/Ply->numAtoms)*r_new[0];
        r_com[1] += (1.0/Ply->numAtoms)*r_new[1];
        r_com[2] += (1.0/Ply->numAtoms)*r_new[2];
        placeAtom( i, r, r_new );
      }
      printf("%.3f %.3f %.3f\n", meanx, meany, meanz);

      meanx = 0.0;
      for( int i=0; i<Ply->numAtoms; ++i){

        FOR_ALL_K r_shifted[k] = r[DIM*i + k] - (r_com[k]-r0[k]);
        meanx += r_shifted[2];
        placeAtom(i,r,r_shifted);

      }
      printf("%.3f\n", meanx);
      break;


    case TRACER:
      
      for (int n = 0; n < (int) sqrt(Ply->numAtoms); ++n)
      { 
        for (int m = 0; m < (int) sqrt(Ply->numAtoms); ++m)
        {
          r_new[0] = (double) m - sqrt(Ply->numAtoms)/2;
          r_new[1] = (double) n - sqrt(Ply->numAtoms)/2;
          r_new[2] = 19.0;
          placeAtom( ii, r, r_new );
          printf("%d\t%f\n", ii,r_new[0]);
          ++ii;
        }
      }

      break;

    case RANDOM:
      do {
          
      FOR_ALL_K r_new[k] = 10.0*(2.0*(double)rand()/(double)RAND_MAX - 1.0);
          
      }while (overlap( Ply, Ply->firstAtomID, r_new, r ) != 0);
        
      placeAtom( Ply->firstAtomID, r, r_new );
      printf("%f\t%f\t%f\n", r_new[0], r_new[1], r_new[2]);
      
      break; 
    case RANDOM_SAW:

      do {
          
      FOR_ALL_K r_new[k] = 10.0*(2.0*(double)rand()/(double)RAND_MAX - 1.0);
      printf("finding 1\n");
      printf("%f\t%f\t%f\n", r_new[0], r_new[1], r_new[2]);
    
      }while (overlap( Ply, Ply->firstAtomID, r_new, r ) != 0);

      placeAtom( Ply->firstAtomID, r, r_new );

      for( int i=Ply->firstAtomID+1; i<(Ply->firstAtomID + Ply->numAtoms); ++i ){
        do {
          printf("finding 2\n");

          phi = MD_TWOPI * (double)rand()/(double)RAND_MAX;
          theta=0;
          if( DIM==3 ) {
            dr[2] = MD_q0*( 2.0*(double)rand()/(double)RAND_MAX - 1.0);
            theta = asin(dr[2]/MD_q0);
          }
          else printf("DIM must be 2 or 3 for SAW.\n");
          dr[0] = MD_q0*cos(theta)*cos(phi);
          dr[1] = MD_q0*cos(theta)*sin(phi);
          FOR_ALL_K r_new[k] = r[(i-1)*DIM + k] + dr[k];
          printf("%f\t%f\t%f\n", r_new[0], r_new[1], r_new[2]);
        }while (overlap( Ply,  i, r_new, r ) != 0);

        placeAtom( i, r, r_new );
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
    if( sqrt(r2)<1.0*MD_q0 ) result++;

  }
  return result;

}





