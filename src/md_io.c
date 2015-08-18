#include "hdr/md_io.h"


void printAtomR( int atom, double* r ){
  /* Print a single atom's position to console */

  FOR_ALL_K printf("%.9f\t", r[DIM*atom + k]);

  printf("\n");
}

void saveRtofile( double* r, int N_tot, FILE* fp ){
  /*  Save data to a file easy for MATLAB */

  int i;
  double *Q;

  for( i=0; i<N_tot; i++ ){
  
    Q = r + DIM*i;

    //fprintf( fp, "C\t%lf\t%lf\t%lf\n", Q[0], Q[1], Q[2] );
    fprintf( fp, "\t%lf\t%lf\t%lf\n", Q[0], Q[1], Q[2] );
  
  }
  fprintf( fp, "\n" );


}

void saveXYZtofile( Polymers* PlyList, int numPolymers, double* r, int N_tot, FILE* fp ){
  /*  XYZ file format for importing trajectories into VMD
      NOTE Currently using Polymer ID as the molecule name in VMD, this might be messing it up */

  int i, numAtoms, atomID, id;
  double *Q;

  fprintf(fp,"%d\n", N_tot);
  fprintf(fp,"Obligatory comment\n");
  for( i=0; i<numPolymers; i++ ){
    
    id = (PlyList+i)->firstAtomID;
    numAtoms = (PlyList+i)->numAtoms;
    for (int j = 0; j < numAtoms; ++j)
    { 
      atomID = id + j;
      Q = r + DIM*atomID;
      fprintf( fp, "%d\t%lf\t%lf\t%lf\n", i, Q[0], Q[1], Q[2] );
    }

  
  }


}
