#include "hdr/md_io.h"


void printAtomR( int atom, double* r ){
  
  FOR_ALL_K printf("%.9f\t", r[DIM*atom + k]);

  printf("\n");
}

void saveRtofile( double* r, int N_tot, FILE* fp ){
  
  int i;
  double *Q;

  //fprintf(fp,"%d\n", N_tot);
  //fprintf(fp,"Comment haha\n");
  for( i=0; i<N_tot; i++ ){
  
    Q = r + DIM*i;

    //fprintf( fp, "C\t%lf\t%lf\t%lf\n", Q[0], Q[1], Q[2] );
    fprintf( fp, "\t%lf\t%lf\t%lf\n", Q[0], Q[1], Q[2] );
  
  }
  fprintf( fp, "\n" );


}

void saveStatstofile( double t, double* r_com, double Rg, FILE* fp ){

	double r_com_mag = VSqr( r_com );
	fprintf( fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, r_com[0], r_com[1], r_com[2], r_com_mag, Rg );

}