#include "hdr/md_simulation.h"

void simulation( Polymers* Ply, Params parameters, double* r_init, int numPolymers, int simnum ){

	double *r, *f, *dw, *D, *B, *r_ij;
	FILE *fp;


	int total_time = parameters.total_time;
	int hydro = parameters.hydro;
	double temperature = parameters.temperature;
	double size = parameters.size;

	int tsave = 100; // Output/save every tsave steps
	int N_tot=0;
	double t0, t1, t2, t1000; 

	/* Count the total # of monomers in the simulation */
	for (int i = 0; i < numPolymers; ++i) N_tot += Ply[i].numAtoms;
	printf("---N_tot = %d\n", N_tot);
	char filenameR[sizeof "dat/R_XXX_XX.dat"]; 
	sprintf(filenameR, "dat/R_%03d_%02d.xyz", N_tot, simnum);
	fp = fopen(filenameR,"w");
	
	srand(time(NULL));
	
	/* Allocate the memory for the main vectors/matrices */
	allocAll( &r, &f, &dw, &D, &B, &r_ij, N_tot );

	/* Initialise the positions with previously determined initial locations */
	for (int i = 0; i < DIM*N_tot; ++i) r[i] = r_init[i];

	t0 = clock();
	t1 = t0;
	t2 = t0;
	for (int t = 0; t < total_time; ++t)
	{

		/* all the action happens in here */
		update( Ply, numPolymers, r, f, dw, r_ij, D, B, N_tot, hydro, size );
		
		if( t % tsave == 0 ){
			
			t1 = clock();
			t1000 = t1 - t2;
			printf("N: %03d\tSim: %02d\tElapsed: %.3f\tRemaining: %.3f\n", N_tot, simnum, (t1-t0)/1000000, (total_time/tsave - t/tsave)*t1000/1000000);
			t2 = t1;

			saveXYZtofile( Ply, numPolymers, r, N_tot, fp );

		} 
	}

	free(r);
	free(f);
	free(dw);
	free(D);
	free(B);
	free(r_ij);
	fclose(fp);

}

void update( Polymers* Ply, int numPolymers, double* r, double* f, double* dw, double* r_ij, double* D, double* B, int N_tot, int hydro, double L ){
	
	/* These are prefactors to D*f and B*dw respectively. Consider making global if zeta never changes */
	double a1 = MD_dt/MD_zeta, a2 = sqrt( 2.0*MD_dt );
	static int update_t=0;
	int calcB_failure = 0;
	int index;

	/* Bonds, soft-sphere repulsions, and velocity-perscribed forces calculated here */
	computeForces( Ply, numPolymers, r, f, r_ij, N_tot );

	/* Unit noise calculated here (zero for velocity-perscribed monomers) */
	calc_dw(Ply, numPolymers, dw, N_tot);
	

	if (hydro)
	{
		/* NOTE maybe faster, if less mem-efficient to update Df and Bdw rather than create, initialise etc */
		double* Df =  malloc( DIM*N_tot * sizeof(double) );
		double* Bdw = malloc( DIM*N_tot * sizeof(double) );
		
		for (int i = 0; i < DIM*N_tot; ++i)
		{
			Df[i]  = 0.0;
			Bdw[i] = 0.0;
		}
		calcD( r, r_ij, D, N_tot, ROTNE ); // NOHI, OSEEN, ROTNE2, ROTNE3 (ROTNE)
		mat_multiply_A_x( D, DIM*N_tot, f, Df );

		calcB_failure = calcB( D, DIM*N_tot, B );
		mat_multiply_A_x( B, DIM*N_tot, dw, Bdw );
		
		double dxtmp=0;
		for (int p = 0; p < numPolymers; ++p)
		{


			
			switch((Ply+p)->perscription){
				case 0:
					for (int n = 0; n < DIM*((Ply+p)->numAtoms); ++n)
					{			
						index = DIM*((Ply+p)->firstAtomID) + n;		
						dxtmp = a1*Df[index] + a2*Bdw[index]; 
						r[index] += dxtmp;
					}
					break;
				case HELIX:
					for (int n = 0; n < DIM*((Ply+p)->numAtoms); ++n)
					{

						index = DIM*((Ply+p)->firstAtomID) + n;		
						r[index]+= a1*f[index];
					}						
					break;

			}
			
		}
		

		free(Df);
		free(Bdw);
		update_t++;

	}
	else{
		
		/* D is diag(1,1,1,1...) so dr becomes:
		 		dr = f MD_dt + dw
		*/

		for (int n = 0; n < DIM*N_tot; ++n) r[n] += a1*f[n] + a2*dw[n]; 

	}

	/* computeForces() increments f so f needs to be reset */
	for (int n = 0; n < DIM*N_tot; ++n)
	{
		f[n] = 0.0;
	}

	double rn;
	for (int n = 0; n < DIM*N_tot; ++n)
	{
		rn = r[n];
		if( rn < L && rn >=0 ){


		}
		else if ( rn >= L)	
		{
			r[n] -= L;
		}
		else{

			r[n] += L;
		}
	}

}

void calculate_CoM( Polymers Ply, double* r, double* r_com ){

	int N = Ply.numAtoms;
	int n = Ply.firstAtomID;
	double mx, my, mz; // Mean x, y, z
	mx=0.0;
	my=0.0;
	mz=0.0;
	for (int i = 0; i < N; ++i)
	{
		mx += 				(1.0/(double)N) * r[ DIM*(n+i)    ];
		my += 				(1.0/(double)N) * r[ DIM*(n+i) + 1];
		if( DIM == 3 ) mz +=(1.0/(double)N) * r[ DIM*(n+i) + 2];
	}

	r_com[0] = mx;
	r_com[1] = my;
	if( DIM == 3 ) r_com[2] = mz;

}

double calculate_Rg( Polymers Ply, double* r, double* r_com ){


	int N = Ply.numAtoms;
	int n = Ply.firstAtomID;
	double Rg2 = 0.0;
	double r_n[DIM];

	for (int i = 0; i < N; ++i)
	{
		FOR_ALL_K r_n[k] = r[DIM*(n+i) + k] - r_com[k];
		FOR_ALL_K Rg2 += (1.0/(double)N)*r_n[k]*r_n[k];
	}

	return sqrt(Rg2);


}