#include "hdr/md_simulation.h"

void simulation( Polymers* Ply, Params parameters, double* r_init, int numPolymers, int simnum, int poly1, int polyn  ){

	double *r, *f, *dw, *D, *B, *Df, *Bdw, *r_ij;
	FILE *fp;


	int total_time = parameters.total_time;
	int hydro = parameters.hydro;
	double temperature = parameters.temperature;

	int tsave = 10000; // Output/save every tsave steps
	int N_tot=0;
	double t0, t1, t2, t1000; 

	/* Count the total # of monomers in the simulation */
	for (int i = 0; i < numPolymers; ++i) N_tot += Ply[i].numAtoms;
	//printf("---N_tot = %d\n", N_tot);
	char filenameR[sizeof "dat3/R_XXX_XX.dat"]; 
	sprintf(filenameR, "dat3/R_%03d_%02d.xyz", N_tot, simnum);
	fp = fopen(filenameR,"w");
	
	time_t seed = time(NULL);
	srand( seed );
	
	/* Allocate the memory for the main vectors/matrices */
	allocAll( &r, &f, &dw, &D, &B, &Df, &Bdw, &r_ij, N_tot );

	/* Initialise the positions with previously determined initial locations */
	for (int i = 0; i < DIM*N_tot; ++i) {
		r[i] = r_init[i];
	}
	t0 = clock();
	t1 = t0;
	t2 = t0;
	for (int t = 0; t < total_time; ++t)
	{
		/*if ( t > total_time-1000)
		{
			hydro=1;
		}*/
		/* all the action happens in here */
		update( Ply, numPolymers, r, f, dw, r_ij, D, B, Df, Bdw, N_tot, hydro, (double)t*MD_dt, poly1, polyn );
		if( t == 10000 ){

			//flowfield( N_tot, r, f );
		}
		if( t % tsave == 0 ){
			
			t1 = clock();
			t1000 = t1 - t2;
			printf("N: %03d\tSim: %02d\tElapsed: %.3f\tRemaining: %.3f\n", N_tot, simnum, (t1-t0)/1000000, (total_time/tsave - t/tsave)*t1000/1000000);
			t2 = t1;
			meanforce( r, f, Ply[0].numAtoms );

			saveXYZtofile( Ply, numPolymers, r, N_tot, fp );

		} 
	}

	free(r);
	free(f);
	free(dw);
	free(D);
	free(B);
	free(Df);
	free(Bdw);
	free(r_ij);
	fclose(fp);

}

void update( Polymers* Ply, int numPolymers, double* r, double* f, double* dw, double* r_ij, double* D, double* B, double* Df, double* Bdw, int N_tot, int hydro, double t, int poly1, int polyn  ){
	
	/* computeForces() increments f so f needs to be reset */
	for (int n = 0; n < DIM*N_tot; ++n)
	{
		f[n] = 0.0;
		Df[n]  = 0.0;
		Bdw[n] = 0.0;
	}

	/* These are prefactors to D*f and B*dw respectively. Consider making global if zeta never changes */
	double a1 = MD_dt/MD_zeta, a2 = sqrt( 2.0*MD_dt );
	static int update_t=0;
	int calcB_failure = 0;
	int index;

	/* Bonds, soft-sphere repulsions, and velocity-perscribed forces calculated here */
	computeForces( Ply, numPolymers, r, f, r_ij, N_tot, t );

	/* Unit noise calculated here (zero for velocity-perscribed monomers) */
	calc_dw(Ply, numPolymers, dw, N_tot);
	

	if (hydro)
	{
		
		if( (int)(t/MD_dt) % 1000 == 0 ){
			calcD( r, r_ij, D, N_tot, ROTNE ); // NOHI, OSEEN, ROTNE2, ROTNE3 (ROTNE)
			//calcB_failure = calcB( D, B, 0, N_tot );
			calcB_failure = calcB( D, DIM*N_tot, B );
		}
		//mat_multiply_sym_A_x( D, DIM*N_tot, f, Df );
		//mat_multiply_sym_A_x( B, DIM*N_tot, dw, Bdw );
		mat_multiply2( D, f, Df, B, dw, Bdw, DIM*N_tot );
		
		double dxtmp=0;
		for (int p = 0; p < numPolymers; ++p)
		{


			switch( (Ply+p)->perscription ){
				case NORMAL:
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
						r[index]+= a1*Df[index];
					}						
					break;

			}
			
		}
		

		update_t++;

	}
	else{
		
		/* D is diag(1,1,1,1...) so dr becomes:
		 		dr = f MD_dt + dw
		*/

		for (int n = 0; n < DIM*N_tot; ++n) r[n] += a1*f[n] + a2*dw[n]; 

	}



}

void flowfield( int numAtoms, double* r, double* f ){

	int grid_N, grid_M, grid_L;
	FILE* fout = fopen("vfield.txt","w");
	grid_N = 50;
	grid_M = 50;
	grid_L = 150;
	double grid_dx = 1;

	double rx, ry, rz, Rx, Ry, Rz, rr, C1,C2,C3, rf;
	double* vx = (double*)malloc( grid_N*grid_M*sizeof(double) );
	double* vy = (double*)malloc( grid_N*grid_M*sizeof(double) );
	double* vz = (double*)malloc( grid_N*grid_M*sizeof(double) );
	double a = MD_q0/2.0;
	for (int i = 0; i < grid_M*grid_N; ++i)
	{
		vx[i] = 0.0;
		vy[i] = 0.0;
		vz[i] = 0.0;
	}


	rx = -(grid_dx*(double)grid_N)/2;
	ry = -(grid_dx*(double)grid_M)/2;
	rz = -(grid_dx*(double)grid_L)/2;

	printf("\n");
	for (int ix = 0; ix < grid_N; ++ix)
	{			
		
		for (int iy = 0; iy < grid_M; ++iy)
		{
			for (int iz = 0; iz < grid_L; ++iz)
			{
			
			
				for (int n = 0; n < numAtoms; ++n)
				{
					Rx = r[DIM*n];
					Ry = r[DIM*n+1];
					Rz = r[DIM*n+2];

					rr =  (rx-Rx)*(rx-Rx) + (ry-Ry)*(ry-Ry) + (rz-Rz)*(rz-Rz);
					C1 = 1 + 2.0*a*a/(3.0*rr);
					C2 = 1 - 2.0*a*a/rr;
					C3 = 1 - 9.0*sqrt(rr)/(32.0*a);

					rf = (rx - Rx)*f[DIM*n  ] + (ry - Ry)*f[DIM*n+1] + (rz - Rz)*f[DIM*n+2];

					if( rr >= 4.0*a*a ) {
						vx[ ix*grid_M + iy ] += (3.0*a/(4.0*sqrt(rr)) )*( C1*f[DIM*n  ] + C2*rf*(rx-Rx)/rr );				
						vy[ ix*grid_M + iy ] += (3.0*a/(4.0*sqrt(rr)) )*( C1*f[DIM*n+1] + C2*rf*(ry-Ry)/rr );
						vz[ ix*grid_M + iy ] += (3.0*a/(4.0*sqrt(rr)) )*( C1*f[DIM*n+2] + C2*rf*(rz-Rz)/rr );
					}
					else{
						vx[ ix*grid_M + iy ] += C3*f[DIM*n  ] + (3.0/(32.0*a))*rf*(rx-Rx)/sqrt(rr);
						vy[ ix*grid_M + iy ] += C3*f[DIM*n+1] + (3.0/(32.0*a))*rf*(ry-Ry)/sqrt(rr);
						vz[ ix*grid_M + iy ] += C3*f[DIM*n+2] + (3.0/(32.0*a))*rf*(rz-Rz)/sqrt(rr);
				
					}
					
				
				}
				fprintf(fout, "%.3f,%.3f,%.3f\n", vx[ ix*grid_M + iy ],vy[ ix*grid_M + iy ],vz[ ix*grid_M + iy ]);
				rz+=grid_dx;
			}

			ry+=grid_dx;
			rz = -(grid_dx*(double)grid_L)/2;
			//printf("%lf ", vx[ ix*grid_M + iy ]);
		}
		rx+=grid_dx;
		ry = -(grid_dx*(double)grid_M)/2;
		//printf("\n");
	}
	
	/*printf("\n----vy------\n");
	for (int ix = 0; ix < grid_N; ++ix)
	{			
		
		for (int iy = 0; iy < grid_M; ++iy)
		{
			printf("%.3f ", vy[ ix*grid_M + iy ]);
		}
		printf("\n");
	}

	printf("\n----vz------\n");
	for (int ix = 0; ix < grid_N; ++ix)
	{			
		
		for (int iy = 0; iy < grid_M; ++iy)
		{
			printf("%.3f ", vz[ ix*grid_M + iy ]);
		}
		printf("\n");
	}*/

	free(vx);
	free(vy);
	free(vz);
	fclose(fout);

}

void meanforce( double* r, double* f, int numAtoms ){

double mean_radius = 0.0;
double torque_z = 0.0;
double force_z = 0.0;

for (int n = 0; n < numAtoms; ++n)
{	
	mean_radius += sqrt( r[DIM*n]*r[DIM*n] + r[DIM*n+1]*r[DIM*n+1] );
	torque_z += f[DIM*n]*r[DIM*n+1] - f[DIM*n+1]*r[DIM*n];
	force_z += f[DIM*n+2]; 
}

printf("%f\t%f\t%f\n", mean_radius/numAtoms, torque_z, force_z);

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
