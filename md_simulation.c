#include "hdr/md_simulation.h"



Results simulation( Polymers* Ply, Params parameters, double* r_init, int numPolymers, int simnum ){

	int total_time = parameters.total_time;
	int HYDRO = parameters.HYDRO;
	double temperature = parameters.temperature;
	double mean_Rg = 0.0;
	double mean_r_com = 0.0;
	double *r, *f, *dw, *D, *B, *r_ij;
	int N_tot=0;
	FILE *fp, *fp2;
	double t0, t1, t2, t1000;

	for (int i = 0; i < numPolymers; ++i)
	{
		N_tot += Ply[i].numAtoms;
	}
	
	char filenameR[sizeof "dat/R_XXX_XX.dat"]; 
	sprintf(filenameR, "dat/R_%03d_%02d.dat", N_tot, simnum);

	fp = fopen(filenameR,"w");
	fp2 = fopen("dat/stats.dat","w"); 
	
	srand(time(NULL));
	//printf("About to allocate...\n");

	allocAll( &r, &f, &dw, &D, &B, &r_ij, N_tot );

	//printf("about to initialise r\n");
	for (int i = 0; i < DIM*N_tot; ++i)
	{
		//printf("r_init: %lf\t", r_init[i]);
		r[i] = r_init[i];
		//printf("r: %lf\n", r[i]);	
		//fflush(stdout);

	}
	//printf("about to run simulation...\n");

	t0 = clock();
	t1 = t0;
	t2 = t0;
	for (int t = 0; t < total_time; ++t)
	{
		//printf("%d\n",t);



		update( Ply, numPolymers, r, f, dw, r_ij, D, B, N_tot, HYDRO );
		if( t % 1000 == 0 ){
			
			t1 = clock();

			t1000 = t1 - t2;
			printf("N: %03d\tSim: %02d\tElapsed: %.3f\tRemaining: %.3f\n", N_tot, simnum, (t1-t0)/1000000, (total_time/1000 - t/1000)*t1000/1000000);
			
			t2 = t1;

			//printf("%d\n", t);

			//double* r_com = malloc( 3*sizeof(double) );
			//double Rg = 0.0;
			//FOR_ALL_K r_com[k] = 0.0;
			//calculate_CoM( Ply[0], r, r_com );
			//Rg = calculate_Rg(  Ply[0], r, r_com );
			saveRtofile( r, N_tot, fp );
			//saveStatstofile( t, r_com, Rg, fp2 );
			
			//mean_Rg += (1000/(double)total_time)*Rg;
			//mean_r_com += (1000.0/(double)total_time)*( sqrt(r_com[0]*r_com[0] + r_com[1]*r_com[1]) + r_com[2]*r_com[2]);

			//free(r_com);
		//printf("%lf\t%lf\n", (double)t/(double)total_time, Rg);

		} 
	}
	printf("got here at least.\n");
	Results new_results = { .mean_Rg = mean_Rg, .mean_r_com = mean_r_com };
	printf("and here.\n");
	free(r);
	free(f);
	free(dw);
	free(D);
	free(B);
	free(r_ij);
	printf("freed.\n");
	fclose(fp);
	fclose(fp2);

	return new_results;
}

void update( Polymers* Ply, int numPolymers, double* r, double* f, double* dw, double* r_ij, double* D, double* B, int N_tot, int HYDRO ){
	
	double a1, a2;
	int failure = 0;
	computeForces( Ply, numPolymers, r, f, r_ij, N_tot );
	calc_dw(dw, N_tot);
	static int update_t=0;
	a1 = dt/zeta;
	a2 = sqrt( 2.0*dt );
	//a2 = 0.0;
	
	if (HYDRO)
	{
		/* Need to calculate D, and B=D^1/2 s.t.
		 		dr = D.f dt + B.dw
		*/
		
		double* Df =  malloc( DIM*N_tot * sizeof(double) );
		double* Bdw = malloc( DIM*N_tot * sizeof(double) );
		for (int i = 0; i < DIM*N_tot; ++i)
		{
			Df[i]  = 0.0;
			Bdw[i] = 0.0;
		}
		calcD( r, r_ij, D, N_tot, ROTNE ); // OSEEN or ROTNE
		mat_multiply_A_x( D, DIM*N_tot, f, Df );

		failure = calcB( D, DIM*N_tot, B );

		if(failure){

				FILE *cholerror;
				char filename[sizeof "dat/BD100000.txt"];

    			sprintf(filename, "dat/BD%06d.txt", update_t);
			

				cholerror = fopen( filename, "w");
				//printf("%f\t%f\t%f\t%f\t%f\n", f[n], Df[n], dw[n], Bdw[n], dxtmp);
				
				for (int xx = 0; xx < DIM*N_tot; ++xx)
				{
					fprintf(cholerror, "%f\t", r[xx]);
				}
				fprintf(cholerror, "\n\n");

				for (int xx = 0; xx < DIM*N_tot; ++xx)
				{
					fprintf(cholerror, "%f\t", f[xx]);
				}
				fprintf(cholerror, "\n\n");

				for( int xx=0; xx<DIM*N_tot; xx++ ){

					for (int yy = 0; yy < DIM*N_tot; ++yy)
					{
						fprintf(cholerror, "%f\t", D[squ(xx,yy,DIM*N_tot)]);
					}
					fprintf(cholerror, "\n");
				}
				fprintf(cholerror, "\n");
				for( int xx=0; xx<DIM*N_tot; xx++ ){

					for (int yy = 0; yy < DIM*N_tot; ++yy)
					{
						fprintf(cholerror, "%f\t", B[squ(xx,yy,DIM*N_tot)]);
					}
					fprintf(cholerror, "\n");
				}
				fclose(cholerror);

				exit(0);

		}

		mat_multiply_A_x( B, DIM*N_tot, dw, Bdw );
		
		double dxtmp=0;
		for (int n = 0; n < DIM*N_tot; ++n)
		{
			
			dxtmp = a1*Df[n] + a2*Bdw[n]; // dont forget prefactors
			if(dxtmp > 1.0) {
				/*FILE *cholerror;
				char filename[sizeof "dat/BD100.txt"];

    			sprintf(filename, "dat/BD%03d.txt", n%100);
			

				cholerror = fopen( filename, "w");
				printf("%f\t%f\t%f\t%f\t%f\n", f[n], Df[n], dw[n], Bdw[n], dxtmp);
				
				for( int xx=0; xx<DIM*N_tot; xx++ ){

					for (int yy = 0; yy < DIM*N_tot; ++yy)
					{
						fprintf(cholerror, "%f\t", D[squ(xx,yy,DIM*N_tot)]);
					}
					fprintf(cholerror, "\n");
				}
				fprintf(cholerror, "\n");
				for( int xx=0; xx<DIM*N_tot; xx++ ){

					for (int yy = 0; yy < DIM*N_tot; ++yy)
					{
						fprintf(cholerror, "%f\t", B[squ(xx,yy,DIM*N_tot)]);
					}
					fprintf(cholerror, "\n");
				}
				fclose(cholerror);*/



			}
			r[n]+=dxtmp;

		}

		free(Df);
		free(Bdw);
		update_t ++;

	}
	else{
		
		/* D is diag(1,1,1,1...) so dr becomes:
		 		dr = f dt + dw
		*/


		for (int n = 0; n < DIM*N_tot; ++n)
		{
			//printf("%lf\t%lf\n", a1*f[n], a2*dw[n]);
			//printf("%f\n", f[n]);
			r[n] += a1*f[n] + a2*dw[n]; 
			//printf("%f\t%f\n", a1*f[n], a2*dw[n]);
		}
	}

	for (int n = 0; n < DIM*N_tot; ++n)
	{
		f[n] = 0.0;
	}

}

void calculate_CoM( Polymers Ply, double* r, double* r_com ){

	/* ONLY WORKS FOR DIM=3 */
	int N = Ply.numAtoms;
	int n = Ply.firstAtomID;
	double mx, my, mz; // Mean x, y, z
	mx=0.0;
	my=0.0;
	mz=0.0;
	for (int i = 0; i < N; ++i)
	{
		mx += (1.0/(double)N) * r[ DIM*(n+i)    ];
		my += (1.0/(double)N) * r[ DIM*(n+i) + 1];
		mz += (1.0/(double)N) * r[ DIM*(n+i) + 2];
	}

	r_com[0] = mx;
	r_com[1] = my;
	r_com[2] = mz;

}

double calculate_Rg( Polymers Ply, double* r, double* r_com ){

	/* ONLY WORKS FOR DIM=3 */
	int N = Ply.numAtoms;
	int n = Ply.firstAtomID;
	double Rg2 = 0.0;
	double r_n[3];

	for (int i = 0; i < N; ++i)
	{
		FOR_ALL_K r_n[k] = r[3*(n+i) + k] - r_com[k];
		Rg2 += (1.0/(double)N)*VSqr( r_n );
	}

	return sqrt(Rg2);


}