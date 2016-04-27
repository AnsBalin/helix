#include "hdr/md_simulation.h"

void simulation( Polymers* Ply, Params parameters, double* r_init, int numPolymers, int simnum ){

	double *r, *f, *dw, *D, *B, *Df, *Bdw, *r_ij, Rg, sqDisp;
	FILE *fp1, *fp2, *fp3, *fp4;

	double *r_com = malloc( 3*sizeof(double) );

	int total_time = parameters.total_time;
	int hydro = parameters.hydro;
	double temperature = parameters.temperature;

	int tsave = 100; // Output/save every tsave steps
	int N_tot=0;
	double t0, t1, t2, t1000; 

	/* Count the total # of monomers in the simulation */
	for (int i = 0; i < numPolymers; ++i) N_tot += Ply[i].numAtoms;
	//printf("---N_tot = %d\n", N_tot);
	char filenameR[sizeof "flowfield_exp1/R_XXX_XXX.dat"]; 
	sprintf(filenameR, "flowfield_exp1/R_%03d_%03d.xyz", N_tot, simnum);
	
	char filenameF[sizeof "poly_exp1/F_XXX_XXX.dat"]; 
	sprintf(filenameF, "flowfield_exp1/F_%03d_%03d.xyz", N_tot, simnum);

	
	/*char filenameSqDisp[sizeof "poly_exp1/sqDisp_XXX_XXX.dat"]; 
	sprintf(filenameSqDisp, "poly_exp1/sqDisp_%03d_%03d.xyz", N_tot, simnum);
	
	char filenameRg[sizeof "poly_exp1/Rg_XXX_XXX.dat"]; 
	sprintf(filenameRg, "poly_exp1/Rg_%03d_%03d.xyz", N_tot, simnum);
	*/
	fp1 = fopen(filenameR,"w");
	fp2 = fopen(filenameF,"w");
	//fp3 = fopen(filenameSqDisp,"w");
	//fp4 = fopen(filenameRg,"w");

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

	int N, M, L;
	N = 50;
	M = 50;
	L = 150;
	double* vx0 = malloc( sizeof(double)*N*M*L );
	double* vy0 = malloc( sizeof(double)*N*M*L );
	double* vz0 = malloc( sizeof(double)*N*M*L );

	for(int i=0; i<M*L; ++i){

		vx0[i] = 0.0;
		vy0[i] = 0.0;
		vz0[i] = 0.0;
	}

	//calculate_CoM(Ply, r, r_com);
	//shift( r, r_com, N_tot );


	for (int t = 0; t < total_time; ++t)
	{
		/*if ( t > total_time-1000)
		{
			hydro=1;
		}*/
		/* all the action happens in here */
		//printf("%d\n", t);
		update( Ply, numPolymers, r, f, dw, r_ij, D, B, Df, Bdw, N_tot, hydro, (double)t*MD_dt);
		if( t >= 50000 && (t%10==0) ){

			flowfield3( N_tot, r, f, vx0, vy0, vz0, N,M, L );
		}

		//calculate_CoM(Ply, r, r_com);
		//Rg = calculate_Rg(Ply, r, r_com);
		if ( t % tsave == 0){
			//fprintf( fp4, "%.3f\n", Rg );
		}
		//sqDisp = calculate_SqDisplacement( r_com );
		//fprintf( fp3, "%.3f\n", sqDisp );
		if( t % tsave == 0 ){
			
			t1 = clock();
			t1000 = t1 - t2;
			printf("N: %03d\tSim: %02d\tElapsed: %.3f\tRemaining: %.3f\n", N_tot, simnum, (t1-t0)/1000000, (total_time/tsave - t/tsave)*t1000/1000000);
			t2 = t1;
			meanforce( r, f, Ply[0].numAtoms );

			saveXYZtofile( Ply, numPolymers, r, N_tot, fp1 );
			saveXYZtofile( Ply, numPolymers, f, N_tot, fp2 );

		} 
	}

	FILE* flowfield = fopen("vfield.txt","w");
	for(int i=0; i<N*M*L; ++i){
		fprintf(flowfield, "%.3f,%.3f,%.3f\n",  vx0[i], vy0[i], vz0[i]);
	}

	free(vx0);
	free(vy0);
	free(vz0);

	free(r);
	free(f);
	free(dw);
	free(D);
	free(B);
	free(Df);
	free(Bdw);
	free(r_ij);
	fclose(fp1);
	fclose(fp2);
	//fclose(fp3);
	//close(fp4);
	fclose(flowfield);

}

void simulation2( Polymers* Ply, Params parameters, double* r_init, int numPolymers, int simnum, double* sqDisp, double* Rg ){

	double *r, *f, *dw, *D, *B, *Df, *Bdw, *r_ij;
	FILE *fp1, *fp2, *fp3, *fp4;

	double *r_com = malloc( 3*sizeof(double) );

	int total_time = parameters.total_time;
	int hydro = parameters.hydro;
	double temperature = parameters.temperature;

	int tsave = 100; // Output/save every tsave steps
	int N_tot=0;
	double t0, t1, t2, t1000; 

	/* Count the total # of monomers in the simulation */
	for (int i = 0; i < numPolymers; ++i) N_tot += Ply[i].numAtoms;
	//printf("---N_tot = %d\n", N_tot);
	char filenameR[sizeof "flowfield_exp1/R_XXX_XXX.dat"]; 
	sprintf(filenameR, "poly_exp2/R_%03d_%03d.xyz", N_tot, simnum);
	
	char filenameF[sizeof "poly_exp1/F_XXX_XXX.dat"]; 
	sprintf(filenameF, "poly_exp2/F_%03d_%03d.xyz", N_tot, simnum);


	fp1 = fopen(filenameR,"w");
	fp2 = fopen(filenameF,"w");


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


	calculate_CoM(Ply, r, r_com);
	shift( r, r_com, N_tot );


	for (int t = 0; t < total_time; ++t)
	{

		update( Ply, numPolymers, r, f, dw, r_ij, D, B, Df, Bdw, N_tot, hydro, (double)t*MD_dt);


		calculate_CoM(Ply, r, r_com);
		Rg[t] += calculate_Rg(Ply, r, r_com);
		sqDisp[t] += calculate_SqDisplacement( r_com );
		
		if( t % tsave == 0 ){
			
			t1 = clock();
			t1000 = t1 - t2;
			printf("N: %03d\tSim: %02d\tElapsed: %.3f\tRemaining: %.3f\n", N_tot, simnum, (t1-t0)/1000000, (total_time/tsave - t/tsave)*t1000/1000000);
			t2 = t1;
			meanforce( r, f, Ply[0].numAtoms );

			//saveXYZtofile( Ply, numPolymers, r, N_tot, fp1 );
			//saveXYZtofile( Ply, numPolymers, f, N_tot, fp2 );

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
	fclose(fp1);
	fclose(fp2);

}


void update( Polymers* Ply, int numPolymers, double* r, double* f, double* dw, double* r_ij, double* D, double* B, double* Df, double* Bdw, int N_tot, int hydro, double t  ){
	
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
					/*if (r[0] > 5.5){
						for (int n = 0; n < ((Ply+p)->numAtoms); ++n)
						{			
							r[n] -= 10.0;
						}
					}
					else if (r[0] < -5.5){
						for (int n = 0; n < ((Ply+p)->numAtoms); ++n)
						{			
							r[n] += 10.0;
						}
					}
					else if (r[1] > 5.5){
						for (int n = 0; n < ((Ply+p)->numAtoms); ++n)
						{			
							r[n+1] -= 10.0;
						}
					}
					else if (r[1] < -5.5){
						for (int n = 0; n < ((Ply+p)->numAtoms); ++n)
						{			
							r[n+1] += 10.0;
						}
					}
					else if (r[2] > 5.5){
						for (int n = 0; n < ((Ply+p)->numAtoms); ++n)
						{			
							r[n+2] -= 10.0;
						}
					}
					else if (r[2] < -5.5){
						for (int n = 0; n < ((Ply+p)->numAtoms); ++n)
						{			
							r[n+2] += 10.0;
						}
					}*/
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
				
		double dxtmp=0;
		for (int p = 0; p < numPolymers; ++p)
		{


		
					for (int n = 0; n < DIM*((Ply+p)->numAtoms); ++n)
					{			
						index = DIM*((Ply+p)->firstAtomID) + n;		
						dxtmp = a1*f[index] + a2*dw[index]; 
						
						r[index] += dxtmp; //+0.00000001;
					}
					
					/*if (r[DIM*((Ply+p)->firstAtomID)] > 10.5){
						for (int n = 0; n < ((Ply+p)->numAtoms); ++n)
						{			
							r[DIM*((Ply+p)->firstAtomID+n)] -= 20.0;
						}
					}
					else if (r[DIM*((Ply+p)->firstAtomID)] < -10.5){
						for (int n = 0; n < ((Ply+p)->numAtoms); ++n)
						{			
							r[DIM*((Ply+p)->firstAtomID+n)] += 20.0;
						}
					}
					else if (r[DIM*((Ply+p)->firstAtomID)+1] > 10.5){
						for (int n = 0; n < ((Ply+p)->numAtoms); ++n)
						{			
							r[DIM*((Ply+p)->firstAtomID+n)+1] -= 20.0;
						}
					}
					else if (r[DIM*((Ply+p)->firstAtomID) + 1] < -10.5){
						for (int n = 0; n < ((Ply+p)->numAtoms); ++n)
						{			
							r[DIM*((Ply+p)->firstAtomID+n)+1] += 20.0;
						}
					}
					else if (r[DIM*((Ply+p)->firstAtomID)+2] > 10.5){
						for (int n = 0; n < ((Ply+p)->numAtoms); ++n)
						{			
							r[DIM*((Ply+p)->firstAtomID+n)+2] -= 20.0;
						}
					}
					else if (r[DIM*((Ply+p)->firstAtomID)+2] < -10.5){
						for (int n = 0; n < ((Ply+p)->numAtoms); ++n)
						{			
							r[DIM*((Ply+p)->firstAtomID+n)+2] += 20.0;
						}
					}*/
			
		}
		

		update_t++;
		/* D is diag(1,1,1,1...) so dr becomes:
		 		dr = f MD_dt + dw
		*/

		//for (int n = 0; n < DIM*N_tot; ++n) r[n] += a1*f[n] + a2*dw[n]; 

	}



}

void flowfield( int numAtoms, double* r, double* f, double* vx0, double* vy0, double* vz0, int M, int L  ){

	int grid_N, grid_M, grid_L;
	//FILE* fout = fopen("vfield.txt","w");
	grid_N = 50;
	grid_M = 50;
	grid_L = 150;
	double grid_dx = 0.5*MD_q0;

	double rx, ry, rz, Rx, Ry, Rz, rr, C1,C2,C3, rf;
	//double* vx = (double*)malloc( grid_N*grid_M*sizeof(double) );
	//double* vy = (double*)malloc( grid_N*grid_M*sizeof(double) );
	//double* vz = (double*)malloc( grid_N*grid_M*sizeof(double) );
	double a = MD_q0/2.0;
	/*for (int i = 0; i < grid_M*grid_N; ++i)
	{
		vx[i] = 0.0;
		vy[i] = 0.0;
		vz[i] = 0.0;
	}*/


	rx = -(grid_dx*(double)grid_N)/2;
	ry = 0.0;//-(grid_dx*(double)grid_M)/2;
	rz = -(grid_dx*(double)grid_L)/2;

	printf("\n");

	for (int ix = 0; ix < grid_N; ++ix)
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
					vx0[ ix*grid_L + iz ] += (3.0*a/(4.0*sqrt(rr)) )*( C1*f[DIM*n  ] + C2*rf*(rx-Rx)/rr );				
					vy0[ ix*grid_L + iz ] += (3.0*a/(4.0*sqrt(rr)) )*( C1*f[DIM*n+1] + C2*rf*(ry-Ry)/rr );
					vz0[ ix*grid_L + iz ] += (3.0*a/(4.0*sqrt(rr)) )*( C1*f[DIM*n+2] + C2*rf*(rz-Rz)/rr );
					//_vx += (3.0*a/(4.0*sqrt(rr)) )*( C1*f[DIM*n  ] + C2*rf*(rx-Rx)/rr );				
					//_vy += (3.0*a/(4.0*sqrt(rr)) )*( C1*f[DIM*n+1] + C2*rf*(ry-Ry)/rr );
					//_vz += (3.0*a/(4.0*sqrt(rr)) )*( C1*f[DIM*n+2] + C2*rf*(rz-Rz)/rr );
				}
				else{
					vx0[ ix*grid_L + iz ] += C3*f[DIM*n  ] + (3.0/(32.0*a))*rf*(rx-Rx)/sqrt(rr);
					vy0[ ix*grid_L + iz ] += C3*f[DIM*n+1] + (3.0/(32.0*a))*rf*(ry-Ry)/sqrt(rr);
					vz0[ ix*grid_L + iz ] += C3*f[DIM*n+2] + (3.0/(32.0*a))*rf*(rz-Rz)/sqrt(rr);
					//_vx += C3*f[DIM*n  ] + (3.0/(32.0*a))*rf*(rx-Rx)/sqrt(rr);
					//_vy += C3*f[DIM*n+1] + (3.0/(32.0*a))*rf*(ry-Ry)/sqrt(rr);
					//_vz += C3*f[DIM*n+2] + (3.0/(32.0*a))*rf*(rz-Rz)/sqrt(rr);
		
				}
				
			
			}
			//fprintf(fout, "%.3f,%.3f,%.3f\n", vx[ ix*grid_M + iy ],vy[ ix*grid_M + iy ],vz[ ix*grid_M + iy ]);
			//fprintf(fout, "%.3f,%.3f,%.3f\n", _vx, _vy, _vz);
			rz+=grid_dx;
		}

		//ry+=grid_dx;
		rz = -(grid_dx*(double)grid_L)/2;
		//printf("%lf ", vx[ ix*grid_M + iy ]);

		rx+=grid_dx;
		//ry = -(grid_dx*(double)grid_M)/2;
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

	//free(vx);
	//free(vy);
	//free(vz);
	//fclose(fout);

}
void flowfield3( int numAtoms, double* r, double* f, double* vx0, double* vy0, double* vz0, int N, int M, int L  ){

	int grid_N, grid_M, grid_L;
	grid_N = N;
	grid_M = M;
	grid_L = L;
	double grid_dx = 0.5*MD_q0;

	double rx, ry, rz, Rx, Ry, Rz, rr, C1,C2,C3, rf;

	double a = MD_q0/2.0;


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
						vx0[ (ix*grid_M + iy)*grid_L + iz ] += (3.0*a/(4.0*sqrt(rr)) )*( C1*f[DIM*n  ] + C2*rf*(rx-Rx)/rr );				
						vy0[ (ix*grid_M + iy)*grid_L + iz ] += (3.0*a/(4.0*sqrt(rr)) )*( C1*f[DIM*n+1] + C2*rf*(ry-Ry)/rr );
						vz0[ (ix*grid_M + iy)*grid_L + iz ] += (3.0*a/(4.0*sqrt(rr)) )*( C1*f[DIM*n+2] + C2*rf*(rz-Rz)/rr );

					}
					else{
						vx0[ (ix*grid_M + iy)*grid_L + iz ] += C3*f[DIM*n  ] + (3.0/(32.0*a))*rf*(rx-Rx)/sqrt(rr);
						vy0[ (ix*grid_M + iy)*grid_L + iz ] += C3*f[DIM*n+1] + (3.0/(32.0*a))*rf*(ry-Ry)/sqrt(rr);
						vz0[ (ix*grid_M + iy)*grid_L + iz ] += C3*f[DIM*n+2] + (3.0/(32.0*a))*rf*(rz-Rz)/sqrt(rr);

					}
					
				}
				rz+=grid_dx;
			}

			rz = -(grid_dx*(double)grid_L)/2;
			ry+=grid_dx;
		}
		ry = -(grid_dx*(double)grid_M)/2;
		rx+=grid_dx;


	}

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

void calculate_CoM( Polymers* Ply, double* r, double* r_com ){

	int N = Ply->numAtoms;
	int n = Ply->firstAtomID;
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

double calculate_Rg( Polymers* Ply, double* r, double* r_com ){


	int N = Ply->numAtoms;
	int n = Ply->firstAtomID;
	double Rg2 = 0.0;
	double r_n[DIM];

	for (int i = 0; i < N; ++i)
	{
		FOR_ALL_K r_n[k] = r[DIM*(n+i) + k] - r_com[k];
		FOR_ALL_K Rg2 += (1.0/(double)N)*r_n[k]*r_n[k];
	}

	return sqrt(Rg2);


}

void shift( double* r, double* r_com, int N ){

	for(int i=0; i<N; ++i){

		for(int k=0; k<3; ++k){
			r[3*i+k] -= r_com[k];
		}

	}
}

double calculate_SqDisplacement( double* r_com ){

	return r_com[0]*r_com[0]+r_com[1]*r_com[1]+r_com[2]*r_com[2];
}
