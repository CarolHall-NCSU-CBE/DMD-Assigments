//============================================================================
// Name        : squarewell.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "swell.h"

int main() {
	srand((unsigned int)time(NULL));
	//srand48((unsigned short)time(NULL));
	//clock_t start = clock();

	double density;
	double temperature;
	int ncoll;

	// dynamically allocate memory for arrays
	double *(rx), *(ry), *(rz);
	double *(vx), *(vy), *(vz);
	double *(tmprx), *(tmpry), *(tmprz);
	int *(colltype);

	rx = (double*)malloc(n*sizeof(double));
	if(rx == NULL){
		printf("Out of memoryr\n");
		exit(0);
	}
	ry = (double*)malloc(n*sizeof(double));
	if(ry == NULL){
		printf("Out of memoryry\n");
		exit(0);
	}
	rz = (double*)malloc(n*sizeof(double));
	if(rz == NULL){
		printf("Out of memoryrz\n");
		exit(0);
	}

	vx = (double*)malloc(n*sizeof(double));
	if(vx == NULL){
		printf("Out of memoryvx\n");
		exit(0);
	}
	vy = (double*)malloc(n*sizeof(double));
	if(vy == NULL){
		printf("Out of memoryvy\n");
		exit(0);
	}
	vz = (double*)malloc(n*sizeof(double));
	if(vz == NULL){
		printf("Out of memoryvz\n");
		exit(0);
	}

	tmprx = (double*)malloc(n*sizeof(double));
	if(tmprx == NULL){
		printf("Out of memoryr\n");
		exit(0);
	}

	tmpry = (double*)malloc(n*sizeof(double));
	if(tmpry == NULL){
		printf("Out of memoryr\n");
		exit(0);
	}

	tmprz = (double*)malloc(n*sizeof(double));
	if(tmprz == NULL){
		printf("Out of memoryr\n");
		exit(0);
	}

	colltype = (int*)malloc(n*sizeof(int));
	if(colltype == NULL){
		printf("Out of memoryr\n");
		exit(0);
	}
	// end of dynamic allocation;

	double e = 0.0;
	unsigned char overlap;

	double *coltime;
	int *partner;
	coltime = (double*)malloc((n+2)*sizeof(double));
	if(coltime == NULL){
		printf("Out of memoryc\n");
		exit(0);
	}
	partner = (int*)malloc(n*sizeof(int));
	if(partner == NULL){
		printf("Out of memoryp\n");
		exit(0);
	}

	double w;
	double acw;
	double tij;
	int hist[maxbin];

	for(int i = 0; i < maxbin; i++)
	{
		hist[i] = 0;
	}
	for(int j = 0; j < n; j++)
	{
		colltype[j] = 0;
	}
	double penergy = 1;

	printf("Enter Reduced Density (N/V) * Sigma ^ 3: ");
	fflush(stdout);
	scanf("%lf", &density);

	printf("Enter Reduced Temperature: ");
	fflush(stdout);
	scanf("%lf", &temperature);

	printf("Enter Number of Collisions Required: ");
	fflush(stdout);
	scanf("%d", &ncoll);

	double sigma = pow(density/((double)n), 1.0/3.0);
	double sigma2 = (1 + 0.5) * sigma;

	if(sigma2 >= 0.5)
	{
		printf("Square well diameter too large, artificats may develop\n");
		exit(0);
	}

	fcc(rx, ry, rz);

	//comvel(temperature, vx, vy, vz);
	comvel(10.0, vx, vy, vz);

	overlap = check(sigma, rx, ry, rz, vx, vy, vz);
	if(overlap == 1)
	{
		printf("Particle overlap in initial configuration");
		// is return 0 how you're supposed to do this?
		return 0;
	}
	int kecntr = 50;

	e = kineticenergy(vx,vy,vz);
	double totale = e + potentialenergy(rx,ry,rz,sigma2,penergy);
	printf("Initial Energy: %f\n", totale);
	double en = e / ((double)n);
	double calctemp = 2.0 * en / 3.0;
	printf("Initial Temperature: %f\n", calctemp);
	en = e /((double)n);
	double enkt = en/calctemp;
	printf("Initial e/nkt: %f\n", enkt);

	for(int i = 0; i < n; i++)
	{
		coltime[i] = timbig;
		partner[i] = n;
	}

	for(int j = 0; j < n; j++)
	{
		uplist(sigma, sigma2, j, rx, ry, rz, vx, vy, vz, &(coltime[j]), &(partner[j]), &(colltype[j]));
	}

	coltime[n] = 5.0;

	double ghosttime = ghostcoeff * fabs(sqrt(-2.0 * log(((rand()%100)+0.5) / 100.0)) * cos(2.0 * pi * (rand()%100)/100.0));
	printf("initialgt: %f\n", ghosttime);
	coltime[n+1] = ghosttime;

	acw = 0.0;

	printf("***Start of Dynamics***\n");
	fflush(stdout);

	//Begining of Main Loop

	int steps = 0;
	double t = 0.0;
	int uplistcntr = 10;
	int burnin = 100000;
	int grcount = 2000000;
	int nextsphere;
	int temppartner;
	const double deltawrite = 0.005;
	int tcountwrite = 0;
	

	FILE *writeout;
	writeout = fopen("runconfigs", "w+");
	fprintf(writeout, "REMARK    deltat: %f\n", deltawrite);

	//print out initial configuration to runconfigs file
	fprintf(writeout, "MODEL%*d\n", 9, tcountwrite);
	for(int i = 0; i < n; i++)
	{
		double rxi = rx[i];
		double ryi = ry[i];
		double rzi = rz[i];

		fprintf(writeout, "ATOM %*d  N", 6, i);
		fprintf(writeout, "%*.3f%*.3f%*.3f\n", 24, rxi,8,ryi,8,rzi);
	}
	fprintf(writeout, "ENDMDL\n");
	tcountwrite += 1;

	int tmpcoll = -1;
	double burnintime = 0.0;

	// File to store temperature recording
	FILE *temprec;
	temprec = fopen("temprec", "w+");

	unsigned int ghoston = 1;

	for(int coll = 0; coll < ncoll; coll++)
	{
		/*if(coll == 2000000)
		{
			printf("Ghost collisions turned off at %d collisions\n", coll);
			ghoston = 0;
		}*/

		//records the temperature before each collision
		e = kineticenergy(vx,vy,vz);
		en = e / ((double)n);
		calctemp = 2.0 * en / 3.0;
		fprintf(temprec, "%d %f\n", coll, calctemp);

		//switches to user defined (lower) temperature after the higher, burnin temperature
		if(coll == burnin)
		{
			printf("Burn in period complete at %d collisions\n", coll);
			comvel(temperature, vx, vy, vz);
			for(int i = 0; i < n; i++)
			{
				coltime[i] = timbig;
				partner[i] = n;
			}

			for(int j = 0; j < n; j++)
			{
				uplist(sigma, sigma2, j, rx, ry, rz, vx, vy, vz, &(coltime[j]), &(partner[j]), &(colltype[j]));
			}

			// reseting variable for the main part of the loop after burnin
			coltime[n] = 5.0;
			e = kineticenergy(vx,vy,vz);
			totale = e + potentialenergy(rx,ry,rz,sigma2,penergy);
			acw = 0;
			burnintime = t;
			ghosttime = ghostcoeff * fabs(sqrt(-2.0 * log(((rand()%100)+0.5) / 100.0)) * cos(2.0 * pi * (rand()%100)/100.0));
			coltime[n+1] = ghosttime;
			printf("New total energy: %f\n", totale);
		}

		tij = timbig;
		for(int k = 0; k < n+1; k++)
		{
			if(coltime[k] < tij && k != (n-1))
			{
				tij = coltime[k];
				nextsphere = k;
			}
		}

		// checks to see if the ghost collision time is smaller than the normal particles' times
		if(coltime[n+1] < tij && coll > burnin)
		{
			tij = coltime[n+1];
			nextsphere = n+1;
		}

		// check to see if there were any negative collisions
		if(tij <= 0)
		{
			printf("\ntij <= 0\n");
			printf("tij: %f\n", tij);
			printf("nextsphere: %d\n", nextsphere);
			printf("coll: %d\n", coll);
			exit(0);
		}

		// attempt to write out to a file
		while((t+tij) >= (tcountwrite*deltawrite))
		{
			double temptime = tcountwrite * deltawrite - t;
			fprintf(writeout, "MODEL%*d\n", 9, tcountwrite);
			for(int i = 0; i < n; i++)
			{
				double rxi = rx[i];
				double ryi = ry[i];
				double rzi = rz[i];

				rxi = rxi + vx[i]*temptime;
				ryi = ryi + vy[i]*temptime;
				rzi = rzi + vz[i]*temptime;
				rxi = rxi - (1.0*(double)((int)(rxi/(double)(1.0) + (0.5*copysign(1.0,rxi)))));
				ryi = ryi - (1.0*(double)((int)(ryi/(double)(1.0) + (0.5*copysign(1.0,ryi)))));
				rzi = rzi - (1.0*(double)((int)(rzi/(double)(1.0) + (0.5*copysign(1.0,rzi)))));
				
				tmprx[i] = rxi;
				tmpry[i] = ryi;
				tmprz[i] = rzi;

				fprintf(writeout, "ATOM %*d  N", 6, i);
				fprintf(writeout, "%*.3f%*.3f%*.3f\n", 24, rxi,8,ryi,8,rzi);
			}
			if(coll > grcount)
			{
				grsort(tmprx, tmpry, tmprz, hist);
				grcount = grcount + 20;
				coltime[n] = t + 5.0;
				steps = steps + 1;
			}
			tmpcoll = coll;

			fprintf(writeout, "ENDMDL\n");
			tcountwrite = tcountwrite + 1;
		}

		/*if(grcount == coll || nextsphere == n)
		{
			grsort(rx, ry, rz, hist);
			grcount = grcount + 20;
			coltime[n] = t + 5.0;
			steps = steps + 1;
		}
		*/
		
		t = t + tij;
		for(int k = 0; k < n; k++)
		{
			coltime[k] = coltime[k] - tij;
			rx[k] = rx[k] + vx[k]*tij;
			ry[k] = ry[k] + vy[k]*tij;
			rz[k] = rz[k] + vz[k]*tij;
			rx[k] = rx[k] - (1.0*(double)((int)(rx[k]/(double)(1.0) + (0.5*copysign(1.0,rx[k])))));
			ry[k] = ry[k] - (1.0*(double)((int)(ry[k]/(double)(1.0) + (0.5*copysign(1.0,ry[k])))));
			rz[k] = rz[k] - (1.0*(double)((int)(rz[k]/(double)(1.0) + (0.5*copysign(1.0,rz[k])))));
		}

		// dont ever use coltime[n]
		//coltime[n] = coltime[n] - tij;
		if(ghoston == 1)
		{
			coltime[n+1] = coltime[n+1] - tij;
		}

		if(nextsphere < n)
		{
			temppartner = partner[nextsphere];
			
			w = bump(sigma, sigma2, penergy, nextsphere, temppartner, rx, ry, rz, vx, vy, vz, &(colltype[nextsphere]));

			acw = acw + w;

			for(int l = 0; l < n; l++)
			{
				if((l == nextsphere) || (partner[l] == nextsphere) || (l == temppartner) || (partner[l] == temppartner))
				{
					uplist(sigma, sigma2, l, rx, ry, rz, vx, vy, vz, &(coltime[l]), &(partner[l]), &(colltype[l]));
				}
			}

			dnlist(sigma, sigma2, nextsphere, rx, ry, rz, vx, vy, vz, coltime, partner, colltype);
			dnlist(sigma, sigma2, temppartner, rx, ry, rz, vx, vy, vz, coltime, partner, colltype);

			// why is this here?
			if((coll + 1) == (uplistcntr - 1))
			{
				for(int i = 0; i < n; i++)
				{
					uplist(sigma, sigma2, i, rx, ry, rz, vx, vy, vz, &(coltime[i]), &(partner[i]), &(colltype[i]));
				}
				uplistcntr = uplistcntr + 10;
			}
		}
		else if(coll > burnin && nextsphere == n+1)
		{
			// is this the right range of values for rndsphere?
			int rndsphere = (int)(ranf() * n);
			//printf("rndsphere: %d\n", rndsphere);
			if(rndsphere == n)
			{
				printf("not the right particle\n");
			}

			ghostcoll(temperature, &(vx[rndsphere]), &(vy[rndsphere]), &(vz[rndsphere]));
	
			ghosttime = ghostcoeff * fabs(sqrt(-2.0 * log(((rand()%100)+0.5) / 100.0)) * cos(2.0 * pi * (rand()%100)/100.0));
			
			// turns off ghost collisions
			/*if(ghoston == 0)
			{
				ghosttime = 10000;
			}*/

			while (ghosttime <= 0)
			{
				printf("Entered\n");
				ghosttime = ghostcoeff * fabs(sqrt(-2.0 * log(((rand()%100)+1) / 100.0)) * cos(2.0 * pi * (rand()%100)/100.0));
			}

			coltime[n+1] = ghosttime;
			for(int l = 0; l < n; l++)
			{
				if((l == rndsphere) || (partner[l] == rndsphere))
				{
					uplist(sigma, sigma2, l, rx, ry, rz, vx, vy, vz, &(coltime[l]), &(partner[l]), &(colltype[l]));
				}
			}
			dnlist(sigma, sigma2, rndsphere, rx, ry, rz, vx, vy, vz, coltime, partner, colltype);
		}

		/*if ((coll + 1) == kecntr)
		{
			// maybe do this with energy() function later?
			for(int i = 0; i < n; i++)
			{
				e = e + pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2);
			}

			e = 0.5 * e;
			en = e / (double)n;
			enkt = en/calctemp;

			kecntr = kecntr + 50;
		}*/

		// Move forward particles by time tij
	}

	fclose(writeout);
	fclose(temprec);

	// End of Main Loop

	printf("***End of Dynamics***\n");
	printf("Final Colliding Pair: %d %d\n", nextsphere, temppartner);

	overlap = check(sigma, rx, ry, rz, vx, vy, vz);
	if(overlap == 1)
	{
		printf("Particle overlap in final configuration");
		// is return 0 how you're supposed to do this?
		return 0;
	}

	//outputs the ending configuration to a file called endconfig
	FILE *endfile;
	endfile = fopen("endconfig", "w+");
	for(int i = 0; i < n; i++)
	{
		fprintf(endfile, "ATOM %*d  N", 6, i);
		fprintf(endfile, "%*.3f%*.3f%*.3f\n", 24, rx[i],8,ry[i],8,rz[i]);
	}

	fclose(endfile);
	//should this be here?
	e = kineticenergy(vx,vy,vz);
	double pvnkt1 = (acw / ((double)n *3.0 * (t - burnintime) * calctemp)) + 1;
	en = e / ((double)n);
	calctemp = 2.0 * en / 3.0;

	totale = e + potentialenergy(rx,ry,rz,sigma2,penergy);
	
	enkt = en / calctemp;
	t = t * sqrt(calctemp) / sigma;
	double rate = ((double)ncoll) / (t - burnintime);
	double tbc = ((double)n)/rate/2.0;
	double frac = pi * density / 6;

	printf("The final n is: %d\n", n);
	printf("The final e is: %f\n", totale);
	printf("The final temp is: %f\n", calctemp);
	printf("The final volume fraction is: %f\n", frac);
	printf("The final acw is: %f\n", acw);
	printf("Final time is: %f\n", t);
	printf("Collision rate is: %f\n", rate);
	printf("Mean collision time: %f\n", tbc);
	printf("Final e/nkt is: %f\n", enkt);
	printf("PV/nkt is: %f\n", pvnkt1);
	printf("The final grcount is: %d\n", grcount);
	fflush(stdout);

	// calculate g(r) from data from grsort
	double grconst = 4.0*pi*(density/(pow(sigma, 3)))/3.0;
	double gr[maxbin];
	double f[maxbin];

	FILE *grout;
	grout = fopen("gr.out", "w+");

	for(int bin = 0; bin < maxbin; bin++)
	{
		double rlower = (double)(bin)*delr;
		double rupper = rlower + delr;
		double nideal = grconst*(pow(rupper,3) - pow(rlower,3));
		gr[bin] = (double)(hist[bin])/(double)steps/double(n)/nideal;
		f[bin] = rlower + delr/2.0;

		fprintf(grout, "%lf %lf\n", (f[bin]/sigma),(gr[bin]));
	}
	fclose(grout);

	FILE *endingout;
	endingout = fopen("output", "w+");
	fprintf(endingout, "The final n is: %d\n", n);
	fprintf(endingout,"The final e is: %lf\n", e);
	fprintf(endingout,"The final temp is: %lf\n", calctemp);
	fprintf(endingout,"The final volume fraction is: %lf\n", frac);
	fprintf(endingout,"The final acw is: %lf\n", acw);
	fprintf(endingout,"Final time is: %lf\n", t);
	fprintf(endingout,"Collision rate is: %lf\n", rate);
	fprintf(endingout,"Mean collision time: %lf\n", tbc);
	fprintf(endingout,"Final e/nkt is: %lf\n", enkt);
	fprintf(endingout,"PV/nkt is: %lf\n", pvnkt1);
	fprintf(endingout,"The final grcount is: %d\n", grcount);
	fclose(endingout);

	//printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	return(0);
}

/* Creates FCC lattice for n linear molecules
 *
 * int n should be an integer of the form 4 * nc^3
 * where n is defined in hsphere.h
 *
 * int nc				number of fcc unit cells in each direction
 * double rroot3		1.0 / sqrt( 3.0 )
 * ex(n), ey(n), ez(n) 	unit vectors giving orientations
 * */
void fcc(double rx[], double ry[], double rz[])
{
	double (*ex), (*ey), (*ez);
	ex = (double*)malloc(n*sizeof(double));
	ey = (double*)malloc(n*sizeof(double));
	ez = (double*)malloc(n*sizeof(double));

	if((ex == NULL) || (ey == NULL) || (ez == NULL)){
		printf("Out of memory\n");
		exit(0);
	}
	double rxtemp[4], rytemp[4], rztemp[4];

	// calculates nc from n in header file
	int nc = ceil(pow(n / 4.0, 1.0/3.0));
	double rroot3 = 0.5773503;

	// calculates side of unit cell
	double cell = 1.0 / nc;
	double cell2 = 0.5 * cell;

	// Sub-lattice A

	rxtemp[0] = 0.0;
	rytemp[0] = 0.0;
	rztemp[0] = 0.0;
	ex[0] = rroot3;
	ey[0] = rroot3;
	ez[0] = rroot3;

	// Sub-lattice B

	rxtemp[1] = cell2;
	rytemp[1] = cell2;
	rztemp[1] = 0.0;
	ex[1] = rroot3;
	ey[1] = -rroot3;
	ez[1] = -rroot3;

	// Sub-lattice C

	rxtemp[2] = 0.0;
	rytemp[2] = cell2;
	rztemp[2] = cell2;
	ex[2] = -rroot3;
	ey[2] = rroot3;
	ez[2] = -rroot3;

	// Sub-lattice D

	rxtemp[3] = cell2;
	rytemp[3] = 0.0;
	rztemp[3] = cell2;
	ex[3] = -rroot3;
	ey[3] = -rroot3;
	ez[3] = rroot3;

	// Construct lattice from unit cell
	int M = 0;

	FILE *pdbfile;
	pdbfile = fopen("initialconfig", "w+");

	for(int iz = 0; iz < nc; iz++)
	{
		for(int iy = 0; iy < nc; iy++)
		{
			for(int ix = 0; ix < nc; ix++)
			{
				for(int iref = 0; iref < 4; iref++)
				{
					int nextpos  = iref + M;

					rx[nextpos] = rxtemp[iref] + cell * (double)(ix);
					ry[nextpos] = rytemp[iref] + cell * (double)(iy);
					rz[nextpos] = rztemp[iref] + cell * (double)(iz);

					rx[nextpos] = rx[nextpos] - 0.5;
					ry[nextpos] = ry[nextpos] - 0.5;
					rz[nextpos] = rz[nextpos] - 0.5;

					// prints out position data to a file called pdb

					fprintf(pdbfile, "ATOM %*d  N", 6, iref+M);
					fprintf(pdbfile, "%*.3f%*.3f%*.3f\n", 24, rx[nextpos],8,ry[nextpos],8,rz[nextpos]);

					ex[nextpos] = ex[iref];
					ey[nextpos] = ey[iref];
					ez[nextpos] = ez[iref];
				}
				M = M + 4;
			}
		}
	}
	fclose(pdbfile);
}

double bump(double sigma, double sigma2, double energy, int nextsphere, int temppartner, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], int *colltype)
{
	double sigsq = pow(sigma, 2.0);
	double sig2sq = pow(sigma2, 2.0);

	double rxij = rx[nextsphere] - rx[temppartner];
	double ryij = ry[nextsphere] - ry[temppartner];
	double rzij = rz[nextsphere] - rz[temppartner];

	rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0) + (0.5*copysign(1.0,rxij)))));
	ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0) + (0.5*copysign(1.0,ryij)))));
	rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0) + (0.5*copysign(1.0,rzij)))));

	double vxij = vx[nextsphere] - vx[temppartner];
	double vyij = vy[nextsphere] - vy[temppartner];
	double vzij = vz[nextsphere] - vz[temppartner];

	double bij = rxij*vxij + ryij*vyij + rzij*vzij;
	double bijsq = pow(bij, 2.0);

	double delvx = 0;
	double delvy = 0;
	double delvz = 0;
	
	double smdist = 5.0/pow(10.0,11.0);
	double smbump = 0;

	// Hard Sphere Collision
	if(*colltype == 1)
	{
		//printf("Hard sphere collision\n");
		double factor = (bij)/sigsq;

		delvx = -factor*rxij;
		delvy = -factor*ryij;
		delvz = -factor*rzij;
		
		smbump = 1;
	}

	// Attractive Collision
	else if(*colltype == 2)
	{
		double compare = 4*sig2sq*energy;
		
		// has enough energy to get out of well
		if(bijsq > compare)
		{
			double parenterm = -sqrt(-4*sig2sq*energy + bijsq) + bij;
			double factor = 1/(2.0*sig2sq) * parenterm;

			//printf("Dissociation\n");
			delvx = -factor*rxij;
			delvy = -factor*ryij;
			delvz = -factor*rzij;			
			
			smbump = 1;
		}

		//doesnt have enough energy to get out of well
		else if(bijsq <= compare)
		{
			double factor = (bij)/sig2sq;

			//printf("Square Well Bounce\n");
			delvx = -factor*rxij;
			delvy = -factor*ryij;
			delvz = -factor*rzij;
			
			smbump = -1;
		}
	}
	else if(*colltype == 3)
	{
		//printf("Capture\n");
		double root = 4.0*sig2sq*energy + bijsq;
		double parenterm = sqrt(root) + bij;
		double factor = 1.0/(2.0*sig2sq) * parenterm;

		delvx = -factor*rxij;
		delvy = -factor*ryij;
		delvz = -factor*rzij;
	
		smbump = -1;		
	}

	vx[nextsphere] = vx[nextsphere] + delvx;
	vx[temppartner] = vx[temppartner] - delvx;
	vy[nextsphere] = vy[nextsphere] + delvy;
	vy[temppartner] = vy[temppartner] - delvy;
	vz[nextsphere] = vz[nextsphere] + delvz;
	vz[temppartner] = vz[temppartner] - delvz;

	vxij = vx[nextsphere] - vx[temppartner];
	vyij = vy[nextsphere] - vy[temppartner];
	vzij = vz[nextsphere] - vz[temppartner];

	double rij = sqrt( pow(rxij,2.0) + pow(ryij,2.0) + pow(rzij,2.0));

	rx[nextsphere] = rx[nextsphere] + smbump*smdist*rxij/rij;
	ry[nextsphere] = ry[nextsphere] + smbump*smdist*ryij/rij;
	rz[nextsphere] = rz[nextsphere] + smbump*smdist*rzij/rij;

	rx[temppartner] = rx[temppartner] - smbump*smdist*rxij/rij;
	ry[temppartner] = ry[temppartner] - smbump*smdist*ryij/rij;
	rz[temppartner] = rz[temppartner] - smbump*smdist*rzij/rij;

	double w = delvx*rxij + delvy*ryij + delvz*rzij;

	return w;
}

void dnlist(double sigma, double sigma2, int j, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], double coltime[], int partner[], int colltype[])
{
	if(j == 0)
	{
		return;
	}

	double sigsq = pow(sigma,2.0);
	double sig2sq = pow(sigma2, 2.0);
	double rxj, ryj, rzj, vxj, vyj, vzj;
	double rxij, ryij, rzij, vxij, vyij, vzij;
	double bij;

	rxj = rx[j];
	ryj = ry[j];
	rzj = rz[j];
	vxj = vx[j];
	vyj = vy[j];
	vzj = vz[j];

	double tij = timbig;
	int tempcolltype;

	for(int i = 0; i < j; i++)
	{
		rxij = rx[i] - rxj;
		ryij = ry[i] - ryj;
		rzij = rz[i] - rzj;
		rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0) + (0.5*copysign(1.0,rxij)))));
		ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0) + (0.5*copysign(1.0,ryij)))));
		rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0) + (0.5*copysign(1.0,rzij)))));

		vxij = vx[i] - vxj;
		vyij = vy[i] - vyj;
		vzij = vz[i] - vzj;
		bij = rxij*vxij + ryij*vyij + rzij*vzij;

		tij = timbig;

		// CASE I
		if(bij < 0.0)
		{
			double rijsq = pow(rxij,2) + pow(ryij,2) + pow(rzij,2);


			// CASE (a)
			if((rijsq - sig2sq) <= 0)
			{
				double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
				double discr = pow(bij,2) - vijsq*(rijsq - sigsq);

				// CASE (1)
				if(discr > 0)
				{
					//is this minus supposed to be here?
					tij = (-bij - sqrt(discr))/vijsq;
					tempcolltype = 1;
				}

				// CASE(2)
				else if(discr <= 0)
				{
					double discr2 = pow(bij,2) - vijsq*(rijsq - sig2sq);
					tij = (-bij + sqrt(discr2))/vijsq;
					tempcolltype = 2;
				}
			}

			// CASE (b)
			else if((rijsq-sig2sq) > 0)
			{
				double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
				double discr2 = pow(bij,2) - vijsq*(rijsq - sig2sq);

				// CASE (1)
				if(discr2 >= 0)
				{
					tij = (-bij - sqrt(discr2))/vijsq;
					tempcolltype = 3;
				}
			}
		}

		// CASE II
		else if(bij > 0.0)
		{
			double rijsq = pow(rxij,2) + pow(ryij,2) + pow(rzij,2);

			// CASE a
			if((rijsq - sig2sq) < 0)
			{
				double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
				double discr2 = pow(bij,2) - vijsq*(rijsq - sig2sq);

				tij = (-bij + sqrt(discr2))/vijsq;
				tempcolltype = 2;
			}

		}

		if(tij < coltime[i])
		{
			coltime[i] = tij;
			partner[i] = j;
			colltype[i] = tempcolltype;
		}
	}

}

void uplist(double sigma, double sigma2, int i, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], double *coltime, int *partner, int *colltype)
{
	if(i == (n-1))
	{
		return;
	}

	double rxi, ryi, rzi, vxi, vyi, vzi;
	double rxij, ryij, rzij, vxij, vyij, vzij;
	double bij;

	double sigsq = pow(sigma, 2.0);
	double sig2sq = pow(sigma2, 2.0);

	*coltime =  timbig;
	rxi = rx[i];
	ryi = ry[i];
	rzi = rz[i];
	vxi = vx[i];
	vyi = vy[i];
	vzi = vz[i];

	double tij = timbig;

	for(int j = i+1; j < n; j++)
	{
		rxij = rxi - rx[j];
		ryij = ryi - ry[j];
		rzij = rzi - rz[j];

		rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0) + (0.5*copysign(1.0,rxij)))));
		ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0) + (0.5*copysign(1.0,ryij)))));
		rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0) + (0.5*copysign(1.0,rzij)))));

		vxij = vxi - vx[j];
		vyij = vyi - vy[j];
		vzij = vzi - vz[j];
		bij = rxij*vxij + ryij*vyij + rzij*vzij;
		// was missing a j on vzij last one

		tij = timbig;
		int tempcolltype;

		// CASE I
		if(bij < 0.0)
		{
			double rijsq = pow(rxij,2) + pow(ryij,2) + pow(rzij,2);

			// CASE (a), equals because after bump they should be considered inside the well
			if((rijsq - sig2sq) <= 0)
			{
				double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
				double discr = pow(bij,2) - vijsq*(rijsq - sigsq);

				// CASE (1)
				if(discr > 0)
				{
					tij = (-bij - sqrt(discr))/vijsq;
					tempcolltype = 1;
				}

				// CASE(2)
				else
				{
					double discr2 = pow(bij,2) - vijsq*(rijsq - sig2sq);
					/*if(discr2 < 0)
					{
						printf("discr2 < 0 in uplist\n");
						exit(0);
					}*/
					tij = (-bij + sqrt(discr2))/vijsq;
					tempcolltype = 2;
				}
			}

			// CASE (b)
			else
			{
				double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
				double discr2 = pow(bij,2) - vijsq*(rijsq - sig2sq);

				// CASE (1)
				if(discr2 >= 0)
				{
					tij = (-bij - sqrt(discr2))/vijsq;
					tempcolltype = 3;
				}
			}
		}

		// CASE II
		else if(bij > 0.0)
		{
			double rijsq = pow(rxij,2) + pow(ryij,2) + pow(rzij,2);

			// CASE a, NOT equal to since want to consider it out of the well after it collides
			if((rijsq - sig2sq) < 0)
			{
				double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
				double discr2 = pow(bij,2) - vijsq*(rijsq - sig2sq);
				tij = (-bij + sqrt(discr2))/vijsq;
				tempcolltype = 2;
			}

		}

		if(tij < *coltime)
		{
			*colltype = tempcolltype;
			*coltime = tij;
			*partner = j;
		}
	}
}

/*
 * Sorts something?
 */
void grsort(double rx[], double ry[], double rz[], int hist[])
{
	for(int i = 0; i < n-1; i++)
	{
		double rxi = rx[i];
		double ryi = ry[i];
		double rzi = rz[i];

		double rxij;
		double ryij;
		double rzij;
		for(int j = i + 1; j < n; j++)
		{
			rxij = rxi - rx[j];
			ryij = ryi - ry[j];
			rzij = rzi - rz[j];

			rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0) + (0.5*copysign(1.0,rxij)))));
			ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0) + (0.5*copysign(1.0,ryij)))));
			rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0) + (0.5*copysign(1.0,rzij)))));

			double rijsq = pow(rxij, 2) + pow(ryij, 2) + pow(rzij, 2);
			double rij = sqrt(rijsq);
			int bin = (int)(rij/delr);

			if(bin < maxbin)
				hist[bin] = hist[bin] + 2;
		}
	}
}

/*
 * Checks to see if there are any overlaps in the configuration
 */
unsigned char check(double sigma, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[])
{
	double sigsq = pow(sigma,2);
	unsigned char overlap = 0;
	double rij;

	for(int i = 0; i < (n - 1); i++)
	{
		double rxi = rx[i];
		double ryi = ry[i];
		double rzi = rz[i];

		double rxij;
		double ryij;
		double rzij;
		for(int j = i + 1; j < n; j++)
		{
			rxij = rxi - rx[j];
			ryij = ryi - ry[j];
			rzij = rzi - rz[j];

			rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0) + (0.5*copysign(1.0,rxij)))));
			ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0) + (0.5*copysign(1.0,ryij)))));
			rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0) + (0.5*copysign(1.0,rzij)))));

			double rijsq = pow(rxij, 2) + pow(ryij, 2) + pow(rzij, 2);

			if (rijsq < sigsq)
			{
				rij = sqrt(rijsq / sigsq);
				if((1.0 - rij) > tol)
				{
					printf("r%d %d: %f, tol: %f\n", i,j,rij, tol);
					overlap = 1;
				}
			}
		}
	}
	return overlap;
}

/*
 * Determines kinetic energy of the structure from velocity data
 */
double kineticenergy(double vx[], double vy[], double vz[])
{
	double e = 0.0;

	for(int i = 0; i < n; i++)
	{
		e = e + pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2);
	}

	e = 0.5*e;
	return e;
}

double potentialenergy(double rx[], double ry[], double rz[], double sigma2, double energy)
{
	double rxij, ryij, rzij;
	double rijsq;
	
	double pe = 0.0;

	for(int i = 0; i <= (n-2); i++)
	{
		for(int j = (i+1); j < n; j++)
		{
			rxij = rx[i] - rx[j];
			ryij = ry[i] - ry[j];
			rzij = rz[i] - rz[j];

			rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0) + (0.5*copysign(1.0,rxij)))));
			ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0) + (0.5*copysign(1.0,ryij)))));
			rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0) + (0.5*copysign(1.0,rzij)))));
			
			rijsq = pow(rxij, 2.0) + pow(ryij, 2.0) + pow(rzij, 2.0);

			if(((rijsq) - pow(sigma2, 2.0)) <= 0)
			{
				pe = pe - energy;
			}
		}
	}

	return pe;			
}
/*
 * Determines translational velocities from Maxwell-Boltzmann distribution
 */
void comvel(double temp, double vx[], double vy[], double vz[])
{
	double rtemp = sqrt(temp);

	double sumx = 0.0;
	double sumy = 0.0;
	double sumz = 0.0;

	for(int i = 0; i < n; i++)
	{
		vx[i] = rtemp * gauss();
		vy[i] = rtemp * gauss();
		vz[i] = rtemp * gauss();

		sumx = sumx + vx[i];
		sumy = sumy + vy[i];
		sumz = sumz + vz[i];
	}

	sumx = sumx / n;
	sumy = sumy / n;
	sumz = sumz / n;

	for(int j = 0; j < n; j++)
	{
		vx[j] = vx[j] - sumx;
		vy[j] = vy[j] - sumy;
		vz[j] = vz[j] - sumz;
	}
}

void ghostcoll(double usrtemp, double *vx, double *vy, double *vz)
{
	double rtemp = sqrt(usrtemp);
	*vx = rtemp * gauss();
	*vy = rtemp * gauss();
	*vz = rtemp * gauss();
}
	

/*
 * Random variate from gaussian distribution
 * Zero mean and unit variance
 */
double gauss()
{
	float a1 = 3.949846138, a3 = 0.252408784, a5 = 0.076542912, a7 = 0.00835596, a9 = 0.029899776;
	double r, r2;
	double sum = 0.0;

	for(int i = 0; i < 12; i++)
	{
		sum = sum + ranf();
	}

	r = (sum - 6.0) / 4.0;
	r2 = r * r;

	double gauss = (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 +a1 ) * r;

	return gauss;

}

/*
 * Change this random number generator, this one is not good
 */
double ranf()
{
	double dra;
	dra = drand48();
	return dra;

}
