//============================================================================
// Name        : hardsphere.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "hsphere.h"

int main() {
	srand((unsigned int)time(NULL));
	//srand48((unsigned int)time(NULL));

	//srand(401);
	//char title[80];
	double density;
	double temperature;
	int ncoll;
	double rx[n], ry[n], rz[n];
	double vx[n], vy[n], vz[n];
	double e = 0.0;
	unsigned char overlap;
	double coltime[n+1];
	int partner[n];
	double w;
	double acw;
	double tij;
	int hist[maxbin];
	double frac;

	for(int i = 0; i < maxbin; i++)
	{
		hist[i] = 0;
	}

	//printf("Enter Reduced Density (N/V) * Sigma ^ 3: ");
	//fflush(stdout);
	//scanf("%lf", &density);

	printf("Enter Packing Fraction (V/V0): ");
	fflush(stdout);
	scanf("%lf", &frac);

	density = (6 * frac) / pi;

	printf("Enter Number of Collisions Required: ");
	fflush(stdout);
	scanf("%d", &ncoll);

	printf("Enter Reduced Temperature: ");
	fflush(stdout);
	scanf("%lf", &temperature);

	//temperature = floor(double(ncoll)/100000.0);
	cout << "Initial Temp: " << temperature << endl;

	double sigma = pow(density/((double)n), 1.0/3.0);
	//printf("%lf\n", sigma);

	fcc(rx, ry, rz);
	comvel(temperature, vx, vy, vz);

	overlap = check(sigma, rx, ry, rz, vx, vy, vz);
	if(overlap == 1)
	{
		printf("Particle overlap in initial configuration");
		// is return 0 how you're supposed to do this?
		return 0;
	}
	int kecntr = 50;

	e = energy(vx,vy,vz);
	double en = e / ((double)n);
	double calctemp = 2.0 * en / 3.0;
	printf("Initial Temperature: %f\n", calctemp);
	printf("Initial e: %f\n", e);
	double enkt = en/calctemp;
	printf("Initial e/nkt: %f\n", enkt);

	for(int i = 0; i < n; i++)
	{
		coltime[i] = timbig;
		partner[i] = n;
	}

	for(int j = 0; j < n; j++)
	{
		uplist(sigma, j, rx, ry, rz, vx, vy, vz, &(coltime[j]), &(partner[j]));
	}

	coltime[n] = timbig;
	acw = 0.0;

	//deletes old runconfigs file with earlier data
	if(remove("runconfigs") != 0)
		perror("Could not delete runconfigs");

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

	printf("***Start of Dynamics***\n");
	fflush(stdout);


	//Begining of Main Loop

	int steps = 0;
	double t = 0.0;
	int uplistcntr = 10;
	int grcount = 20;
	int nextsphere;
	int temppartner;

	for(int coll = 0; coll < ncoll; coll++)
	{
		tij = timbig;
		for(int k = 0; k < n+1; k++)
		{
			if(coltime[k] < tij)
			{
				tij = coltime[k];
				nextsphere = k;
			}
		}

		if(coltime[nextsphere] < -0.000001)
		{
			printf("negative collision\n");
			printf("coll: %d\n", coll);
			printf("nextsphere: %d, partner: %d\n", nextsphere, partner[nextsphere]);
			printf("%f\n", coltime[nextsphere]);
			exit(0);
		}

		if(coll % 100000 == 0)
		{
			printf("coll: %d\n", coll);
		}
		// attempt to write out position data to a file
		while((t + tij) >= (tcountwrite*deltawrite))
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
				rxi = rxi - rint(rxi);
				ryi = ryi - rint(ryi);
				rzi = rzi - rint(rzi);

				fprintf(writeout, "ATOM %*d  N", 6, i);
				fprintf(writeout, "%*.3f%*.3f%*.3f\n", 24, rxi,8,ryi,8,rzi);
			}

			fprintf(writeout, "ENDMDL\n");
			tcountwrite = tcountwrite + 1;
		}

		if((grcount == coll) || (nextsphere == n))
		{
			if(grcount <= coll)
			{
				grsort(rx, ry, rz, hist);
				grcount = grcount + 20;
			}

			//coltime[n] = ((1.0-sigma)/100.0) / maxvelocity;
			steps = steps + 1;
			if(nextsphere == n)
			{
				printf("nxtsphere was n, coll: %d\n", coll);
				coltime[n] = timbig;

				t = t + tij;
				for(int k = 0; k < n; k++)
				{
					coltime[k] = coltime[k] - tij;
					rx[k] = rx[k] + vx[k]*tij;
					ry[k] = ry[k] + vy[k]*tij;
					rz[k] = rz[k] + vz[k]*tij;
					rx[k] = rx[k] - rint(rx[k]);
					ry[k] = ry[k] - rint(ry[k]);
					rz[k] = rz[k] - rint(rz[k]);
				}

				for(int i = 0; i < n; i++)
				{
					uplist(sigma, i, rx, ry, rz, vx, vy, vz, &(coltime[i]), &(partner[i]));
				}
					//printf("tij: %f\n", tij);
			}

		}
		else
		{
			temppartner = partner[nextsphere];
			t = t + tij;
			for(int k = 0; k < n; k++)
			{
				coltime[k] = coltime[k] - tij;
				rx[k] = rx[k] + vx[k]*tij;
				ry[k] = ry[k] + vy[k]*tij;
				rz[k] = rz[k] + vz[k]*tij;
				rx[k] = rx[k] - rint(rx[k]);
				ry[k] = ry[k] - rint(ry[k]);
				rz[k] = rz[k] - rint(rz[k]);
			}

			coltime[n] = coltime[n] - tij;

			w = bump(sigma, nextsphere, temppartner, rx, ry, rz, vx, vy, vz);

			acw = acw + w;

			for(int l = 0; l < n; l++)
			{
				if((l == nextsphere) || (partner[l] == nextsphere) || (l == temppartner) || (partner[l] == temppartner))
				{
					uplist(sigma, l, rx, ry, rz, vx, vy, vz, &(coltime[l]), &(partner[l]));
				}
			}

			dnlist(sigma, nextsphere, rx, ry, rz, vx, vy, vz, coltime, partner);
			dnlist(sigma, temppartner, rx, ry, rz, vx, vy, vz, coltime, partner);

			/*if((coll + 1) == (uplistcntr - 1))
			{
				for(int i = 0; i < n; i++)
				{
					uplist(sigma, i, rx, ry, rz, vx, vy, vz, &(coltime[i]), &(partner[i]));
				}
				uplistcntr = uplistcntr + 10;
			}*/
		}

		if ((coll + 1) == kecntr)
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
		}

		// Move forward particles by time tij
	}
	fclose(writeout);

	// End of Main Loop



	printf("\n***End of Dynamics***\n");
	printf("Final Colliding Pair: %d %d\n", nextsphere, temppartner);

	overlap = check(sigma, rx, ry, rz, vx, vy, vz);
	if(overlap == 1)
	{
		printf("Particle overlap in final configuration");
		// is return 0 how you're supposed to do this?
		exit(0);
	}

	//prints out ending configuration to a file called endconfig
	FILE *endfile;
	endfile = fopen("endconfig", "w+");
	for(int i = 0; i < n; i++)
	{
		fprintf(endfile, "ATOM %*d  N", 6, i);
		fprintf(endfile, "%*.3f%*.3f%*.3f\n", 24, rx[i],8,ry[i],8,rz[i]);
	}

	fclose(endfile);
	//should this be here?
	e = energy(vx,vy,vz);
	double pvnkt1 = (acw / ((double)n *3.0 * t * calctemp)) + 1;
	en = e / ((double)n);
	calctemp = 2.0 * en / 3.0;
	enkt = en / calctemp;
	t = t* sqrt(calctemp) / sigma;
	double rate = ((double)ncoll) / t;
	double tbc = ((double)n)/rate/2.0;
	frac = pi * density / 6;

	printf("The final n is: %d\n", n);
	printf("The final e is: %f\n", e);
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
	double grconst = 4.0*pi*(density/(pow(sigma, 3.0)))/3.0;
	double grtemp;
	double ftemp;

	FILE *grout;
	grout = fopen("gr.out", "w+");

	for(int bin = 0; bin < maxbin; bin++)
	{
		double rlower = (double)(bin)*delr;
		double rupper = rlower + delr;
		double nideal = grconst*(pow(rupper,3.0) - pow(rlower,3.0));

		grtemp = ((double)hist[bin])/((double)steps)/((double)n)/nideal;
		ftemp = rlower + delr/2.0;

		fprintf(grout, "%lf %lf\n", (ftemp/ sigma),(grtemp));
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
	double ex[n], ey[n], ez[n];
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
	pdbfile = fopen("pdbtest", "w+");

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

double bump(double sigma, int nextsphere, int temppartner, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[])
{
	double sigsq = pow(sigma, 2.0);

	double rxij = rx[nextsphere] - rx[temppartner];
	double ryij = ry[nextsphere] - ry[temppartner];
	double rzij = rz[nextsphere] - rz[temppartner];

	rxij = rxij - rint(rxij);
	ryij = ryij - rint(ryij);
	rzij = rzij - rint(rzij);


	double factor = (rxij*(vx[nextsphere] - vx[temppartner]) + ryij*(vy[nextsphere] - vy[temppartner]) + rzij*(vz[nextsphere] - vz[temppartner]))/sigsq;

	double delvx = -factor*rxij;
	double delvy = -factor*ryij;
	double delvz = -factor*rzij;


	vx[nextsphere] = vx[nextsphere] + delvx;
	vx[temppartner] = vx[temppartner] - delvx;
	vy[nextsphere] = vy[nextsphere] + delvy;
	vy[temppartner] = vy[temppartner] - delvy;
	vz[nextsphere] = vz[nextsphere] + delvz;
	vz[temppartner] = vz[temppartner] - delvz;

	double w = delvx*rxij + delvy*ryij + delvz*rzij;

	return w;
}

void dnlist(double sigma, int j, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], double coltime[], int partner[])
{
	if(j == 0)
	{
		return;
	}

	double sigsq = pow(sigma,2);
	double rxj, ryj, rzj, vxj, vyj, vzj;
	double rxij, ryij, rzij, vxij, vyij, vzij;
	double bij;

	rxj = rx[j];
	ryj = ry[j];
	rzj = rz[j];
	vxj = vx[j];
	vyj = vy[j];
	vzj = vz[j];

	for(int i = 0; i < j; i++)
	{
		rxij = rx[i] - rxj;
		ryij = ry[i] - ryj;
		rzij = rz[i] - rzj;
		rxij = rxij - rint(rxij);
		ryij = ryij -  rint(ryij);
		rzij = rzij - rint(rzij);
		vxij = vx[i] - vxj;
		vyij = vy[i] - vyj;
		vzij = vz[i] - vzj;
		bij = rxij*vxij + ryij*vyij + rzij*vzij;

		if(bij < 0.0)
		{
			double rijsq = pow(rxij,2) + pow(ryij,2) + pow(rzij,2);
			double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
			double discr = pow(bij,2) - vijsq * (rijsq - sigsq);

			if(discr > 0.0)
			{
				double tij = (-bij - sqrt(discr)) / vijsq;
				if(tij < coltime[i])
				{
					coltime[i] = tij;
					partner[i] = j;
				}
			}
		}
	}

}

void uplist(double sigma, int i, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], double *coltime, int *partner)
{
	if(i == (n-1))
	{
		return;
	}

	double rxi, ryi, rzi, vxi, vyi, vzi;
	double rxij, ryij, rzij, vxij, vyij, vzij;
	double bij;

	double sigsq = pow(sigma, 2);
	*coltime =  timbig;
	rxi = rx[i];
	ryi = ry[i];
	rzi = rz[i];
	vxi = vx[i];
	vyi = vy[i];
	vzi = vz[i];

	for(int j = i+1; j < n; j++)
	{
		rxij = rxi - rx[j];
		ryij = ryi - ry[j];
		rzij = rzi - rz[j];

		rxij = rxij - rint(rxij);
		ryij = ryij - rint(ryij);
		rzij = rzij - rint(rzij);

		vxij = vxi - vx[j];
		vyij = vyi - vy[j];
		vzij = vzi - vz[j];
		bij = rxij*vxij + ryij*vyij + rzij*vzij;
		// was missing a j on vzij last one

		if(bij < 0.0)
		{
			double rijsq = pow(rxij,2) + pow(ryij,2) + pow(rzij,2);
			double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
			double discr = pow(bij,2) - vijsq*(rijsq - sigsq);

			if(discr > 0.0)
			{
				double tij = (-bij - sqrt(discr))/vijsq;
				if(tij < *coltime)
				{
					*coltime = tij;
					*partner = j;
				}
			}
		}

	}
}

/*
 * Does something with the radial distribution function
 */
void grsort(double rx[], double ry[], double rz[], int hist[])
{
	for(int i = 0; i < (n-1); i++)
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


			double rijsq = pow(rxij, 2.0) + pow(ryij, 2.0) + pow(rzij, 2.0);
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

			rxij = rxij - rint(rxij);
			ryij = ryij - rint(ryij);
			rzij = rzij - rint(rzij);

			//printf("%f %f %f\n", rxij, ryij, rzij);
			double rijsq = pow(rxij, 2) + pow(ryij, 2) + pow(rzij, 2);

			if (rijsq < sigsq)
			{
				rij = sqrt(rijsq / sigsq);
				if((1.0 - rij) > tol)
				{
					//printf("r%d %d: %f, tol: %f\n", i,j,rij, tol);
					overlap = 1;
				}
			}
		}
	}
	return overlap;
}

/*
 * Determines energy of the structure from velocity data
 */
double energy(double vx[], double vy[], double vz[])
{
	double e = 0.0;
	for(int i = 0; i < n; i++)
	{
		e = e + pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2);
	}
	e = 0.5*e;
	return e;
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
	//printf("gauss: %f\n", gauss);
	return gauss;
}

/*
 * Change this random number generator, this one is not good
 */
double ranf()
{
	double dra;
	dra = rand() / (double)(RAND_MAX);
	return dra;

}
