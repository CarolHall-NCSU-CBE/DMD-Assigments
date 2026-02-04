//
// Square-well Code to calculate compressibility factor and radial distribution function
//
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <time.h>
using namespace std;

// Define necessary global parameters
int const n = 108;
double const timbig = 1.0E10;
double const pi = 3.14159265359;
double const tol = 1.0E-4;
int const maxbin = 7000;
double const delr = 0.0001;
double const ghostcoeff = 0.001;

// Create initial configuration from FCC lattice
void createcoord(double rx[], double ry[], double rz[])
{
    // Error message for memory allocation in initial configuration
    double (*ex), (*ey), (*ez);
    ex = (double*)malloc(n*sizeof(double));
	ey = (double*)malloc(n*sizeof(double));
	ez = (double*)malloc(n*sizeof(double));

    if ((ex == NULL) || (ey == NULL) || (ez == NULL))
    {
        cout << "Out of memory in initial coordinate creation" << endl;
        exit(0);
    }

    double temprx[4], tempry[4], temprz[4];                 // Declare temporary position variables -- using r conflicts with previously declared variables

    int const nc = ceil(pow(n / 4.0, 1.0/3.0));             // Define number of cell types
    double cell, cell2;                                     // Declare cell variables
    int m, ix, iy, iz, iref;                                // Declare counters
    cell = 1.0 / double(nc);
    cell2 = 0.5 * cell;

    // Build unit cell
    temprx[0] = 0.0;                                            // Sub-lattice A
    tempry[0] = 0.0;
    temprz[0] = 0.0;

    temprx[1] = cell2;                                          // Sub-lattice B
    tempry[1] = cell2;
    temprz[1] = 0.0;

    temprx[2] = 0.0;                                            // Sub-lattice C
    tempry[2] = cell2;
    temprz[2] = cell2;

    temprx[3] = cell2;                                          // Sub-lattice D
    tempry[3] = 0.0;
    temprz[3] = cell2;

    // Construct lattice from unit cell
    m = 0;                                                      // Initialize counter m
    for (iz = 0; iz < nc; ++iz)
    {
        for (iy = 0; iy < nc; ++iy)
        {
            for (ix = 0; ix < nc; ++ix)
            {
                for (iref = 0; iref < 4; ++iref)
                {
                    rx[iref + m] = temprx[iref] + cell * double(ix); // Place copy of unit cell adjacent to original cell on x-axis
                    ry[iref + m] = tempry[iref] + cell * double(iy); // Place copy of unit cell adjacent to original cell on y-axis
                    rz[iref + m] = temprz[iref] + cell * double(iz); // Place copy of unit cell adjacent to original cell on z-axis

                    // Shift center of box to origin
                    rx[iref + m] = rx[iref + m] - 0.5;
                    ry[iref + m] = ry[iref + m] - 0.5;
                    rz[iref + m] = rz[iref + m] - 0.5;
                }
                m = m + 4;                                     // Now move to next spot to fill
            }
        }
    }
}

// drandm - Random number generator
double drandm()
{
    double drandm;                                          // Declare drandm as variable for random distribution
    drandm = rand() / (double)RAND_MAX;                     // Give a random number between 0 and 1
    return drandm;
}

// Gauss - real function returning a uniform random normal variate
// from a distribution with zero mean and unit variance (reference: pg. 348 in Allen and Tildesley)
double gauss()
{
	double a1 = 3.949846138, a3 = 0.252408784, a5 = 0.076542912, a7 = 0.00835596, a9 = 0.029899776;  // Gauss constants
	double r, r2;
	double rsum = 0.0;

	for(int i = 0; i < 12; i++)
	{
		rsum = rsum + drandm();
	}

	r = (rsum - 6.0) / 4.0;
	r2 = r * r;

	double gauss = (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 +a1 ) * r;

	return gauss;
}
// Assign velocities to particles from Maxwell-Boltzmann distribution
void createvel(double temp, double vx[], double vy[], double vz[])
{
    double rtemp = sqrt(temp);
    double sumx = 0.0;
    double sumy = 0.0;
    double sumz = 0.0;

    for (int i = 0; i < n; ++i)
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

    for (int j = 0; j < n; ++j)
    {
        vx[j] = vx[j] - sumx;
        vy[j] = vy[j] - sumy;
        vz[j] = vz[j] - sumz;
    }
}

//	writePeriodicCoordinates - writes out the periodic coordinates to the specified file in the format necessary for vmd pdb files
void writePeriodicCoordinates(FILE *PDBFile, double *rx, double *ry, double *rz)
{
	for(int i = 0; i < n; i++)
	{
		double rxi = rx[i];
		double ryi = ry[i];
		double rzi = rz[i];
		double boxlx, boxly, boxlz;
		boxlx = 1.0;
		boxly = 1.0;
		boxlz = 1.0;

		rxi = rxi - (boxlx*(double)((int)(rxi/(double)(boxlx) + (0.5*copysign(1.0,rxi)))));
		ryi = ryi - (boxly*(double)((int)(ryi/(double)(boxly) + (0.5*copysign(1.0,ryi)))));
		rzi = rzi - (boxlz*(double)((int)(rzi/(double)(boxlz) + (0.5*copysign(1.0,rzi)))));
		fprintf(PDBFile, "ATOM %*d  N%*d", 6, i, 6, 0);
		fprintf(PDBFile, "%*.3f%*.3f%*.3f\n", 17, 100.0*rxi, 8, 100.0*ryi, 8, 100.0*rzi);
		// Note that the coordinates written out are 100x that of actual values
	}
}

// Check - checks for overlaps in the configuration
unsigned char check(double sigma, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[])
{
	double sigmasq = pow(sigma, 2.0);
	double rij, rxij, ryij, rzij;
	unsigned char overlap = 0;                              // Default value of overlap set at 0 for no overlap

	for(int i = 0; i < (n - 1); i++)
	{
		double rxi = rx[i];
		double ryi = ry[i];
		double rzi = rz[i];

		for(int j = i + 1; j < n; j++)
		{
			rxij = rxi - rx[j];
			ryij = ryi - ry[j];
			rzij = rzi - rz[j];

			rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0) + (0.5*copysign(1.0,rxij)))));
			ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0) + (0.5*copysign(1.0,ryij)))));
			rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0) + (0.5*copysign(1.0,rzij)))));

			double rijsq = pow(rxij, 2) + pow(ryij, 2) + pow(rzij, 2);

			if (rijsq < sigmasq)
			{
				rij = sqrt(rijsq / sigmasq);
				if((1.0 - rij) > tol)
				{
                    //cout << "i, j, rij/sigma = " << i << ", " << j << ", " << rij << endl;
					overlap = 1.0;                          // If there is an overlap, return this value
				}
			}
		}
	}
	return overlap;
}

// Kineticen - determine kinetic energy of particles from velocities
double kineticen(double vx[], double vy[], double vz[])
{
    double e = 0.0;

    for (int i = 0; i < n; ++i)
    {
        e = e + pow(vx[i], 2.0) + pow(vy[i], 2.0) + pow(vz[i], 2.0);
    }

    e = 0.5*e;
    return e;
}

// Potentialen - determine potential energy of system from positions
double potentialen(double rx[], double ry[], double rz[], double sigma2, double energy)
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

// Ghostcoll - perform ghost collision
void ghostcoll(double usrtemp, double *vx, double *vy, double *vz)
{
    double rtemp = sqrt(usrtemp);
    *vx = rtemp * gauss();
    *vy = rtemp * gauss();
    *vz = rtemp * gauss();
}

// Bump - compute collision dynamics for particles i and j (assuming i and j are in contact) and collision Virial w
double bump(double sigma, double sigma2, double energy, int nextsphere, int temppartnr, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], int *colltype)
{
    double sigmasq = pow(sigma, 2.0);
    double sigma2sq = pow(sigma2, 2.0);

    double rxij = rx[nextsphere] - rx[temppartnr];
    double ryij = ry[nextsphere] - ry[temppartnr];
    double rzij = rz[nextsphere] - rz[temppartnr];

    rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0) + (0.5*copysign(1.0,rxij)))));
	ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0) + (0.5*copysign(1.0,ryij)))));
	rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0) + (0.5*copysign(1.0,rzij)))));

	double vxij = vx[nextsphere] - vx[temppartnr];
	double vyij = vy[nextsphere] - vy[temppartnr];
	double vzij = vz[nextsphere] - vz[temppartnr];

	double bij = rxij*vxij + ryij*vyij + rzij*vzij;
	double bijsq = pow(bij, 2.0);

	double delvx = 0;
	double delvy = 0;
	double delvz = 0;

	double smdist = 5.0/pow(10.0,11.0);
	double smbump = 0;                                  // Adjust based on what kind of collision occurs

	// Hard Sphere Collision (Core Collision)
	if(*colltype == 1)
	{
		double factor = (bij)/sigmasq;
        delvx = -factor*rxij;
        delvy = -factor*ryij;
        delvz = -factor*rzij;
		smbump = 1;
	}

	// Attractive Collision (Bounce or Disassociation)
	else if(*colltype == 2)
	{
		double compare = 4*sigma2sq*energy;

		// If particle has enough energy to get out of well
		if(bijsq > compare)
		{
			double parenterm = -sqrt(-4*sigma2sq*energy + bijsq) + bij;
			double factor = 1/(2.0*sigma2sq) * parenterm;
            delvx = -factor*rxij;
            delvy = -factor*ryij;
            delvz = -factor*rzij;

			//cout << "Dissociation" << endl;
			smbump = 1;
		}

		// If particle does not have enough energy to get out of well
		else if(bijsq <= compare)
		{
			double factor = (bij)/sigma2sq;
            delvx = -factor*rxij;
            delvy = -factor*ryij;
            delvz = -factor*rzij;

			//cout << "Square Well Bounce" << endl;
			smbump = -1;
		}
	}
	// Capture Collision
	else if(*colltype == 3)
	{
		double root = 4.0*sigma2sq*energy + bijsq;
		double parenterm = sqrt(root) + bij;
		double factor = 1.0/(2.0*sigma2sq) * parenterm;

        delvx = -factor*rxij;
        delvy = -factor*ryij;
        delvz = -factor*rzij;

		smbump = -1;
	}

	vx[nextsphere] = vx[nextsphere] + delvx;
	vx[temppartnr] = vx[temppartnr] - delvx;
	vy[nextsphere] = vy[nextsphere] + delvy;
	vy[temppartnr] = vy[temppartnr] - delvy;
	vz[nextsphere] = vz[nextsphere] + delvz;
	vz[temppartnr] = vz[temppartnr] - delvz;

	vxij = vx[nextsphere] - vx[temppartnr];
	vyij = vy[nextsphere] - vy[temppartnr];
	vzij = vz[nextsphere] - vz[temppartnr];

	double rij = sqrt( pow(rxij,2.0) + pow(ryij,2.0) + pow(rzij,2.0));

	rx[nextsphere] = rx[nextsphere] + smbump*smdist*rxij/rij;
	ry[nextsphere] = ry[nextsphere] + smbump*smdist*ryij/rij;
	rz[nextsphere] = rz[nextsphere] + smbump*smdist*rzij/rij;

	rx[temppartnr] = rx[temppartnr] - smbump*smdist*rxij/rij;
	ry[temppartnr] = ry[temppartnr] - smbump*smdist*ryij/rij;
	rz[temppartnr] = rz[temppartnr] - smbump*smdist*rzij/rij;

	double w = delvx*rxij + delvy*ryij + delvz*rzij;

	return w;
}

// Dnlist - Locate other possible collision in list of particle partners
void dnlist(double sigma, double sigma2, int j, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], double coltime[], int partner[], int colltype[])
{
    double sigsq = pow(sigma,2.0);
	double sig2sq = pow(sigma2, 2.0);
	double rxj, ryj, rzj, vxj, vyj, vzj;
	double rxij, ryij, rzij, vxij, vyij, vzij;
	double bij;

	// When we get to the last item in the list, stop routine
	if(j == 0)
	{
		return;
	}

	rxj = rx[j];
	ryj = ry[j];
	rzj = rz[j];
	vxj = vx[j];
	vyj = vy[j];
	vzj = vz[j];

	double tij = timbig;
	int tempcolltype;

	for(int i = 0; i < j; ++i)                                          // Calculate bijs for all particles in system
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

		// Case 1 - particles approaching each other
		if(bij < 0.0)
		{
			double rijsq = pow(rxij,2) + pow(ryij,2) + pow(rzij,2);

			// CASE 1A - particles approaching each other inside sigma2 reach
			if((rijsq - sig2sq) <= 0)
			{
				double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
				double discr = pow(bij,2) - vijsq*(rijsq - sigsq);

				// Case 1Ai - particles undergo hard sphere / core collision
				if(discr > 0)
				{
					tij = (-bij - sqrt(discr))/vijsq;
					tempcolltype = 1;                                   // Assign temporary collision type to hard sphere collision
				}

				// Case 1Aii - particles bounce or dissociate
				else if(discr <= 0)
				{
					double discr2 = pow(bij,2) - vijsq*(rijsq - sig2sq);
					tij = (-bij + sqrt(discr2))/vijsq;
					tempcolltype = 2;                                   // Assign temporary collision type to bounce/dissociate collision
				}
			}

			// Case 1B - particles approaching each outside of sigma reach
			else if((rijsq-sig2sq) > 0)
			{
				double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
				double discr2 = pow(bij,2) - vijsq*(rijsq - sig2sq);

				// Case 1Bi - particles undergo capture collision
				if(discr2 >= 0)
				{
					tij = (-bij - sqrt(discr2))/vijsq;
					tempcolltype = 3;                                   // Assign temporary collision type to capture collision
				}
			}
		}

		// Case 2 - particles receding from each other, but still inside sigma2 reach
		else if(bij > 0.0)
		{
			double rijsq = pow(rxij,2) + pow(ryij,2) + pow(rzij,2);

			// Case 2A - particles undergo bounce/dissociate collision
			if((rijsq - sig2sq) < 0)
			{
				double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
				double discr2 = pow(bij,2) - vijsq*(rijsq - sig2sq);

				tij = (-bij + sqrt(discr2))/vijsq;
				tempcolltype = 2;
			}

		}

		// Now change list details after collision
		if(tij < coltime[i])
		{
			coltime[i] = tij;
			partner[i] = j;
			colltype[i] = tempcolltype;
		}
	}
}

// Uplist - Update the other partners and collisions after a collision is made
void uplist(double sigma, double sigma2, int i, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], double *coltime, int *partner, int *colltype)
{
    double rxi, ryi, rzi, vxi, vyi, vzi;
	double rxij, ryij, rzij, vxij, vyij, vzij;
	double bij;

	// When we get to the last item in the list, stop routine
	if(i == (n-1))
	{
		return;
	}

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

	for(int j = i+1; j < n; ++j)
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

		tij = timbig;
		int tempcolltype;

		// Case 1 - particles approached each other
		if(bij < 0.0)
		{
			double rijsq = pow(rxij,2) + pow(ryij,2) + pow(rzij,2);

			// Case 1A - particles approached each other within sigma2 reach
			if((rijsq - sig2sq) <= 0)                                       // Less than or equal to, because after bump, particles are considered to be inside the well
			{
				double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
				double discr = pow(bij,2) - vijsq*(rijsq - sigsq);

				// Case 1Ai - particles underwent hard sphere / core collision
				if(discr > 0)
				{
					tij = (-bij - sqrt(discr))/vijsq;
					tempcolltype = 1;                                       // Assign temporary collision type to hard sphere / core collision
				}

				// Case 1Aii - particles underwent dissociation/bounce collision
				else
				{
					double discr2 = pow(bij,2) - vijsq*(rijsq - sig2sq);
					/*if(discr2 < 0)
					{
						cout << "discr2 < 0 in uplist" << endl;
						exit(0);
					}*/
					tij = (-bij + sqrt(discr2))/vijsq;
					tempcolltype = 2;                                       // Assign temporary collision type to dissociation/bounce collision
				}
			}

			// Case 1B - particles approaching each other outside of sigma reach
			else
			{
				double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
				double discr2 = pow(bij,2) - vijsq*(rijsq - sig2sq);

				// Case 1Bi - particles underwent capture collision
				if(discr2 >= 0)
				{
					tij = (-bij - sqrt(discr2))/vijsq;
					tempcolltype = 3;                                       // Assign temporary collision type to capture collision
				}
			}
		}

		// Case 2 - particles receded from each other, but within sigma2 reach
		else if(bij > 0.0)
		{
			double rijsq = pow(rxij,2) + pow(ryij,2) + pow(rzij,2);

			// Case 2A - particles underwent dissociation/bounce collision
			if((rijsq - sig2sq) < 0)                                        // Less than because particles are out of the well after collision
			{
				double vijsq = pow(vxij,2) + pow(vyij,2) + pow(vzij,2);
				double discr2 = pow(bij,2) - vijsq*(rijsq - sig2sq);
				tij = (-bij + sqrt(discr2))/vijsq;
				tempcolltype = 2;                                           // Assign temporary collision type to dissociation/bounce collision
			}

		}

		// Now change list details after collision
		if(tij < *coltime)
		{
			*colltype = tempcolltype;
			*coltime = tij;
			*partner = j;
		}
	}
}

// Grsort - Set up radial distribution calculations - sorts radial distances into bins based on a scale of rlower to rupper with maxbin divisions
void grsort(double rx[], double ry[], double rz[], int hist[])
{
 	for(int i = 0; i < (n-1); ++i)
	{
		double rxi = rx[i];
		double ryi = ry[i];
		double rzi = rz[i];

		double rxij;
		double ryij;
		double rzij;
		for(int j = (i+1); j < n; ++j)
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

int main()
{
    // Activate random number generator
    srand((unsigned int)time(NULL));

    // Declare input variables
    double density;
    double temp;
    int ncoll;

    // Dynamically allocate memory for position and velocity arrays because this is about to get crazy~
    double *(rx), *(ry), *(rz);
    double *(vx), *(vy), *(vz);
    double *(temprx), *(tempry), *(temprz);
    int *(colltype);

    // Error messages for memory
    rx = (double*)malloc(n*sizeof(double));
    if (rx == NULL)
    {
        cout << "Out of memory in rx array" << endl;
        exit(0);
    }

    ry = (double*)malloc(n*sizeof(double));
    if (ry == NULL)
    {
        cout << "Out of memory in ry array" << endl;
        exit(0);
    }

    rz = (double*)malloc(n*sizeof(double));
    if (rz == NULL)
    {
        cout << "Out of memory in rz array" << endl;
        exit(0);
    }

    vx = (double*)malloc(n*sizeof(double));
    if (vx == NULL)
    {
        cout << "Out of memory in vx array" << endl;
        exit(0);
    }

    vy = (double*)malloc(n*sizeof(double));
    if (vy == NULL)
    {
        cout << "Out of memory in vy array" << endl;
        exit(0);
    }

    vz = (double*)malloc(n*sizeof(double));
    if (vz == NULL)
    {
        cout << "Out of memory in vz array" << endl;
        exit(0);
    }

    temprx = (double*)malloc(n*sizeof(double));
    if (temprx == NULL)
    {
        cout << "Out of memory in temprx array" << endl;
        exit(0);
    }

    tempry = (double*)malloc(n*sizeof(double));
    if (tempry == NULL)
    {
        cout << "Out of memory in tempry array" << endl;
        exit(0);
    }

    temprz = (double*)malloc(n*sizeof(double));
    if (temprz == NULL)
    {
        cout << "Out of memory in temprz array" << endl;
        exit(0);
    }

    colltype = (int*)malloc(n*sizeof(int));
    if (colltype == NULL)
    {
        cout << "Out of memory in colltype array" << endl;
        exit(0);
    }

    // Declare dynamic arrays for collision times and partners
    double e = 0.0;
    unsigned char overlap;
    double *coltime;
    int *partner;

    // Error messages for memory
    coltime = (double*)malloc((n+2)*sizeof(double));
    if (coltime == NULL)
    {
        cout << "Out of memory in coltime array" << endl;
        exit(0);
    }

    partner = (int*)malloc(n*sizeof(int));
    if (partner == NULL)
    {
        cout << "Out of memory in partner list" << endl;
        exit(0);
    }

    // Set up hist and colltype lists
    double w, acw;
    double tij;
    int hist[maxbin];

    for (int i = 0; i < maxbin; ++i)
    {
        hist[i] = 0;
    }

    for (int j = 0; j < n; ++j)
    {
        colltype[j] = 0;
    }

    double initpotentialen = 1.0;

    // Prompt and receive inputs
    cout << "Program Square Well" << endl;
    cout << "Molecular Dynamics of Square Well Spheres" << endl;
    cout << "Results in Units kt = sigma = 1" << endl;

    cout << "Enter reduced density [(n/v)*(sigma^3)]: ";
    cin >> density;

    cout << "Enter initial reduced temperature: ";
    cin >> temp;

    cout << "Enter number of collisions required: ";
    cin >> ncoll;

    // Set up initial parameters
    double sigma = pow((density/(double)n), 1.0/3.0);
    double sigma2 = (1 + 0.5) * sigma;

    if (sigma2 >= 0.5)
    {
        cout << "Square well diameter too large" << endl;
        exit(0);
    }

    // Build initial configuration of system
    createcoord(rx, ry, rz);
    createvel(10.0, vx, vy, vz);

    // Create file to look at initial configuration
    FILE *writeoutinit;
    char buf[1000];
    sprintf(buf, "initialconfig");
    writeoutinit = fopen(buf, "w+");
    fprintf(writeoutinit, "MODEL%*d\n", 9, 0);
    fprintf(writeoutinit, "CRYST1%*.3f%*.3f%*.3f%*.2f%*.2f%*.2f%*s%*d\n", 9, 100.0, 9, 100.0, 9, 100.0, 7, 90.00, 7, 90.00, 7, 90.00, 11, "P", 4, 1);
    writePeriodicCoordinates(writeoutinit, rx, ry, rz);
    fprintf(writeoutinit, "ENDMDL\n");

    // Check for overlaps
    overlap = check(sigma, rx, ry, rz, vx, vy, vz);
    if (overlap == 1)
    {
        cout << "Particle overlap in initial configuration" << endl;
        return 0;
    }

    // Calculate initial energy of system
    e = kineticen(vx, vy, vz);
    double totalen = e + potentialen(rx, ry, rz, sigma2, initpotentialen);
    cout << "Initial energy of system: " << totalen << endl;

    double en = e / ((double)n);
    double calctemp = 2.0 * en / 3.0;
    cout << "Initial temperature: " << calctemp << endl;

    double enkt = en/calctemp;
    cout << "Initial e/nkt: " << enkt << endl;

    for (int i = 0; i < n; ++i)
    {
        coltime[i] = timbig;
        partner[i] = n;
    }

    for (int j = 0; j < n; ++j)
    {
        uplist(sigma, sigma2, j, rx, ry, rz, vx, vy, vz, &(coltime[j]), &(partner[j]), &(colltype[j]));
    }

    // Set up ghost particle to prompt calculations
    coltime[n] = 5.0;

    double ghosttime = ghostcoeff * fabs(sqrt(-2.0*log(((rand()%100)+0.5)/100.0)) * cos(2.0*pi*(rand()%100)/100.0));
    cout << "Initial ghost time: " << ghosttime << endl;
    coltime[n+1] = ghosttime;

    // Zero Virial Accumulator
    acw = 0.0;

    // Beginning of Main Loop
    cout << "**Starting Dynamics...**" << endl;

    int steps = 0;
    double t = 0.0;
    int uplistcntr = 10;
    int burnin = 100000;
    int grcount = 2000000;
    int nextsphere;
    int temppartnr;
    const double deltawrite = 0.005;
    int tcountwrite = 0;

    // Print out file configurations
    FILE *writeout;
    writeout = fopen("runconfigs", "w+");
    fprintf(writeout, "REMARK  deltat: %f\n", deltawrite);

    fprintf(writeout, "MODEL%*d\n", 9, tcountwrite);
    for (int i = 0; i < n; i++)
    {
        double rxi = rx[i];
        double ryi = ry[i];
        double rzi = rz[i];

        fprintf(writeout, "ATOM %*d  N", 6, i);
		fprintf(writeout, "%*.3f%*.3f%*.3f\n", 24, rxi, 8, ryi, 8, rzi);
    }
    fprintf(writeout, "ENDMDL\n");
    tcountwrite += 1;

    double burnintime = 0.0;

    // Print out temperature recordings
    FILE *temprec;
    temprec = fopen("temprec", "w+");
    unsigned int ghoston = 1;
    for (int coll = 0; coll < ncoll; ++coll)
    {
        // Record temperature before each collision
        e = kineticen(vx, vy, vz);
        en = e / ((double)n);
        calctemp = 2.0 * en / 3.0;
        fprintf(temprec, "%d %f\n", coll, calctemp);

        // Switches to user defined temperature (lower) after higher, burnin temperature
        if (coll == burnin)
        {
            cout << "Burn in period complete at " << coll << " collisions" << endl;
            createvel(temp, vx, vy, vz);
            for (int i = 0; i < n; ++i)
            {
                coltime[i] = timbig;
                partner[i] = n;
            }

            for (int j = 0; j < n; ++j)
            {
                uplist(sigma, sigma2, j, rx, ry, rz, vx, vy, vz, &(coltime[j]), &(partner[j]), &(colltype[j]));
            }

            // Reset variable for main part of loop after burn in
            coltime[n] = 5.0;
            e = kineticen(vx, vy, vz);
            totalen = e + potentialen(rx, ry, rz, sigma2, initpotentialen);
            acw = 0;
            burnintime = t;
            ghosttime = ghostcoeff * fabs(sqrt(-2.0*log(((rand()%100)+0.5)/100.0))*cos(2.0*pi*(rand()%100)/100.0));
            coltime[n+1] = ghosttime;
            cout << "New total energy: " << totalen << endl;
        }

        tij = timbig;
        for (int k = 0; k < (n+1); ++k)
        {
            if (coltime[k] < tij && k != (n-1))
            {
                tij = coltime[k];
                nextsphere = k;
            }
        }

        // Check if ghost collision time smaller than normal particle collision time
        if (coltime[n+1] < tij && coll > burnin)
        {
            tij = coltime[n+1];
            nextsphere = n+1;
        }

        // Check if there were collision times less than allowed (negative collisions)
        if (tij <= 0)
        {
            cout << "Negative collision at t = " << tij << ", with sphere number" << nextsphere << "and collision " << coll << endl;
            exit(0);
        }

        // Write out configurations to file
        while ((t+tij) >= (tcountwrite*deltawrite))
        {
            double temptime = tcountwrite * deltawrite - t;
            fprintf(writeout, "MODEL%*d\n", 9, tcountwrite);
            for (int i = 0; i < n; ++i)
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

				temprx[i] = rxi;
				tempry[i] = ryi;
				temprz[i] = rzi;

				fprintf(writeout, "ATOM %*d  N", 6, i);
				fprintf(writeout, "%*.3f%*.3f%*.3f\n", 24, rxi,8,ryi,8,rzi);
            }
            if (coll > grcount)
            {
                grsort(temprx, tempry, temprz, hist);
                grcount = grcount + 20;
                coltime[n] = t + 5.0;
                steps = steps + 1;
            }
            //int tempcoll = coll;

            fprintf(writeout, "ENDMDL\n");
            tcountwrite = tcountwrite + 1;
        }

        // Build radial distribution function
        // if (grcount == coll || nextsphere == n)
        //{
        //    grsort(rx, ry, rz, hist);
        //    grcount = grcount + 20;
        //    coltime[n] = t + 5.0;
        //    steps = steps + 1;
        //}

        // Move particles forward by time tij and reduce collision times
        t = t + tij;
        for (int k = 0; k < n; ++k)
        {
            coltime[k] = coltime[k] - tij;
            rx[k] = rx[k] + vx[k]*tij;
			ry[k] = ry[k] + vy[k]*tij;
			rz[k] = rz[k] + vz[k]*tij;
			rx[k] = rx[k] - (1.0*(double)((int)(rx[k]/(double)(1.0) + (0.5*copysign(1.0,rx[k])))));
			ry[k] = ry[k] - (1.0*(double)((int)(ry[k]/(double)(1.0) + (0.5*copysign(1.0,ry[k])))));
			rz[k] = rz[k] - (1.0*(double)((int)(rz[k]/(double)(1.0) + (0.5*copysign(1.0,rz[k])))));
        }

        if (ghoston == 1)
        {
            coltime[n+1] = coltime[n+1] - tij;
        }

        if (nextsphere < n)
        {
            temppartnr = partner[nextsphere];
            w = bump(sigma, sigma2, initpotentialen, nextsphere, temppartnr, rx, ry, rz, vx, vy, vz, &(colltype[nextsphere]));
            acw = acw + w;

            for (int l = 0; l < n; ++l)
            {
                if ((l == nextsphere) || (partner[l] == nextsphere) || (l == temppartnr) || (partner[l] == temppartnr))
                {
                    uplist(sigma, sigma2, l, rx, ry, rz, vx, vy, vz, &(coltime[l]), &(partner[l]), &(colltype[l]));
                }
            }

            dnlist(sigma, sigma2, nextsphere, rx, ry, rz, vx, vy, vz, coltime, partner, colltype);
            dnlist(sigma, sigma2, temppartnr, rx, ry, rz, vx, vy, vz, coltime, partner, colltype);

            // Update uplist if necessary
            if ((coll + 1) == (uplistcntr - 1))
            {
                for (int i = 0; i < n; ++i)
                {
                    uplist(sigma, sigma2, i, rx, ry, rz, vx, vy, vz, &(coltime[i]), &(partner[i]), &(colltype[i]));
                }
                uplistcntr = uplistcntr + 10;
            }
        }
        else if (coll > burnin && nextsphere == n+1)
        {
            int rndsphere = (int)(drandm()*n);
            if (rndsphere == n)
            {
                cout << "Wrong particle chosen" << endl;
            }

            ghostcoll(temp, &(vx[rndsphere]), &(vy[rndsphere]), &(vz[rndsphere]));

            ghosttime = ghostcoeff*fabs(sqrt(-2.0*log(((rand()%100)+0.5)/100.0))*cos(2.0*pi*(rand()%100)/100.0));

            // Turn off ghost collisions
            // if (ghoston == 0)
            //{
            //    ghosttime = 10000;
            //}

            while (ghosttime <= 0)
            {
                ghosttime = ghostcoeff * fabs(sqrt(-2.0*log(((rand()%100)+1)/100.0))*cos(2.0*pi*(rand()%100)/100.0));
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
    }
    fclose(writeout);
    fclose(temprec);

    // End main loop

    cout << "**End of Dynamics**" << endl;
    cout << "Final colliding pair: " << nextsphere << ", " << temppartnr << endl;

    overlap = check(sigma, rx, ry, rz, vx, vy, vz);
    if (overlap == 1)
    {
        cout << "Particle overlap in final configuration" << endl;
        return 0;
    }

    // Write out ending configuration
    FILE *endfile;
    endfile = fopen("endconfig", "w+");
    for (int i =0; i < n; ++i)
    {
        fprintf(endfile, "ATOM %*d  N", 6, i);
		fprintf(endfile, "%*.3f%*.3f%*.3f\n", 24, rx[i],8,ry[i],8,rz[i]);
	}
	fclose(endfile);

    e = kineticen(vx, vy, vz);
    double pvnkt1 = (acw / ((double)n*3.0*(t-burnintime)*calctemp)) + 1;
    en = e/((double)n);
    calctemp = 2.0 * en / 3.0;

    totalen = e + potentialen(rx, ry, rz, sigma2, initpotentialen);

    enkt = en / calctemp;
    t = t * sqrt(calctemp) / sigma;
    double rate = ((double)ncoll) / (t - burnintime);
    double tbc = ((double)n)/rate/2.0;
    double frac = pi * density / 6;

    cout << "The final n is: " << n << endl;
    cout << "The final e is: " << totalen << endl;
    cout << "The final temperature is: " << calctemp << endl;
    cout << "The final volume fraction is: " << frac << endl;
    cout << "The final Virial accumulator count is: " << acw << endl;
    cout << "The final time is: " << t << endl;
    cout << "The collision rate is: " << rate << endl;
    cout << "The mean collision time is: " << tbc << endl;
    cout << "The final e/nkt is: " << enkt << endl;
    cout << "Z = PV/nkt = " << pvnkt1 << endl;
    cout << "The final grcount is: " << grcount << endl;
    //fflush(stdout);

    // Calculate g(r) from grsort data
    double grconst = 4.0*pi*(density/(pow(sigma, 3.0)))/3.0;
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

    return 0;
    system("PAUSE");
}

