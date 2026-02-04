//
// Hard Sphere Code to compute compressibility factor and radial distribution function
//
#include <iostream>
#include <iomanip>
//#include <cmath>
//#include <cstdlib>
//#include <ctime>
#include <math.h>
//#include <cctype>
//#include <array>
#include <fstream>
#include <stdlib.h>

using namespace std;

// Define necessary parameters
int const n = 500;                                          // Define sample size and number of bins
double const timbig = 1.0e10;                               // Define time
double const pi = 3.14159265359;                            // Define mathematical constants
double rx[n], ry[n], rz[n], vx[n], vy[n], vz[n];            // Declare position and velocity variables
double coltim[n+1];                                         // Declare collision time array
double sigma;                                               // Declare energy parameter
int partnr[n];                                              // Declare partner list
unsigned int i, j, k, ncoll, coll;                          // Declare counter and collision variables
double density, dij, tij, t, rate;                          // Declare density and time variables
double e, en, enkt, w, pvnkt1, acw, temp, tbc;              // Declare thermo variables
int grconst, grcount, steps;                                // Declare radial distribution variables
int bin;                                                    // Declare bin variable
int const maxbin = 7000;                                    // Define max number of bins
int uplistcntr, kecntr;                                     // Declare uplist center and kinetic energy center
int hist[maxbin];                                           // Declare hist array
double rlower, rupper, nideal, rsigma;                      // Declare limit variables
double gr[maxbin], f[maxbin];                               // Declare radial distribution variables
double grtime;                                              // Declare radial distribution time
double const delr = 0.0001;                                 // Define delr = half box length divided by maxbin
string title;                                               // Declare title string
int ovrlap;                                                 // Declare variable for overlap result

// Createcoord - create initial position configuration (refer to pg. 169 of Allen and Tildesley)
void createcoord(double *rx, double *ry, double *rz)
{
    int const nc = ceil(pow(n / 4.0, 1.0/3.0));             // Define number of cell types
    double cell, cell2;                                     // Declare cell variables
    int m, ix, iy, iz, iref;                                // Declare counters
    cell = 1.0 / double(nc);
    cell2 = 0.5 * cell;

    // Build unit cell
    rx[0] = 0.0;                                            // Sub-lattice A
    ry[0] = 0.0;
    rz[0] = 0.0;

    rx[1] = cell2;                                          // Sub-lattice B
    ry[1] = cell2;
    rz[1] = 0.0;

    rx[2] = 0.0;                                            // Sub-lattice C
    ry[2] = cell2;
    rz[2] = cell2;

    rx[3] = cell2;                                          // Sub-lattice D
    ry[3] = 0.0;
    rz[3] = cell2;

    // Construct lattice from unit cell
    m = 0;
    for (iz = 0; iz < nc; ++iz)
    {
        for (iy = 0; iy < nc; ++iy)
        {
            for (ix = 0; ix < nc; ++ix)
            {
                for (iref = 0; iref < 4; ++iref)
                {
                    rx[iref + m] = rx[iref] + cell * double(ix); // Place copy of unit cell adjacent to original cell on x-axis
                    ry[iref + m] = ry[iref] + cell * double(iy); // Place copy of unit cell adjacent to original cell on y-axis
                    rz[iref + m] = rz[iref] + cell * double(iz); // Place copy of unit cell adjacent to original cell on z-axis

                    // Shift center of box to origin
                    rx[iref + m] = rx[iref + m] - 0.5;
                    ry[iref + m] = ry[iref + m] - 0.5;
                    rz[iref + m] = rz[iref + m] - 0.5;
                }
                m = m + 4;                                      // Now move to next spot to fill
            }
        }
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
    double const a1 = 3.949846138;                              // Assign coefficients for random
    double const a3 = 0.252408784;                              // nonuniform Gaussian distribution
    double const a5 = 0.076542912;
    double const a7 = 0.008355968;
    double const a9 = 0.029899776;
    double rsum, r, r2;                                         // Declare Knuth variables
    double gauss;                                               // Declare gauss variable
    int i;                                                      // Declare integer counter i

    rsum = 0.0;                                                 // Initialize summation of random variates at 0
    for (i = 0; i < 12; ++i)
    {
        rsum = rsum + drandm();                                 // Go to next number in sequence
    }
    r = (rsum - 6.0) / 4.0;                                     // Calculate Knuth R from random variates
    r2 = r * r;
    gauss = ((((a9 * r2 + a7) * r2 + a5) * r2 + a3) * r2 + a1) * r; // Calculate random number polynomial
    //cout << "gauss = " << gauss << endl;
    return gauss;
}

// Createvel - assign spheres initial velocities based on Gaussian distribution
void createvel(double temp, double *vx, double *vy, double *vz)
{
    double rtemp, sumx, sumy, sumz;                                 // Declare temperature and position variables
    int i;                                                          // Declare counter integer i

    rtemp = sqrt(temp);
    //cout << "rtemp = " << rtemp << endl;
    for (i = 0; i < n; ++i)
    {
        vx[i] = rtemp * gauss();                                      // Assign random initial velocity
        vy[i] = rtemp * gauss();                                      // to all of the particles
        vz[i] = rtemp * gauss();                                      // in the box
    }

    // Remove net momentum from system
    sumx = 0.0;
    sumy = 0.0;
    sumz = 0.0;
    for (i = 0; i < n; ++i)
    {
        sumx = sumx + vx[i];
        sumy = sumy + vy[i];
        sumz = sumz + vz[i];
    }

    sumx = sumx / (double)n;
    sumy = sumy / (double)n;
    sumz = sumz / (double)n;

    for (i = 0; i < n; ++i)
    {
        vx[i] = vx[i] - sumx;
        vy[i] = vy[i] - sumy;
        vz[i] = vz[i] - sumz;
    }
}

double check(double sigma, double *rx, double *ry, double *rz, double *vx, double *vy, double *vz)
{
    double const tol = 1.0e-4;                                              // Define tolerance
    int i, j;                                                               // Declare integers
    double rxi, ryi, rzi, rxij, ryij, rzij, rij;                            // Declare variables
    double rijsq, sigsq;                                                    // Declare more variables

    sigsq = pow(sigma, 2.0);                                                // Define sigma squared
    ovrlap = 0;                                                             // Default: false (0)

    for (i = 0; i < (n-1); ++i)
    {
        rxi = rx[i];                                                        // Find position of i-th
        ryi = ry[i];                                                        // sphere for comparison
        rzi = rz[i];
        for (j = (i+1); j < n; ++j)
            {
                rxij = rxi - rx[j];                                             // Find position from j-th
                ryij = ryi - ry[j];                                             // sphere to original
                rzij = rzi - rz[j];                                             // i-th molecule
                rxij = rxij - floor(rxij + 0.5);                                // By convention,
                ryij = ryij - floor(ryij + 0.5);                                // use the nearest image
                rzij = rzij - floor(rzij + 0.5);

                rijsq = (pow(rxij, 2.0)) + (pow(ryij, 2.0)) + (pow(rzij, 2.0)); // Calculate distance

                if (rijsq < sigsq)
                    {                                                               // Really mostly interested in r/sigma
                        rij = sqrt(rijsq / sigsq);                                  // Calculate and check tolerance
                        if ((1.0 - rij) > tol)
                        {
                            //cout << "i, j, rij/sigma = " << i << ", " << j << ", " << rij << endl;
                            ovrlap = 1.0;                                          // If tolerance is exceeded, there is overlap.
                        }
                    }
            }
        }
    return ovrlap;
}

double energy(double *vx, double *vy, double *vz)
{
    double e = 0.0;

        for (i = 0; i < n; ++i)
        {
            e = e + pow(vx[i], 2.0) + pow(vy[i], 2.0) + pow(vz[i], 2.0);
        }
    e = 0.5*e;
    return e;
}

// Uplist - Locate the next possible collision in list of particle partners
void uplist(double sigma, int i, double *rx, double *ry, double *rz, double *vx, double *vy, double *vz, double *coltim, int *partnr)
{
    //int const n = 500;                                        // Define number of particles
    double const timbig = 1.0e10;                               // Define collision time limit
    int j;                                                      // Declare counter variables j
    double rxi, ryi, rzi, rxij, ryij, rzij;                     // Declare position variables
    double vxi, vyi, vzi, vxij, vyij, vzij;                     // Declare velocity variables
    double rijsq, vijsq, bij, tij, discr, sigsq;                // Declare time and other necessary variables

    // Look for collisions with atoms j > i
    if (i == (n-1))
    {
        return;                                                 // Return to the main function once all particles checked
    }

    sigsq = pow(sigma, 2.0);                                    // Square sigma
    *coltim = timbig;                                           // Set collision time limit
    rxi = rx[i];
    ryi = ry[i];
    rzi = rz[i];
    vxi = vx[i];
    vyi = vy[i];
    vzi = vz[i];

    for (j = (i+1); j < n; ++j)
    {
        rxij = rxi - rx[j];
        ryij = ryi - ry[j];
        rzij = rzi - rz[j];
        rxij = rxij - floor(rxij + 0.5);
        ryij = ryij - floor(ryij + 0.5);
        rzij = rzij - floor(rzij + 0.5);
        vxij = vxi - vx[j];
        vyij = vyi - vy[j];
        vzij = vzi - vz[j];
        bij = (rxij*vxij) + (ryij*vyij) + (rzij*vzij);          // Refer to pg. 102 in Allen and Tildesley
    }

    if (bij < 0.0)                                              // Collision may occur between i and j
    {
        rijsq = pow(rxij, 2.0) + pow(ryij, 2.0) + pow(rzij, 2.0);
        vijsq = pow(vxij, 2.0) + pow(vyij, 2.0) + pow(vzij, 2.0);
        discr = pow(bij, 2.0) - vijsq*(rijsq - sigsq);

        if (discr > 0.0)                                        // Must be positive b/c can't take sqrt of negative
        {
            tij = (-bij - sqrt(discr))/vijsq;                   // Calculate time impact
            if (tij < *coltim)                                  // If time impact smaller than collision time,
            {                                                   // update collision time and partner list
                *coltim = tij;
                *partnr = j;
            }
        }
    }
}

// Dnlist - update the other partners and collisions after a collision is made
void dnlist(double sigma, int j, double *rx, double *ry, double *rz, double *vx, double *vy, double *vz, double *coltim, int *partnr)
{
    int i;                                                      // Declare counter integer i
    double rxj, ryj, rzj, rxij, ryij, rzij;                     // Declare position variables
    double vxj, vyj, vzj, vxij, vyij, vzij;                     // Declare velocity variables
    double rijsq, vijsq, bij, tij, discr, sigsq;                // Declare time and other necessary variables

    // Look for collisions with atoms i < j
    if (j == 0)
    {
        return;                                                 // Return to main function once last particle is checked
    }
    sigsq = pow(sigma, 2.0);                                    // Square sigma
    rxj = rx[j];
    ryj = ry[j];
    rzj = rz[j];
    vxj = vx[j];
    vyj = vy[j];
    vzj = vz[j];

    for (i = 0; i < j; ++i)
    {
        rxij = rx[i] - rxj;
        ryij = ry[i] - ryj;
        rzij = rz[i] - rzj;
        rxij = rxij - floor(rxij + 0.5);
        ryij = ryij - floor(ryij + 0.5);
        rzij = rzij - floor(rzij + 0.5);
        vxij = vx[i] - vxj;
        vyij = vy[i] - vyj;
        vzij = vz[i] - vzj;
        bij = (rxij*vxij) + (ryij*vyij) + (rzij*vzij);          // Refer to pg. 105 in Allen and Tildesley
    }

    if (bij < 0.0)                                              // Collision imminent
    {
        rijsq = pow(rxij, 2.0) + pow(ryij, 2.0) + pow(rzij, 2.0);
        vijsq = pow(vxij, 2.0) + pow(vyij, 2.0) + pow(vzij, 2.0);
        discr = pow(bij, 2.0) - vijsq*(rijsq-sigsq);

        if (discr > 0.0)                                        // Must be positive b/c can't take sqrt of negative
        {
            tij = (-bij - sqrt(discr)) / vijsq;                 // Calculate time impact
            if (tij < coltim[i])                                // If time impact is less than the scheduled collision time,
            {                                                   // update the collision time and partner lists
                coltim[i] = tij;
                partnr[i] = j;
            }
        }
    }
}

// Bump - compute collision dynamics for particles i and j (assuming i and j are in contact) and collision Virial w
double bump(double sigma, int nextsphere, int temppartnr, double *rx, double *ry, double *rz, double *vx, double *vy, double *vz)
{
    double rxij, ryij, rzij, factor;                            // Declare variables
    double delvx, delvy, delvz, sigsq;                          // Declare more variables
    double w;

    sigsq = pow(sigma, 2.0);
    rxij = rx[nextsphere] - rx[temppartnr];
    ryij = ry[nextsphere] - ry[temppartnr];
    rzij = rz[nextsphere] - rz[temppartnr];
    rxij = rxij - floor(rxij + 0.5);
    ryij = ryij - floor(ryij + 0.5);
    rzij = rzij - floor(rzij + 0.5);

    factor = (rxij*(vx[nextsphere] - vx[temppartnr]) + ryij*(vy[nextsphere] - vy[temppartnr]) + rzij*(vz[nextsphere] - vz[temppartnr]))/sigsq;

    delvx = -factor*rxij;
    delvy = -factor*ryij;
    delvz = -factor*rzij;

    vx[nextsphere] = vx[nextsphere] + delvx;
    vx[temppartnr] = vx[temppartnr] - delvx;
    vy[nextsphere] = vy[temppartnr] + delvy;
    vy[temppartnr] = vy[temppartnr] - delvy;
    vz[nextsphere] = vz[nextsphere] + delvz;
    vz[temppartnr] = vz[temppartnr] - delvz;

    w = delvx*rxij + delvy*ryij + delvz*rzij;                   // Calculate Virial w
    return w;
}

// Grsort - Set up radial distribution calculations
void grsort(double *rx, double *ry, double *rz, int hist[])
{
    int const n = 500;                                          // Define number of particles
    unsigned int i, j;                                          // Declare integer counters i and j
    double rxi, ryi, rzi;                                       // Declare position variables
    double rxij, ryij, rzij, rijsq, rij;                        // More position variables
    int bin;                                                    // Declare bin variable
    int const maxbin = 7000;                                    // Define maximum number of bins
    double const delr = 0.0001;                                 // Define delr

    for (i = 0; i < (n-1); i++)
    {
        rxi = rx[i];
        ryi = ry[i];
        rzi = rz[i];

        for (j = (i+1); j < n; j++)
        {
            rxij = rxi - rx[j];
            ryij = ryi - ry[j];
            rzij = rzi - rz[j];
            //rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0) + (0.5*copysign(1.0,rxij)))));
			//ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0) + (0.5*copysign(1.0,ryij)))));
			//rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0) + (0.5*copysign(1.0,rzij)))));
            rxij = rxij - floor(rxij + 0.5);
            ryij = ryij - floor(ryij + 0.5);
            rzij = rzij - floor(rzij + 0.5);

            // Calculate minimum image distances
            rijsq = pow(rxij, 2.0) + pow(ryij, 2.0) + pow(rzij, 2.0);
            rij = sqrt(rijsq);
            bin = (int)(rij/delr);

            if (bin < maxbin)
            {
                //cout << "i = " << i << endl;
                //cout << "j = " << j << endl;

				//cout << "bin = " << bin << endl;
				if (i == 0 && j == 359)
                {

                    cout << "rij = " << rij << endl;

                    cout << "rxi address: " << &rx << endl;
                    cout << "ryi address: " << &ry << endl;
                    cout << "rzi address: " << &rz << endl;
                }
				hist[bin] = hist[bin] + 2;
				//cout << "hist[bin] = " << hist[bin] << endl;
            }
        }
    }
}

int main()
{
    //double eta;                                               // Declare eta (packing fraction) variable
    //srand((unsigned int)time(NULL));                            // Seed random variable with clock time
    srand(10);

    // Prompt and receive inputs
    cout << "Program Hard Sphere" << endl;
    cout << "Molecular Dynamics of Hard Spheres" << endl;
    cout << "Results in Units kt = sigma = 1" << endl;

    cout << "Enter title of run: ";
    cin >> title;

    //cout << "Enter packing fraction (V/V0): ";                // If using papers that report packing fraction variation,
    //cin >> eta;                                               // Use these two lines of input plus the density conversion below.
    cout << "Enter reduced density [(n/v)*(sigma^3)]: ";        // If using papers that report density variation,
    cin >> density;                                             // Use these two lines of input.

    cout << "Enter initial reduced temperature: ";              // T_red = T/T_crit, generally, but this changes for different molecules.
    cin >> temp;                                                // Can use something around 1.5, won't make too much of difference.

    //density = (6*eta)/pi;                                     // Convert between volume fraction with sphere size?
    //density = (sqrt(2))/eta;                                  // Convert between volume fraction with (V0 = close packing) and Bannerman, Lue, Woodcock paper
    //cout << "Reduced density = " << density << endl;          // Display reduced density value

    cout << "Enter number of collisions required: ";            // Generally, the higher the density, the more collisions required.
    cin >> ncoll;                                               // Use benchmark of density = 0.1 requiring about 1,000 collisions.

    // Set up initial parameters
    //temp = floor(double(ncoll)/200.0);                        // This initialization of temperature did not work.
    sigma = pow((density/(double)n), 1.0/3.0);                  // Calculate sigma for kinetic energy calculations.

    // Build initial configuration of system
    createcoord(rx, ry, rz);                                    // Builds FCC lattice of particles with specific position coordinates for each particle.
    createvel(temp, vx, vy, vz);                                // Randomly assigns velocities to each particle.

    // Create file to look at initial configuration             // Received from Ryan Maloney
    FILE *writeout;
    char buf[1000];
    sprintf(buf, "initialconfig");
    writeout = fopen(buf, "w+");
    fprintf(writeout, "MODEL%*d\n", 9, 0);
    fprintf(writeout, "CRYST1%*.3f%*.3f%*.3f%*.2f%*.2f%*.2f%*s%*d\n", 9, 100.0, 9, 100.0, 9, 100.0, 7, 90.00, 7, 90.00, 7, 90.00, 11, "P", 4, 1);
    writePeriodicCoordinates(writeout, rx, ry, rz);
    fprintf(writeout, "ENDMDL\n");

    // Check for particle overlaps
    ovrlap = check(sigma, rx, ry, rz, vx, vy, vz);
    if (ovrlap == 1)
    {
        cout << "Particle overlap in initial configuration" << endl;
        return 0;
    }

    // Calculate initial energy of system
    e = energy(vx, vy, vz);
    cout << "Starting energy = " << e << endl;

    kecntr = 50;
    en = e/((double)n);
    temp = 2.0*en/3.0;
    enkt = en/temp;
    cout << "Initial Temperature = " << temp << endl;
    cout << "Initial e/nkt = " << enkt << endl;

    // Set up initial collision lists hist, coltim, and partnr
    for (i = 0; i < maxbin; ++i)
    {
        hist[i] = 0;
    }

    for (i = 0; i < n; ++i)
    {
        coltim[i] = timbig;
        partnr[i] = n;
    }

    for (j = 0; j < n; ++j)
    {
        uplist(sigma, j, rx, ry, rz, vx, vy, vz, &(coltim[j]), &(partnr[j]));
    }

    //cout << "After 1st Uplist, i and coltim(i): " << i << ", " << coltim[i] << endl;

    // Set up ghost particle to prompt radial distribution calculation
    coltim[n] = timbig;

    // Zero Virial Accumulator
    acw = 0.0;

    // Start of Dynamics
    cout << "**Starting Dynamics...**" << endl;

    // Main dynamics loop
    steps = 0;
    t = 0.0;
    uplistcntr = 10;
    grcount = 20;                                               // Increment by 20
	int nextsphere, temppartnr;									// Use these to distinguish counting spheres from other things

    for (coll = 0; coll < ncoll; ++coll)
    {
        // Locate minimum collision time
        tij = timbig;
        //cout << "Here1" << endl;
        for (k = 0; k < (n+1); ++k)
        {
            //cout << "Here2" << endl;
            if (coltim[k] < tij)
            {
                //cout << "Here3" << endl;
                tij = coltim[k];
                nextsphere = k;
            }
        }
        // If collision time becomes negative, this is incorrect
        if (coltim[nextsphere] < -0.000001)
        {
            //cout << "Here4" << endl;
            cout << "Negative collision time" << endl;
            cout << "Collision: " << coll << endl;
            cout << "Sphere Pair: " << nextsphere << ", " << partnr[nextsphere] << endl;
            cout << "Collision time: " << coltim[nextsphere] << endl;
            exit(0);
        }

    if ((grcount == coll) || (nextsphere == n))
    {
        //cout << "Here5" << endl;
        if (grcount <= coll)
        {
            //cout << "Here6" << endl;
            grsort(rx, ry, rz, hist);
            grcount = grcount + 20;
        }
        // Move particles forward by time tij and reduce collision times
        // Apply periodic boundary conditions
        steps = steps + 1;
        if (nextsphere == n)
        {
            //cout << "Here7" << endl;
            cout << "Next Sphere was " << coll << endl;
            coltim[n] = timbig;
            t = t + tij;
            for (k = 0; k < n; ++k)
            {
                //cout << "Here8" << endl;
                coltim[k] = coltim[k] - tij;
                rx[k] = rx[k] + vx[k]*tij;
                ry[k] = ry[k] + vy[k]*tij;
                rz[k] = rz[k] + vz[k]*tij;
                rx[k] = rx[k] - floor(rx[k] + 0.5);
                ry[k] = ry[k] - floor(ry[k] + 0.5);
                rz[k] = rz[k] - floor(rz[k] + 0.5);
            }

            for (i = 0; i < n; ++i)
            {
                //cout << "Here9" << endl;
                uplist(sigma, i, rx, ry, rz, vx, vy, vz, &(coltim[i]), &(partnr[i]));
            }
        }
    }
    else
    {
        //cout << "Here10" << endl;
        temppartnr = partnr[nextsphere];
        t = t + tij;
        for (k = 0; k < n; ++k)
        {
            //cout << "Here11" << endl;
            coltim[k] = coltim[k] - tij;
            rx[k] = rx[k] + vx[k]*tij;
            ry[k] = ry[k] + vy[k]*tij;
            rz[k] = rz[k] + vz[k]*tij;
            rx[k] = rx[k] - floor(rx[k] + 0.5);
            ry[k] = ry[k] - floor(ry[k] + 0.5);
            rz[k] = rz[k] - floor(rz[k] + 0.5);
        }

        //cout << "Here12" << endl;
        coltim[n] = coltim[n] - tij;

        // Compute collision dynamics
        w = bump(sigma, nextsphere, temppartnr, rx, ry, rz, vx, vy, vz);

        acw = acw + w;

        // Reset collision lists for those particles that need it
        //cout << "Here13" << endl;
        for (int l = 0; l < n; ++l)
        {
            //cout << "Here14" << endl
            if ((l == nextsphere) || (partnr[l] == nextsphere) || (l = temppartnr) || (partnr[l] == temppartnr))
            {
                //cout << "Here15" << endl;
                uplist(sigma, l, rx, ry, rz, vx, vy, vz, &(coltim[l]), &(partnr[l]));
            }
        }

        //cout << "Here16" << endl;
        dnlist(sigma, nextsphere, rx, ry, rz, vx, vy, vz, coltim, partnr);
        //cout << "Here17" << endl;
        dnlist(sigma, temppartnr, rx, ry, rz, vx, vy, vz, coltim, partnr);
    }

    //if (coll == uplistcntr-1)
    //{
    //    for (i = 1; i <= n; ++i)
    //    {
    //        uplist(sigma, i, rx, ry, rz, vx, vy, vz, coltim[i], partnr[i]);
    //    }
    //    uplistcntr = (uplistcntr + 10);
    //}

    if ((coll+1) == kecntr)
    {
        //cout << "Here18" << endl;
        for (i = 0; i < n; ++i)
            {
                //cout << "Here19" << endl;
                e = e + pow(vx[i], 2.0) + 2*pow(vy[i], 2.0) + 2*pow(vz[i], 2.0); // Sum up kinetic energy of system
            }
        e = 0.5*e;
        en = e/double(n);
        enkt = en/temp;
        kecntr = kecntr + 50;
    }
    //cout << "Here20" << endl;
    }

    // End main loop and dynamics
    cout << "**Ending dynamics...**" << endl;;
    cout << "Final Colliding Pair: " << nextsphere << ", " << temppartnr << endl;

    cout << "Number of collisions = " << coll << endl;
    cout << "Kinetic energy = " << e << endl;
    cout << "Enkt = " << enkt << endl;
    cout << "Virial w = " << w << endl;

    // Check for particle overlaps
    ovrlap = check(sigma, rx, ry, rz, vx, vy, vz);

    if (ovrlap == 1.0)
    {
        cout << "Particle overlap in final configuration" << endl;
    }

    // Write out configuration
    ofstream file1;
    file1.open ("hardsphere_endconfig.out");
    file1 << "Run Title: " << title << endl;
    file1.precision(4);
    file1 << setw(20) << "r(x)" << setw(20) << "r(y)" << setw(20) << "r(z)" << setw(20) << "v(x)" << setw(20) << "v(y)" << setw(20) << "v(z)" << endl;
    for (i = 1; i <= n; ++i)
    {
        file1 << setw(20) << rx[i] << setw(20) << ry[i] << setw(20) << rz[i] << setw(20) << vx[i] << setw(20) << vy[i] << setw(20) << vz[i] << endl;
    }
    file1.close();

    // Write out information of interest
    e = energy(vx, vy, vz);
    pvnkt1 = acw/(double(n)*3.0*t*temp) + 1;
    en = e/((double)n);
    temp = 2.0 * en / 3.0;
    enkt = en/temp;
    t = t*sqrt(temp)/sigma;
    rate = double(ncoll)/t;
    tbc = double(n)/rate/2.0;

    cout << "The final sample size n = " << n << endl;
    cout << "The final kinetic energy e = " << e << endl;
    cout << "The final temperature temp = " << temp << endl;
    cout << "The final Virial accumulation acw = " << acw << endl;
    cout << "The final time t = " << t << endl;
    cout << "The mean collision rate = " << rate << endl;
    cout << "The mean collision time = " << tbc << endl;
    cout << "The final e/nkt = " << enkt << endl;
    cout << "Z = PV/nkt = " << pvnkt1 << endl;
    cout << "The final grcount = " << grcount << endl;

    // Calculate g(r) - reference: pg. 184 in Allen and Tildesley

    grconst = 4.0*pi*n/3.0;

    for (bin = 0; bin < maxbin; ++bin)
    {
        rlower = double(bin)*delr;
        rupper = rlower + delr;
        nideal = grconst*(pow(rupper, 3.0) - (pow(rlower, 3.0)));
        gr[bin] = ((double)hist[bin])/((double)steps)/((double)n)/nideal;
        f[bin] = rlower + delr/2.0;
    }

    // Write radial distribution file
    ofstream file2;
    file2.open ("gr_radial.out");
    file2 << setw(30) << "f(bin)/sigma" << setw(30) << "gr(bin)" << endl;
    for (bin = 1; bin <= maxbin; ++bin)
    {
        file2 << setw(30) << f[bin]/sigma << setw(30) << gr[bin] << endl;
    }
    file2.close();

    // Write out information of interest to external file
    ofstream file3;
    file3.open ("output.out");
    file3.precision(4);
    file3 << setw(30) << "The final sample size n = " << setw(30) << n << endl;
    file3 << setw(30) << "The final kinetic energy e = " << setw(30) << e << endl;
    file3 << setw(30) << "The final temperature = " << setw(30) << temp << endl;
    file3 << setw(30) << "The final Virial accumulation acw = " << setw(30) << acw << endl;
    file3 << setw(30) << "The final time = " << setw(30) << t << endl;
    file3 << setw(30) << "The mean collision rate = " << setw(30) << rate << endl;
    file3 << setw(30) << "The mean collision time = " << setw(30) << tbc << endl;
    file3 << setw(30) << "The final e/nkt = " << setw(30) << enkt << endl;
    file3 << setw(30) << "Z = PV/nkt = " << setw(30) << pvnkt1 << endl;
    file3 << setw(30) << "The final grcount = " << setw(30) << grcount << endl;
    file3.close();

    return 0;
    system("PAUSE");
    }
