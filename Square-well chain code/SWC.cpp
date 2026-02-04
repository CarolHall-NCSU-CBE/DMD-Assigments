#include "swchain.h"

int main()
{
    srand((unsigned int)time(NULL));    //seed random function
//    srand48(1);

//Variable declaration
    //inputs
    double *rx, *ry, *rz, *vx, *vy, *vz, *coltim;
	double density, temperature, ctemp, sigma1, sigma2, d1, d2, rate, tbc;
    double acw, t, ke, pe, en, totalenergy, enkt, tij, w, grconst, burnintime;
    int ncoll, chk;
	int *coltyp, *partnr;
    char fname [] = "chain";
    int down, up, i, j, bin;
    int coll, burnin, uplistcntr, grcount, steps, ghostcntr, ghost;
    int hist[maxbin];
    double rlower, rupper, nideal, gr[maxbin], f[maxbin], pvnkt, z;
	int **linkedspheres;

    FILE *outfile;
    FILE *grfile;
	FILE *compressibility, *tt;

	rx = (double*)malloc((n+1)*sizeof(double));
	if (rx == NULL)
	{
		printf("Out of memory: rx\n");
		exit(0);
	}
	ry = (double*)malloc((n+1)*sizeof(double));
	if (ry == NULL)
	{
		printf("Out of memory: ry\n");
		exit(0);
	}
	rz = (double*)malloc((n+1)*sizeof(double));
	if (rz == NULL)
	{
		printf("Out of memory: rz\n");
		exit(0);
	}
	vx = (double*)malloc((n+1)*sizeof(double));
	if (vx == NULL)
	{
		printf("Out of memory: vx\n");
		exit(0);
	}
	vy = (double*)malloc((n+1)*sizeof(double));
	if (vy == NULL)
	{
		printf("Out of memory: vy\n");
		exit(0);
	}
	vz = (double*)malloc((n+1)*sizeof(double));
	if (vz == NULL)
	{
		printf("Out of memory: vz\n");
		exit(0);
	}
	coltim = (double*)malloc((n+2)*sizeof(double));
	if (coltim == NULL)
	{
		printf("Out of memory: coltim\n");
		exit(1);
	}
	coltyp = (int*)malloc((n+1)*sizeof(int));
	if (coltyp == NULL)
	{
		printf("Out of memory: coltyp\n");
		exit(0);
	}
	partnr = (int*)malloc((n+1)*sizeof(int));
	if (partnr == NULL)
	{
		printf("Out of memory: partnr\n");
		exit(0);
	}
	linkedspheres = (int**)malloc((n+1)*sizeof(int*));
	if (linkedspheres == NULL)
	{
		printf("Out of memory: linkedspheres\n");
		exit(1);
	}		
	for(i=0;i<n;i++)
	{
		linkedspheres[i] = (int*)malloc(numoflinks*sizeof(int));
		for(j=0;j<numoflinks;j++)
		{	
			linkedspheres[i][j] = -1;
		}
	}
//Read in basic simulation parameters
    printf("PROGRAM SPHERE\nMolecular Dynamics of Square Well Chains\n");
    printf("Enter reduced density (N/V)*sigma^3: ");
    scanf("%lf", &density);
 

    printf("Enter desired reduced temp: ");
    scanf("%lf", &temperature);

//    printf("Enter number of collisions required: ");
//    scanf("%d", &ncoll);
	ncoll = 20000000;
//    printf("Enter configuration filename: ");
//    scanf("%s", fname);
	
//Create initial configuration
    sigma1 = pow(density/((double)n), 1.0/3.0);
	d1 = sigma1 * (1-del);
	d2 = sigma1 * (1+del);
	sigma2 = sigma1*lambda;
	if(sigma2 >= 0.5)
	{
		printf("Square well diameter too large, artifacts may develop\n");
		exit(1);
	}
	initpos(sigma1, d1, d2, rx, ry, rz, linkedspheres);  
    initvel(500.00, vx, vy, vz);              
	
	chk = check(sigma1, d1, d2, rx, ry, rz, linkedspheres);
    if(chk == 1)
    {
        printf("Overlap in initial configuration");
        exit(1);
    }
	
//Calculate energy and initial temp
    ke = 0.0;
    ke = kecalc(vx, vy, vz);
	pe = pecalc(sigma2, rx, ry, rz, linkedspheres);
    totalenergy = ke + pe;
	en = ke/(double)n;
    ctemp = 2.0*en/3.0;
    enkt = en/ctemp;
    //initialize system variables
    
	acw = 0.0;                              //Virial Accumulator
    t = 0.0;        						//length of simulation
    uplistcntr = 50;                       //used to periodically recall uplist
    grcount = 4000000;                           //used to periodically calculate g(r)
    steps = 0;
    ghostcntr = 0;
	burnin = 3000000;
	int burnin2 = 6000000;
	int burnin3 = 9000000;
	int burnin4	= 15000000;
	double avgke = 0;
	int avgkecntr = 0;
	int nextghostcoll = burnin + round(ranf() * 25) + 25;
	int lastcollision = 200;
	int lastpartner = 200;

	for(i=0;i<maxbin;i++)
        hist[i] = 0;
    printf("\nReduced density: \t\t%.5lf\nReduced Temp: \t\t\t%.5lf\nCollisions required: \t\t%d\nConfiguration filename: \t%s\n",density,ctemp,ncoll,fname);
	
	//Seed initial collision lists
    for(i=0;i<n;i++)
    {
		coltim[i] = TIMBIG;
        partnr[i] = n;
		coltyp[i] = 0;
    }
    coltim[n] = TIMBIG;
//Create initial collision lists
    for(i=0;i<n;i++)
    {
        uplist(i, sigma1, sigma2, d1, d2, rx, ry, rz, vx, vy, vz, &coltim[i], &partnr[i], &coltyp[i], linkedspheres[i]);
    }
//Create output files
    outfile = fopen(fname,"w");
    if(!outfile)
    {
        puts("outfile file error!");
        exit(1);
    }
    fprintf(outfile,"DMD Hard Sphere\nDensity:\t\t\t  %lf\nNumber of Collisions: %d\nSigma1:\t\t\t\t  %lf\nSigma2:\t\t\t\t  %lf\n",density,ncoll,sigma1,sigma2);
    grfile = fopen("g(r).txt","w");
    if(!grfile)
	{
		puts("gr file error");
		exit(1);
	}
	compressibility = fopen("Compressibility factor.txt","w");
    if(!compressibility)
    {
        puts("compressibility file error!");
        exit(1);
    }
    fprintf(outfile,"Initial temp:\t\t\t%lf\nInitial ke:\t\t\t%lf\nInitial pe:\t\t\t%lf\nInitial en\t\t\t%lf\nInitial enkt:\t\t\t%lf\n", ctemp,ke,pe,en,enkt);
	
	FILE *energyrecord = fopen("Energy Record.txt", "w");
	if(!energyrecord)
    {
        puts("Energy Record file error!");
        exit(1);
    }
	fprintf(energyrecord,"Time:\nEnergy:\n");
	tt = fopen("tempvstime.txt", "w");
	if(!tt)
    {
        puts("tt file error!");
        exit(1);
    }
	FILE *gc = fopen("ghost collisions","w");
	if(!gc)
    {
        puts("gc file error!");
        exit(1);
    }


int negcolcntr = 0;
//*************Start of Dynamics*********************************
    for(coll=0;coll<ncoll;coll++)
    {
		if(coll == burnin4)
		{
			acw = 0;
			burnintime = t;
		}
		if(coll == burnin3)
		{
			printf("Burn in period complete at %d collisions\n", coll);
			initvel(temperature, vx, vy, vz);
			for(int i = 0; i < n; i++)
			{
				coltim[i] = TIMBIG;
				partnr[i] = n;
			}
			for(int j = 0; j < n; j++)
			{
				uplist(j, sigma1, sigma2, d1, d2, rx, ry, rz, vx, vy, vz, &(coltim[j]), &(partnr[j]), &(coltyp[j]), linkedspheres[j]);
			}
			// reseting variable for the main part of the loop after burnin
			coltim[n] = 5.0;
		}
		if(coll == burnin2)
		{
			printf("Burn in period complete at %d collisions\n", coll);
			initvel(10.0, vx, vy, vz);
			for(int i = 0; i < n; i++)
			{
				coltim[i] = TIMBIG;
				partnr[i] = n;
			}
			for(int j = 0; j < n; j++)
			{
				uplist(j, sigma1, sigma2, d1, d2, rx, ry, rz, vx, vy, vz, &(coltim[j]), &(partnr[j]), &(coltyp[j]), linkedspheres[j]);
			}
			// reseting variable for the main part of the loop after burnin
			coltim[n] = 5.0;
		}  
        if(coll == burnin)
		{
			printf("Burn in period complete at %d collisions\n", coll);
			initvel(100.00, vx, vy, vz);
			for(int i = 0; i < n; i++)
			{
				coltim[i] = TIMBIG;
				partnr[i] = n;
			}
			for(int j = 0; j < n; j++)
			{
				uplist(j, sigma1, sigma2, d1, d2, rx, ry, rz, vx, vy, vz, &(coltim[j]), &(partnr[j]), &(coltyp[j]), linkedspheres[j]);
			}
			// reseting variable for the main part of the loop after burnin
			coltim[n] = 5.0;
		}  
		
	//locate minimum collision time
		tij =  TIMBIG;   //seed tij
        for(i=0;i<n+1;i++)
        {
            if(coltim[i] < tij)     
            {
                tij = coltim[i];    
                down = i;     
            }                      
        }
	//look for negative collisions
		if(tij <= 0)
		{
            up = partnr[down];
			fprintf(gc,"negative collision!\n********Collision %d***********\nColtim: %lf\ni: %d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\nj: %d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\nColltype: %d\nSigma1: %lf\n",coll, coltim[down],down,rx[down],ry[down],rz[down],vx[down],vy[down],vz[down],up,rx[up],ry[up],rz[up],vx[up],vy[up],vz[up],coltyp[down],sigma1);
			negcolcntr++;
			printf("NEGATIVE Collision");
	        exit(1);
        }
		
	//radial distribution function calculations
		if(grcount == coll)
        {
            grsort(rx, ry, rz, hist);
			grcount = grcount + 20;
			coltim[n] = t + 5.0;
			steps = steps + 1;
		}              
		if(down == n-1)
        {
            printf("next collision sphere was n, coll number %d\n",coll);
            coltim[n] = TIMBIG;
            t = t + tij;
            for(i=0;i<n;i++) //move particles, adjust coltim for each particle
            {
                coltim[i] = coltim[i]-tij;
                rx[i] = rx[i] + vx[i]*tij;
                ry[i] = ry[i] + vy[i]*tij;
                rz[i] = rz[i] + vz[i]*tij;
                rx[i] = rx[i] - (1.0*(double)((int)(rx[i]/(double)(1.0)+(0.5*copysign(1.0,rx[i])))));
                ry[i] = ry[i] - (1.0*(double)((int)(ry[i]/(double)(1.0)+(0.5*copysign(1.0,ry[i])))));
                rz[i] = rz[i] - (1.0*(double)((int)(rz[i]/(double)(1.0)+(0.5*copysign(1.0,rz[i])))));
            }
            for(i=0;i<n;i++)
            {
                uplist(i, sigma1, sigma2, d1, d2, rx, ry, rz, vx, vy, vz, &coltim[i], &partnr[i], &coltyp[i], linkedspheres[i]);
            }
        }
     
	//perform collision
		if(coll == nextghostcoll)
		{
//			fprintf(tt,"%d\n",coll);
			ghost = (int) (ranf()*n);
			while(ghost == lastcollision || ghost == lastpartner)
				ghost = (int) (ranf()*n);
			
			ghostcoll(temperature, &vx[ghost], &vy[ghost], &vz[ghost]);
			ghostcntr++;
			nextghostcoll = round(ranf() * 20) +20 + coll;
			for(i=0;i<n;i++)
				{
					if((i==ghost) || (partnr[i]==ghost))
					{
						uplist(i, sigma1, sigma2, d1, d2, rx, ry, rz, vx, vy, vz, &(coltim[i]), &(partnr[i]), &(coltyp[i]), linkedspheres[i]);
					}
				}
			dnlist(ghost, sigma1, sigma2, d1, d2, rx, ry, rz, vx, vy, vz, coltim, partnr, coltyp, linkedspheres);
//			fprintf(gc, "%d\t%d\n",ghost,coll);		
		} 	
		
		else
		{
			t = t + tij;		
			for(i=0;i<n;i++)
			{
				coltim[i] = coltim[i] - tij;    //reduce collision times
				rx[i] = rx[i] + vx[i]*tij;      //move each particle
				ry[i] = ry[i] + vy[i]*tij;
				rz[i] = rz[i] + vz[i]*tij;
				rx[i] = rx[i] - (1.0*(double)((int)(rx[i]/(double)(1.0)+(0.5*copysign(1.0,rx[i])))));
				ry[i] = ry[i] - (1.0*(double)((int)(ry[i]/(double)(1.0)+(0.5*copysign(1.0,ry[i])))));
				rz[i] = rz[i] - (1.0*(double)((int)(rz[i]/(double)(1.0)+(0.5*copysign(1.0,rz[i])))));
			}
			if(down < n)
			{
				up = partnr[down];
				w = bump(sigma1, sigma2, d1, d2, down, up, rx, ry, rz, vx, vy, vz, &coltyp[down], linkedspheres[down]);                  
				acw = acw + w;
				lastcollision = down;
				lastpartner = up;
				for(i=0;i<n;i++)
				{
					if((i==down) || (partnr[i]==down) || (i==up) || (partnr[i]==up))
					{
						uplist(i, sigma1, sigma2, d1, d2, rx, ry, rz, vx, vy, vz, &coltim[i], &partnr[i], &coltyp[i], linkedspheres[i]);
					}
				}
				dnlist(down, sigma1, sigma2, d1, d2, rx, ry, rz, vx, vy, vz, coltim, partnr, coltyp, linkedspheres);
				dnlist(up, sigma1, sigma2, d1, d2, rx, ry, rz, vx, vy, vz, coltim, partnr, coltyp, linkedspheres);
			}
		}
	
	
	//reset coltim after 100 collisions
		if(coll == (uplistcntr-1))
		{
			for(i=0;i<(n-1);i++)
			{
				uplist(i, sigma1, sigma2, d1, d2, rx, ry, rz, vx, vy, vz, &coltim[i], &partnr[i], &coltyp[i], linkedspheres[i]);
			}
			uplistcntr+=1000;
		}

		if(coll >= (ncoll-2000000))
		{
			ke = kecalc(vx,vy,vz);
			avgke += ke;
			avgkecntr++;
			ctemp = 2.0*avgke/avgkecntr/3.0/(double)n;
			z = chnlen + (acw/(t-burnintime)/(double)(numchn)/3.0/ctemp);
			fprintf(compressibility,"%lf\n",z);
		
		}
		if(coll % 10000 == 0)
		{
			ke = kecalc(vx,vy,vz);
			ctemp = 2.0*ke/3.0/(double)n;
			fprintf(tt,"%lf\t%lf\n",t,ctemp);
			pe = pecalc(sigma2, rx, ry, rz, linkedspheres);
			totalenergy = pe + ke;
			fprintf(energyrecord,"%lf\t%lf\n",t,totalenergy);
		}
		
		if(coll % 1000000 == 0)
		{
			printf("percent complete: %lf\n", (double)coll/ncoll);
		}   
	}
//*************END OF DYNAMICS**************************************************
	
	
	fclose(gc);
/*	FILE *vfinal = fopen("final velocity","w");
	if(!vfinal)
    {
        puts("vfinal file error!");
        exit(1);
    }
	
	
	for(i=0;i<n;i++)
	{
		fprintf(vfinal,"%lf\t%lf\t%lf\n",vx[i],vy[i],vz[i]);
		fprintf(finalconfig, "%lf\t%lf\t%lf\n", rx[i],ry[i],rz[i]);
	}	
	fclose(finalconfig);
	fclose(vfinal);
*/		
	
	chk = check(sigma1, d1, d2, rx, ry, rz, linkedspheres);
    if(chk == 1)
    {
        printf("Overlap in final configuration\n");
        exit(1);
    }
	
	double lastpvnkt;
    //write out interesting information
//    ke = kecalc(vx,vy,vz);
	avgke = avgke/(double) avgkecntr;
	printf("avgke: %lf\navgkecntr: %d\n",avgke,avgkecntr);
	pe = pecalc(sigma2, rx, ry, rz, linkedspheres);
    ctemp = 2.0*avgke/3.0/(double)n;
	lastpvnkt = (acw / ((double)n*3.0*(t-burnintime)*ctemp));
	z = chnlen + (acw/(t-burnintime)/(double)(numchn)/3.0/ctemp);
	totalenergy = ke + pe;
    enkt = ke/(double)n/ctemp;
    t = (t - burnintime)*sqrt(ctemp)/sigma1;
    rate = ((double)ncoll)/t;
    tbc = ((double)n)/rate/2.0;

    fprintf(outfile,"\nFinal n:\t%d\nFinal ke:\t%lf\nFinal pe:\t%lf\nFinal en:\t%lf\nFinal temp:\t%lf\nFinal acw:\t%lf\nFinal time:\t%lf\nCollision rate:\t%lf\n",n,ke,pe,totalenergy,ctemp,acw,t,rate);
    fprintf(outfile,"Mean collision time:\t%lf\nFinal e/nkt:\t%lf\nZ:\t%lf\nFinal grcount:\t%d\nNumber Ghost Coll: \t%d\nNEGCOLL:\t%d\n",tbc,enkt,z,grcount,ghostcntr,negcolcntr);

    //calculate g(r) from pg. 184 in Allen and Tildesley
    grconst = 4.0*pi*n/3.0;
        for(bin=0;bin<maxbin;bin++)
        {
            rlower = ((double)bin)*delr;
            rupper = rlower + delr;
            nideal = grconst*(pow(rupper, 3.0)-pow(rlower, 3.0));
            gr[bin] = (double) hist[bin]/((double) steps)/((double) n)/nideal;
            f[bin] = rlower + delr/2;
            fprintf(grfile, "%lf\t%lf\n", f[bin]/sigma1, gr[bin]);
        }

	fclose(grfile);
    fclose(outfile);
	fclose(compressibility);
	fclose(energyrecord);
	fclose(tt);
	free(rx);
	free(ry);
	free(rz);
	free(vx);
	free(vy);
	free(vz);
	free(coltim);
	free(partnr);
	free(coltyp);
	return 0;
}

//initpos builds chains one at a time by placing a bead, then attaching additional beads until the chain length is met, checking for overlaps at each bead placement

void initpos(double sigma1, double d1, double d2, double rx[], double ry[], double rz[], int **linkedspheres)
{
	double chaindelta, newx, newy, newz, xdir, ydir, nextlayer, tempx, tempy;
	int lengthofarray, chainstart, i, j, l;
	chaindelta = (sigma1 + d2) / 2.0;
	newx = chaindelta - 0.5;
	newy = -chaindelta + 0.5;
	newz = chaindelta - 0.5;
	xdir = 1.0;
	ydir = -1.0;
	lengthofarray = sizeof(linkedspheres[0]);
	//Get the value of the first position of each chain
	for(i = 0; i < numchn; i++)
	{
		chainstart = i * chnlen;
		rx[chainstart] = newx;
		ry[chainstart] = newy;
		rz[chainstart] = newz;
		for(j = 1; j < chnlen; j++)
		{
			for(l = 0; l < lengthofarray; l++)
			{
				if(linkedspheres[i * chnlen + j - 1][l] == -1)
				{
					linkedspheres[i * chnlen + j - 1][l] = i * chnlen + j;
					break;
				}
			}
			for(l = 0; l < lengthofarray; l++)
			{
				if(linkedspheres[i * chnlen + j][l] == -1)
				{
					linkedspheres[i * chnlen + j][l] = i * chnlen + j - 1;
					break;
				}
			}
			
			newx = rx[chainstart + j - 1] + chaindelta * xdir;
			if((float)(newx * xdir) <= (float)(0.5 - chaindelta))  ///put spheres on x dirc, just less than the boundary
			{
				rx[j + chainstart] = newx;
				ry[j + chainstart] = ry[chainstart + j - 1];
				rz[j + chainstart] = rz[chainstart + j - 1];
			}
			else
			{
				newy = ry[chainstart + j - 1] + chaindelta * ydir;
				if((float)(newy * ydir) <= (float)(0.5 - chaindelta))
				{
					xdir = xdir * -1.0;
					rx[j + chainstart] = rx[chainstart + j - 1];
					ry[j + chainstart] = newy;
					rz[j + chainstart] = rz[chainstart + j - 1];
				}
				else
				{
					nextlayer = sqrt(0.5 * pow(chaindelta,2.0));
					tempx = 0.5 * chaindelta * xdir;
					tempy = 0.5 * chaindelta * ydir;
					xdir = xdir * -1.0;
					ydir = ydir * -1.0;
					if((tempx * -1.0 * xdir) > (float)(0.5 - chaindelta))///错误判断
						tempx = 0.5 * chaindelta * xdir;
					if((tempy * -1.0 * ydir) > (float)(0.5 - chaindelta))
						tempy = 0.5 * chaindelta * ydir;
					rx[j + chainstart] = rx[chainstart + j - 1] + tempx;
					ry[j + chainstart] = ry[chainstart + j - 1] + tempy;
					rz[j + chainstart] = rz[chainstart + j - 1] + nextlayer;
					if(rz[chainstart + j - 1] + nextlayer > 0.5)
					{
						printf("Too many chains for the density required. Finished chain %d\n", i-1);
						exit(1);
					}
				}	
			}
		}

		newx = rx[chainstart + chnlen - 1] + chaindelta * xdir;
		if((float)(newx*xdir) <= (float)(0.5 - chaindelta))
		{
			newy = ry[chainstart + chnlen - 1];
			newz = rz[chainstart + chnlen - 1];
		}
		else
		{
			newy = ry[chainstart + chnlen - 1] + chaindelta * ydir;
			if((float)(newy*ydir) <= (float)(0.5 - chaindelta))
			{
				xdir = xdir * -1.0;
				newx = rx[chainstart + chnlen - 1];
				newy = newy;
				newz = rz[chainstart + chnlen - 1];
			}
			else
			{
				nextlayer = sqrt(0.5 * pow(chaindelta,2.0));
				tempx = 0.5 * chaindelta * xdir;
				tempy = 0.5 * chaindelta * ydir;
				xdir = xdir * -1.0;
				ydir = ydir * -1.0;
				if((tempx * -1.0 * xdir) > (float)(0.5 - chaindelta))
					tempx = 0.5 * chaindelta * xdir;
				if((tempy * -1.0 * ydir) > (float)(0.5 - chaindelta))
					tempy = 0.5 * chaindelta * ydir;
				newx = rx[chainstart + chnlen - 1] + tempx;
				newy = ry[chainstart + chnlen - 1] + tempy;
				newz = rz[chainstart + chnlen - 1] + nextlayer;
				if(newz > 0.5)
				{
					printf("Too many chains for the density required. Finished chain %d\n", i-1);
					exit(1);
				}
			
			}
		}
	}
/*	FILE *initialconfig;
	initialconfig = fopen("initialconfig", "w");
	if(!initialconfig)
    {
        puts("configuration file error!");
        exit(1);
    }
	for(i=0;i<n;i++)
		fprintf(initialconfig, "%lf\t%lf\t%lf\n", rx[i],ry[i],rz[i]);
	for(i=0;i<n;i++)
	{
		for(l=0;l<lengthofarray;l++)
		{
			if(linkedspheres[i][l] != -1)
			{
				fprintf(initialconfig, "%*d", 5, linkedspheres[i][l]);
			}
			else
			{
				break;
			}
		}
		fprintf(initialconfig,"\n");
	}
	fclose(initialconfig);				*/
}

//Initial velocities from maxwell-boltzmann distribution
//determined by temp and mass
void initvel(double temperature, double vx[], double vy[], double vz[])
{
	double rtemp, sumx, sumy, sumz;
    int i;
    rtemp = sqrt(temperature);
    sumx = 0.0;
    sumy = 0.0;
    sumz = 0.0;

    for(i=0;i<n;i++)
    {
        vx[i] = rtemp * gauss();
        vy[i] = rtemp * gauss();
        vz[i] = rtemp * gauss();
        sumx = sumx + vx[i];
        sumy = sumy + vy[i];
        sumz = sumz + vz[i];
    }
//remove net momentum
    sumx = sumx / (double)n;
    sumy = sumy / (double)n;
    sumz = sumz / (double)n;

    for(i=0;i<n;i++)
    {
        vx[i] = vx[i] - sumx;
        vy[i] = vy[i] - sumy;
        vz[i] = vz[i] - sumz;
    }
 }

int check(double sig, double d1, double d2, double rx[], double ry[], double rz[], int**linkedspheres)
{
    int i, j, linked, loopkvar;
    double rxi, ryi, rzi, rxij, ryij, rzij, rij, rijsq, sigsq;
    int chk = 0;

	sigsq = sig*sig;
    for(i=0;i<(n-1);i++)
    {
        rxi = rx[i];        //get the ith position for comparison
		ryi = ry[i];
		rzi = rz[i];
		for(j=i+1;j<n;j++)
		{
			linked = 0;
			loopkvar = 0;
			while(loopkvar < sizeof(linkedspheres) && linkedspheres[i][loopkvar] != -1)
			{
				if(linkedspheres[i][loopkvar] == j)
				{
					linked = 1;
					break;
				}
				loopkvar += 1;
			}
			rxij = rxi - rx[j];     
			ryij = ryi - ry[j];    
			rzij = rzi - rz[j];
			rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0)+(0.5*copysign(1.0,rxij)))));
			ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0)+(0.5*copysign(1.0,ryij)))));
			rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0)+(0.5*copysign(1.0,rzij)))));
			rijsq = pow(rxij, 2.0) + pow(ryij, 2.0) + pow(rzij, 2.0);
			if(linked == 1)
			{
				if(rijsq - pow(d1,2.0) < -0.000000000001)
				{
					printf("beads in chain are overlapped!\n");
					printf("r%d %d: %f, tol: %.16f\n", i,j,rijsq, rijsq - pow(d1,2.0));
					printf("%lf\t%lf\t%lf\n",rx[i],ry[i],rz[i]);
					printf("%lf\t%lf\t%lf\n",rx[j],ry[j],rz[j]);
					chk = 1;
				}
				else if(rijsq > pow(d2,2.0))
				{
					printf("chains coming apart!\n");
					printf("r%d %d: %f, dif: %f\n", i,j,rij, TOL);
					chk = 1;
				}	
			}
			else
			{
				if(rijsq < sigsq)
				{
					rij = sqrt(rijsq/sigsq);
					if((1.0-rij) > TOL)
					{
					printf("chain overlap\n %d %d: %f, tol: %f\n", i,j,rij, TOL);
					chk = 1;
					}
				}
			}
		}
	}
    return(chk);
}

//calculate kinetic energy
double kecalc(double vx[], double vy[], double vz[])
{
    double ke = 0.0;
    int d;

    for(d=0;d<n;d++)
        ke = ke + pow(vx[d],2.0) + pow(vy[d],2.0) + pow(vz[d],2.0);
    ke = 0.5 * ke;
    return(ke);
}

//calculate potential energy
double pecalc(double sigma2, double rx[], double ry[], double rz[], int** linkedspheres)
{
    double pe, sig2sq, rxi, ryi, rzi, rxij, ryij, rzij, rijsq;
    int l,k,linked,loopkvar;
	pe = 0.0;
	sig2sq = pow(sigma2, 2.0);

    for(l=0;l<n-1;l++)
	{
		rxi = rx[l];
		ryi = ry[l];
		rzi	= rz[l];
		for(k=l+1;k<n;k++)
		{
			linked = 0;
			loopkvar = 0;
			while(loopkvar < sizeof(linkedspheres) && linkedspheres [l][loopkvar] != -1)
			{
				if(linkedspheres[l][loopkvar] == k)
				{
					linked = 1;
					break;
				}
				loopkvar += 1;
			}
			if(linked == 0)
			{			
				rxij = rxi - rx[k];
				ryij = ryi - ry[k];
				rzij = rzi - rz[k];
				rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0)+(0.5*copysign(1.0,rxij)))));
				ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0)+(0.5*copysign(1.0,ryij)))));
				rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0)+(0.5*copysign(1.0,rzij)))));
				rijsq = pow(rxij, 2.0) + pow(ryij, 2.0) + pow(rzij, 2.0);
				if(rijsq-sig2sq<0.0)
				{
					pe = pe - eps;
				}
			}
		}
	}
	return(pe);
}

void uplist(int i, double sigma1, double sigma2, double d1, double d2, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], double *coltim, int *partnr, int *coltyp, int linkedspheres[])
{
	if(i==(n-1))
        return;
    
	double rxi, ryi, rzi, rxij, ryij, rzij, vxi, vyi, vzi, vxij, vyij, vzij, sig11, sig22, sig1sq, sig2sq, tij, cij1, cij2, bij, bijsq, rijsq, vijsq, disc1, disc2;
    int l, tempcoltyp, linked, loopkvar, chk;
    *coltim = TIMBIG;
    rxi = rx[i];
    ryi = ry[i];
    rzi = rz[i];
    vxi = vx[i];
    vyi = vy[i];
    vzi = vz[i];
    for(l=i+1;l<n;l++)
    {
        rxij = rxi - rx[l];
        ryij = ryi - ry[l];
        rzij = rzi - rz[l];
        rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0)+(0.5*copysign(1.0,rxij)))));
        ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0)+(0.5*copysign(1.0,ryij)))));
        rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0)+(0.5*copysign(1.0,rzij)))));
        vxij = vxi - vx[l];
        vyij = vyi - vy[l];
        vzij = vzi - vz[l];
        bij = rxij * vxij + ryij * vyij + rzij * vzij;
        bijsq = pow(bij, 2.0);
        rijsq = pow(rxij, 2.0) + pow(ryij, 2.0) + pow(rzij, 2.0);
        vijsq = pow(vxij, 2.0) + pow(vyij, 2.0) + pow(vzij, 2.0);
		
		linked = 0;
		loopkvar = 0;
		while(loopkvar < sizeof(linkedspheres) && linkedspheres[loopkvar] != -1)
		{
			if(linkedspheres[loopkvar] == l)
			{
				linked = 1;
				break;
			}
			loopkvar += 1;
		}
		if(linked == 1) 
		{
			sig11 = d1;
			sig22 = d2;
			chk = 0;
		}
		else
		{
			sig11 = sigma1;
			sig22 = sigma2;
			chk = 1;
		}
		sig1sq = sig11 * sig11;
		sig2sq = sig22 * sig22;
		cij1 = rijsq - sig1sq;
		cij2 = rijsq - sig2sq;
		disc1 = bijsq - vijsq * cij1;
		disc2 = bijsq - vijsq * cij2;
		if(bij < 0.0)       //particles approaching
        {
			if(cij2 <= 0.0)		//particles within attractive range, moving towards each other
			{
				if(disc1 > 0.0)		//CORE collision
				{
					tij = (-bij - sqrt(disc1)) / vijsq;
					if(tij < *coltim)
					{
						*coltim = tij;
						*partnr = l;
						*coltyp = 1;
					}
				}
				else				//DISSOCIATE OR BOUNCE, miss core collision, particles collide at sigma2
				{
					tij = (-bij + sqrt(disc2))/vijsq;
					if(tij<*coltim)
					{
						*coltim = tij;
						*partnr = l;
						*coltyp = 2;
					}
				}
			}
			else 	//particles outside attractive range, moving towards each other
			{
				if(disc2 > 0.0 && chk == 1)			//CAPTURE
				{
					tij = (-bij - sqrt(disc2))/vijsq;
					if(tij<*coltim)
					{
						*coltim = tij;
						*partnr = l;
						*coltyp = 3;
					}
				}
			}
        }
		else if(bij > 0.0 && cij2 < 0.0)		//particles receding
		{									//particles within attractive range: BOUNCE/DISSOCIATE
			tij = (-bij + sqrt(disc2))/vijsq;
			if(tij<*coltim)
			{
				*coltim = tij;
				*partnr = l;
				*coltyp = 2;
			}
		}
	}
}

void dnlist(int j, double sigma1, double sigma2, double d1, double d2, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], double coltim[], int partnr[], int coltyp[], int **linkedspheres)
{
    double rxj, ryj, rzj, vxj, vyj, vzj, sig11, sig22, sig1sq, sig2sq;
    double rxij, ryij, rzij, vxij, vyij, vzij;
    double bij, rijsq, vijsq, disc1, disc2, cij1, cij2, tij;
    int l,linked,chk,loopkvar;

    if(j==0)        //if particle number 0 was involved in the collision, dnlist is not needed
        return;
    rxj = rx[j];
    ryj = ry[j];
    rzj = rz[j];
    vxj = vx[j];
    vyj = vy[j];
    vzj = vz[j];

    for(l=0;l<j;l++)
    {
        rxij = rx[l] - rxj;
        ryij = ry[l] - ryj;
        rzij = rz[l] - rzj;
    	rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0)+(0.5*copysign(1.0,rxij)))));
       	ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0)+(0.5*copysign(1.0,ryij)))));
        rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0)+(0.5*copysign(1.0,rzij)))));
        vxij = vx[l] - vxj;
        vyij = vy[l] - vyj;
        vzij = vz[l] - vzj;
        bij = rxij*vxij + ryij*vyij + rzij*vzij;
        rijsq = pow(rxij, 2.0) + pow(ryij, 2.0) + pow(rzij, 2.0);
        vijsq = pow(vxij, 2.0) + pow(vyij, 2.0) + pow(vzij, 2.0);
        linked = 0;
		loopkvar = 0;
		while(loopkvar < sizeof(linkedspheres) && linkedspheres[l][loopkvar] != -1)
		{
			if(linkedspheres[l][loopkvar] == j)
			{
				linked = 1;
				break;
			}
			loopkvar += 1;
		}
		if( linked == 1 )
		{
			sig11 = d1;
			sig22 = d2;
			chk = 0;
		}
		else
		{
			sig11 = sigma1;
			sig22 = sigma2;
			chk = 1;
		}
		sig1sq = sig11 * sig11;
		sig2sq = sig22 * sig22;
		cij1 = rijsq - sig1sq;
		cij2 = rijsq - sig2sq;
		disc1 = pow(bij,2.0) - vijsq * cij1;
		disc2 = pow(bij,2.0) - vijsq * cij2;
		if(bij < 0.0)       //particles approaching
        {
            if(cij2 < 0.0)		//particles within attractive range, moving towards each other
			{

				if(disc1 > 0.0)		//CORE collision
				{
					tij = (-bij - sqrt(disc1)) / vijsq;
					if(tij < coltim[l])
					{
						coltim[l] = tij;
						partnr[l] = j;
						coltyp[l] = 1;
					}
				}
				else				//DISSOCIATE OR BOUNCE, miss core collision, particles collide at sigma2
				{
					tij = (-bij + sqrt(disc2))/vijsq;
					if(tij < coltim[l])
					{
						coltim[l] = tij;
						partnr[l] = j;
						coltyp[l] = 2;
					}
				}
			}
			else			//particles outside attractive range, moving towards each other
			{
				if(disc2 > 0.0 && chk == 1)			//CAPTURE
				{
					tij = (-bij - sqrt(disc2))/vijsq;
					if(tij < coltim[l])
					{
						coltim[l] = tij;
						partnr[l] = j;
						coltyp[l] = 3;
					}
				}
			}
        }
		else if(bij > 0.0 && cij2 < 0.0)		//particles receding
		{
			tij = (-bij + sqrt(disc2))/vijsq;
			if(tij < coltim[l])
			{
				coltim[l] = tij;
				partnr[l] = j;
				coltyp[l] = 2;
			}
		}
    }
}
      

double bump(double sigma1, double sigma2, double d1, double d2, int i, int j, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], int *coltyp, int linkedspheres[])
{
    double sig11, sig22, sig1sq, sig2sq, rxij, ryij, rzij, vxij, vyij, vzij, bij, bijsq, rijsq, vijsq, disc1, disc2, delpe, cij1, cij2, factor, delvx, delvy, delvz, w, smbump, smdist;
	int lengthofarray, linked, l, k, chk;
	lengthofarray = sizeof(linkedspheres[0]);
	l = *coltyp;
	linked = 0;
	for(k=0;k<lengthofarray;k++)
		if(linkedspheres[k] == j)
			{
				linked = 1;
				break;
			}
	if( linked == 1 ) 
		{
			sig11 = sigma1*(1-del);
			sig22 = sigma1*(1+del);
			chk = 0;
		}
		else
		{
			sig11 = sigma1;
			sig22 = sigma2;
			chk = 1;
		}
	
	sig1sq = pow(sig11, 2.0);
    sig2sq = pow(sig22, 2.0);
    rxij = rx[i] - rx[j];
    ryij = ry[i] - ry[j];
    rzij = rz[i] - rz[j];
    rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0)+(0.5*copysign(1.0,rxij)))));
    ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0)+(0.5*copysign(1.0,ryij)))));
    rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0)+(0.5*copysign(1.0,rzij)))));
	vxij = vx[i] - vx[j];
	vyij = vy[i] - vy[j];
	vzij = vz[i] - vz[j];
	bij = rxij*(vxij) + ryij*(vyij) + rzij*(vzij);
	bijsq = pow(bij, 2.0);
	rijsq = pow(rxij, 2.0) + pow(ryij, 2.0) + pow(rzij, 2.0);
    vijsq = pow(vxij, 2.0) + pow(vyij, 2.0) + pow(vzij, 2.0);
	cij1 = rijsq - sig1sq;
	cij2 = rijsq - sig2sq;
	delpe = 4 * sig2sq * eps;
	disc1 = bijsq - vijsq * cij1;
	disc2 = bijsq - delpe;
	factor = 0.0;
	smbump = 0.0;
	smdist = 5e-12;
	if(l == 1)
	{
			factor = bij/sig1sq;
			smbump = 0.0;
	}
	else if(l == 2)
	{	//bounce or dissociate
		if(chk == 0)						//bounce
		{
			factor = bij/sig2sq;
			smbump = -1.0;
		}
		else if(chk == 1 && disc2 < 0)		//bounce
		{
			factor = bij/sig2sq;
			smbump = -1.0;
		}
		else 								//dissocaite
		{
			factor = 0.5*(bij-sqrt(bijsq-delpe))/sig2sq;
			smbump = 1.0;
		}
	}
	else if(l == 3)							//capture
	{
		factor = 0.5*(sqrt(delpe+bijsq)+bij)/sig2sq;
		smbump = -1.0;
	}
	else
	{
	printf("Switch Case fail");
	exit(1);
	}
  
	delvx = -factor*rxij;       //change in velocity
    delvy = -factor*ryij;
    delvz = -factor*rzij;

    vx[i] = vx[i] + delvx;      //calculate new velocity for i and j
    vx[j] = vx[j] - delvx;

    vy[i] = vy[i] + delvy;
	vy[j] = vy[j] - delvy;

    vz[i] = vz[i] + delvz;
    vz[j] = vz[j] - delvz;
	
    w = delvx*rxij + delvy*ryij + delvz*rzij;

	rx[i] = rx[i] + smbump*smdist*rxij*sig22;
	ry[i] = ry[i] + smbump*smdist*ryij*sig22;
	rz[i] = rz[i] + smbump*smdist*rzij*sig22;
	rx[j] = rx[j] - smbump*smdist*rxij*sig22;
	ry[j] = ry[j] - smbump*smdist*ryij*sig22;
	rz[j] = rz[j] - smbump*smdist*rzij*sig22;
	
    return(w);
}

void grsort(double rx[], double ry[], double rz[], int hist[])
{
    double rxi, ryi, rzi, rxij, ryij, rzij, rijsq, rij;
    int i, j, h, bin, chnid[n];
	
	for(i=0;i<n;i++)
		chnid[i] = (int)((i/chnlen)+1);

    for(i=0;i<(n-1);i++)
    {
        rxi = rx[i];
        ryi = ry[i];
        rzi = rz[i];
        for(j=(i+1);j<n;j++)
        {
			if(chnid[i] != chnid[j])
			{
			rxij = rxi - rx[j];
            ryij = ryi - ry[j];
            rzij = rzi - rz[j];
			rxij = rxij - (1.0*(double)((int)(rxij/(double)(1.0) + (0.5*copysign(1.0,rxij)))));
			ryij = ryij - (1.0*(double)((int)(ryij/(double)(1.0) + (0.5*copysign(1.0,ryij)))));
			rzij = rzij - (1.0*(double)((int)(rzij/(double)(1.0) + (0.5*copysign(1.0,rzij)))));
            rijsq = pow(rxij, 2) + pow(ryij, 2) + pow(rzij, 2);
            rij = sqrt(rijsq);
            bin = (int) (rij/delr);
            if(bin<maxbin)
                hist[bin] = hist[bin] + 2;
			}
		}
    }
}

void ghostcoll(double temperature, double *vx, double *vy, double *vz)
{
	double rtemp;
	rtemp = sqrt(temperature);
	*vx = rtemp * gauss();
	*vy = rtemp * gauss();
	*vz = rtemp * gauss();
}

//Generate Gaussian distribution for velocities Allen and Tildesley pg. 347
double gauss(void)
{
    double a1, a3, a5, a7, a9, rsum, r, r2, g;
    int a;

    a1 = 3.949846138;
    a3 = 0.252408784;
    a5 = 0.076542912;
    a7 = 0.008355968;
    a9 = 0.029899776;

    rsum = 0.0;
    for(a=0 ; a<12 ; a++)
    {    
		rsum = rsum + ranf();
	}
	
    r = (rsum - 6.0) / 4.0;
    r2 = r * r;
    g = (((( a9 * r2 + a7) * r2 + a5) * r2 + a3) * r2 + a1) * r;
    return(g);
}

double ranf()
{
	double dra;
	dra = rand()/(double)RAND_MAX;
	return dra;
}