/*
 * swell.h
 *
 *  Created on: May 15, 2012
 *      Author: David
 */

#ifndef SWELL_H_
#define SWELL_H_

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

const int n = 500;
const double timbig = 1.0E10;
const double tol = 1.0E-4;

const double pi = 3.14159265359;
const int maxbin = 7000;
const double delr = 0.0001;

const double ghostcoeff = 0.001;

void fcc(double rx[], double ry[], double rz[]);
void comvel(double temp, double rx[], double ry[], double rz[]);
void gauss(double dummy);
double ranf();
double gauss();

unsigned char check(double sigma, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[]);
double kineticenergy(double vx[], double vy[], double vz[]);
double potentialenergy(double rx[], double ry[], double rz[], double sigma2, double energy);

void grsort(double rx[], double ry[], double rz[], int hist[]);
void uplist(double sigma, double sigma2, int i, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], double *coltime, int *partner, int *colltype);
double bump(double sigma, double sigma2, double energy, int nextsphere, int temppartner, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], int colltype[]);
void dnlist(double sigma, double sigma2, int j, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], double coltime[], int partner[], int *colltype);
void ghostcoll(double usrtemp, double *vx, double *vy, double *vz);

#endif /* SWELL_H_ */
