/*
 * hsphere.h
 *
 *  Created on: May 8, 2012
 *      Author: David
 */

#ifndef HSPHERE_H_
#define HSPHERE_H_

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

const int n = 500;
const double timbig = 1.0E10;
const double tol = 1.0E-4;
const double deltawrite = 0.01;

const double pi = 3.14159265359;
const int maxbin = 7000;
const double delr = 0.0001;

void fcc(double rx[], double ry[], double rz[]);
void comvel(double temp, double rx[], double ry[], double rz[]);
void gauss(double dummy);
double ranf();
double gauss();

unsigned char check(double sigma, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[]);
double energy(double vx[], double vy[], double vz[]);

void grsort(double rx[], double ry[], double rz[], int hist[]);
void uplist(double sigma, int i, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], double *coltime, int *partner);
double bump(double sigma, int nextsphere, int temppartner, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[]);
void dnlist(double sigma, int j, double rx[], double ry[], double rz[], double vx[], double vy[], double vz[], double coltime[], int partner[]);
#endif /* HSPHERE_H_ */
