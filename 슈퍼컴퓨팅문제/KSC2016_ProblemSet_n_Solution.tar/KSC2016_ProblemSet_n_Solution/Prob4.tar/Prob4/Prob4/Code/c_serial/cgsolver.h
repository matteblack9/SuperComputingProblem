#ifndef CGSOLVER_H_
#define CGSOLVER_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void cgsolver(int size, double *matrix, double *rhs, double *solution, int maxiteration, double tolerance);
double innerproduct(double *x, double *y, int size);
void multiply(int size, double *matrix, double *x, double *y);

#endif
