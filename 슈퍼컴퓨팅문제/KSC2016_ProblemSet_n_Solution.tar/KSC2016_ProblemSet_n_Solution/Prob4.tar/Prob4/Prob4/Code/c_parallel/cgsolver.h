#ifndef CGSOLVER_H_
#define CGSOLVER_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void cgsolver(int size, int myrank, int ncpus, double *matrix, double *rhs, double *solution, int maxiteration, double tolerance);
double innerproduct(double *x, double *y, int size);
void multiply(int *rank, int *load, int *sindex, int ncpus, double *matrix, double *x, double *y, double *xtl, double *xtr);

#endif
