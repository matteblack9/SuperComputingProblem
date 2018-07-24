#ifndef _MULTIGRID_H
#define _MULTIGRID_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void restriction(double *uf, double *uc, int nx, int ny);
void interpolation(double *uc, double *uf, int nx, int ny, int myrank, int ncpus);

#endif
