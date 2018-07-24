#ifndef _MULTIGRID_H
#define _MULTIGRID_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void restriction(double *uf, double *uc, int ncx, int ncy);
void interpolation_mpi(double *uc, double *uf, int nfx, int nfy, int my_rank, int prev_rank, int next_rank);
void restriction_mlevel(double *uf, double *uc, int ncx, int ncy, int my_rank, int distance);
void interpolation_mlevel(double *uc, double *uf, int nfx, int nfy, int my_rank, int mpisize, int distance);

#endif
