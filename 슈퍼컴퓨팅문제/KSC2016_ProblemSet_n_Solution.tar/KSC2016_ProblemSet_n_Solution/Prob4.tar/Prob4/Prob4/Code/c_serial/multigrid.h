#ifndef _MULTIGRID_H
#define _MULTIGRID_H

#include <stdio.h>
#include <stdlib.h>

void restriction(double *uf, double *uc, int nc);
void interpolation(double *uc, double *uf, int nf);

#endif
