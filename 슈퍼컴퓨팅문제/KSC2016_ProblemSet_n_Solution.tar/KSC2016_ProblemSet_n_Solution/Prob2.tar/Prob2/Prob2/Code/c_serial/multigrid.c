#include <stdio.h>
#include "multigrid.h"

#define URF(I,J) uf[(I)*(2*nc)+(J)]
#define URC(I,J) uc[(I)*(nc)+(J)]
#define UIC(I,J) ucwithghost[(I)*(nc+2)+(J)]
#define UIF(I,J) uf[(I)*(nf)+(J)]

void restriction(double *uf, double *uc, int nc)
{
    int ic, jc;

    for (ic=0; ic<nc; ic++) 
    {
	for (jc=0; jc<nc; jc++) 
	{
	    URC(ic,jc) = 1.0/4.0*( \
		    URF(2*ic,2*jc)+URF(2*ic+1,2*jc) \
		    +URF(2*ic,2*jc+1)+URF(2*ic+1,2*jc+1));
	}
    }
}

void interpolation(double *uc, double *uf, int nf)
{
    int ic,jc,nc,count=0;
    double *ucwithghost;
    int i,j;
    nc = nf>>1;
    ucwithghost = (double*)malloc((nc+2)*(nc+2)*sizeof(double));

    for (i=1; i<nc+1; i++) 
    {
	for (j=1; j<nc+1; j++) 
	{
	    UIC(i,j) = uc[count];
	    count++;
	}
    }
    for (i=1; i<nc+1; i++) 
    {
	UIC(i,0) = -UIC(i,1);
	UIC(i,nc+1) = -UIC(i,nc);
    }

    for (j=1; j<nc+1; j++) 
    {
	UIC(0,j) = -UIC(1,j);
	UIC(nc+1,j) = -UIC(nc,j);
    }

    UIC(0,0) = -UIC(1,1);
    UIC(0,nc+1) = -UIC(1,nc);
    UIC(nc+1,0) = -UIC(nc,1);
    UIC(nc+1,nc+1) = -UIC(nc,nc);

    for (ic=0; ic<nc; ic++) 
    {
	for (jc=0; jc<nc; jc++) 
	{
	    UIF(2*ic,2*jc) = (9.0*UIC(ic+1,jc+1) \
		    + 3.0*UIC(ic,jc+1) \
		    + 3.0*UIC(ic+1,jc) \
		    + 1.0*UIC(ic,jc))/16.0;

	    UIF(2*ic+1,2*jc) = (9.0*UIC(ic+1,jc+1) \
		    + 3.0*UIC(ic+2,jc+1) \
		    + 3.0*UIC(ic+1,jc) \
		    + 1.0*UIC(ic+2,jc))/16.0;

	    UIF(2*ic,2*jc+1) = (9.0*UIC(ic+1,jc+1) \
		    + 3.0*UIC(ic,jc+1) \
		    + 3.0*UIC(ic+1,jc+2) \
		    + 1.0*UIC(ic,jc+2))/16.0;

	    UIF(2*ic+1,2*jc+1) = (9.0*UIC(ic+1,jc+1) \
		    + 3.0*UIC(ic+2,jc+1) \
		    + 3.0*UIC(ic+1,jc+2) \
		    + 1.0*UIC(ic+2,jc+2))/16.0;
	}
    }

    free(ucwithghost);
}

