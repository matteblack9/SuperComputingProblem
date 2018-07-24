#include <stdio.h>
#include "multigrid.h"

#define URF(I,J) uf[(I)*(2*ny)+(J)]
#define URC(I,J) uc[(I)*(ny)+(J)]
#define UIC(I,J) ucwithghost[(I)*(nc_y+2)+(J)]
#define UIF(I,J) uf[(I)*(ny)+(J)]

void restriction(double *uf, double *uc, int nx, int ny)
{
    int ic, jc;

    for (ic=0; ic<nx; ic++) 
    {
	for (jc=0; jc<ny; jc++) 
	{
	    URC(ic,jc) = 1.0/4.0*( \
		    URF(2*ic,2*jc)+URF(2*ic+1,2*jc) \
		    +URF(2*ic,2*jc+1)+URF(2*ic+1,2*jc+1));
	}
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void interpolation(double *uc, double *uf, int nx, int ny, int myrank, int ncpus)
{
    int ic, jc, nc_x, nc_y, count=0;
    int leftrank, rightrank;
    double *ucwithghost;
    int i,j;
    MPI_Request req[2];
    MPI_Status status[2];

    nc_x = nx/2;
    nc_y = ny/2;
   
    leftrank = (myrank-1+ncpus)%ncpus;
    rightrank = (myrank+1)%ncpus;

    ucwithghost = (double*)malloc((nc_x+2)*(nc_y+2)*sizeof(double));

    for (i=1; i<nc_x+1; i++) 
    {
	for (j=1; j<nc_y+1; j++) 
	{
	    UIC(i,j) = uc[count];
	    count++;
	}
    }
    
    for (i=1; i<nc_x+1; i++) 
    {
	UIC(i,0) = -UIC(i,1);
	UIC(i,nc_y+1) = -UIC(i,nc_y);
    }

    if(ncpus>1) 
    {
	MPI_Irecv(&UIC(nc_x+1,0), nc_y+2, MPI_DOUBLE, rightrank, 1002, MPI_COMM_WORLD, &req[0]);
	MPI_Isend(&UIC(1,0), nc_y+2, MPI_DOUBLE, leftrank, 1002, MPI_COMM_WORLD, &req[1]);
	MPI_Waitall(2, req, status);
	
        MPI_Irecv(&UIC(0,0), nc_y+2, MPI_DOUBLE, leftrank, 1002, MPI_COMM_WORLD, &req[0]);
        MPI_Isend(&UIC(nc_x,0), nc_y+2, MPI_DOUBLE, rightrank, 1002, MPI_COMM_WORLD, &req[1]);
        MPI_Waitall(2, req, status);	
    }

    for (j=1; j<nc_y+1; j++) 
    {
	if(myrank==0)
		UIC(0,j) = -UIC(1,j);
	if(myrank==ncpus-1)
		UIC(nc_x+1,j) = -UIC(nc_x,j);
    }

    if(myrank==0)
    {
	UIC(0,0) = -UIC(1,1);
    	UIC(0,nc_y+1) = -UIC(1,nc_y);
    }

    if(myrank==ncpus-1)
    {
    	UIC(nc_x+1,0) = -UIC(nc_x,1);
    	UIC(nc_x+1,nc_y+1) = -UIC(nc_x,nc_y);
    }

    for (ic=0; ic<nc_x; ic++) 
    {
	for (jc=0; jc<nc_y; jc++) 
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

    MPI_Barrier(MPI_COMM_WORLD);
}

