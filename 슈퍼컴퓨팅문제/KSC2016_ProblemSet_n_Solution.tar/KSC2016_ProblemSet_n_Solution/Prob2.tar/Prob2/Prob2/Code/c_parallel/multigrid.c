#include <stdio.h>
#include "multigrid.h"

#define URF(I,J) uf[(I)*(2*ncy)+(J)]
#define URC(I,J) uc[(I)*(ncy)+(J)]
#define UIC(I,J) uc_ghost[(I)*(ncy+2)+(J)]
#define UIF(I,J) uf[(I)*(nfy)+(J)]
#define UIFD(I,J) uf_double[(I)*(nfy)+(J)]

void restriction(double *uf, double *uc, int ncx, int ncy)
{
	int ic, jc;

	for (ic=0; ic<ncx; ic++) 
	{
		for (jc=0; jc<ncy; jc++) 
		{
			URC(ic,jc) = 1.0/4.0*( \
				URF(2*ic,2*jc)+URF(2*ic+1,2*jc) \
				+URF(2*ic,2*jc+1)+URF(2*ic+1,2*jc+1));
		}
	}
}

void restriction_mlevel(double *uf, double *uc, int ncx, int ncy, int my_rank, int distance)
{
	int ic=0, jc;
	double *uf_next;
	MPI_Status status;

	uf_next = (double*)malloc((ncx)*(ncy*2)*sizeof(double));

	if(my_rank%(2*distance)==0) 
	{
		MPI_Recv(uf_next, ncy*2, MPI_DOUBLE, my_rank+distance, distance, MPI_COMM_WORLD, &status);
		for (jc=0; jc<ncy; jc++) 
		{
			URC(ic,jc) = 1.0/4.0*( uf[2*jc]+uf[2*jc+1]+uf_next[2*jc]+uf_next[2*jc+1]);
		}
	}
	else if(my_rank%(2*distance)==distance)
	{
		MPI_Send(uf, ncy*2, MPI_DOUBLE, my_rank-distance, distance, MPI_COMM_WORLD);
	}

	free(uf_next);
}

void interpolation_mlevel(double *uc, double *uf, int nfx, int nfy, int my_rank, int mpisize, int distance)
{
	int ic,jc,ncx,ncy,count=0;
	double *uc_ghost;
	double *uf_double;
	int i,j;
	MPI_Status status;

	ncx = nfx;
	ncy = nfy>>1;
	uc_ghost = (double*)malloc((ncx+2)*(ncy+2)*sizeof(double));
	uf_double = (double*)malloc((2*nfx)*(nfy)*sizeof(double));

	if(my_rank%(2*distance)==0) 
	{
		for (i=1; i<ncx+1; i++) 
		{
			for (j=1; j<ncy+1; j++) 
			{
				UIC(i,j) = uc[count];
				count++;
			}
		}
		for (i=1; i<ncx+1; i++) 
		{
			UIC(i,0) = -UIC(i,1);
			UIC(i,ncy+1) = -UIC(i,ncy);
		}

		if(my_rank!=0)
		{
			MPI_Recv(&UIC(0,0),ncy+2, MPI_DOUBLE, my_rank-2*distance, 0, MPI_COMM_WORLD, &status);
			MPI_Send(&UIC(1,0),ncy+2, MPI_DOUBLE, my_rank-2*distance, 0, MPI_COMM_WORLD);
		}
		if(my_rank+2*distance!=mpisize)
		{	
			MPI_Send(&UIC(ncx  ,0),ncy+2, MPI_DOUBLE, my_rank+2*distance, 0, MPI_COMM_WORLD);
			MPI_Recv(&UIC(ncx+1,0),ncy+2, MPI_DOUBLE, my_rank+2*distance, 0, MPI_COMM_WORLD, &status);
		}	

		if(my_rank==0)
		{
			for (j=1; j<ncy+1; j++) 
			{
				UIC(0,j) = -UIC(1,j);
			}
			UIC(0,0) = -UIC(1,1);
			UIC(0,ncy+1) = -UIC(1,ncy);
		}
		if(my_rank+2*distance==mpisize)
		{
			for (j=1; j<ncy+1; j++) 
			{
				UIC(ncx+1,j) = -UIC(ncx,j);
			}
			UIC(ncx+1,0) = -UIC(ncx,1);
			UIC(ncx+1,ncy+1) = -UIC(ncx,ncy);
		}

	
		for (ic=0; ic<ncx; ic++) 
		{
			for (jc=0; jc<ncy; jc++) 
			{
				UIFD(2*ic,2*jc) = (9.0*UIC(ic+1,jc+1) \
						+ 3.0*UIC(ic,jc+1) \
						+ 3.0*UIC(ic+1,jc) \
						+ 1.0*UIC(ic,jc))/16.0;
	
				UIFD(2*ic+1,2*jc) = (9.0*UIC(ic+1,jc+1) \
						+ 3.0*UIC(ic+2,jc+1) \
						+ 3.0*UIC(ic+1,jc) \
						+ 1.0*UIC(ic+2,jc))/16.0;
	
				UIFD(2*ic,2*jc+1) = (9.0*UIC(ic+1,jc+1) \
						+ 3.0*UIC(ic,jc+1) \
						+ 3.0*UIC(ic+1,jc+2) \
						+ 1.0*UIC(ic,jc+2))/16.0;
	
				UIFD(2*ic+1,2*jc+1) = (9.0*UIC(ic+1,jc+1) \
						+ 3.0*UIC(ic+2,jc+1) \
						+ 3.0*UIC(ic+1,jc+2) \
						+ 1.0*UIC(ic+2,jc+2))/16.0;
			}
		}
	
		count=0;
		for (i=0; i<nfx; i++) 
		{
			for (j=0; j<nfy; j++) 
			{
				uf[count]=UIFD(i,j);
				count++;
	
			}
		}
		MPI_Send(&UIFD(1,0), nfy, MPI_DOUBLE, my_rank+distance, distance, MPI_COMM_WORLD);
	}
	else if(my_rank%(2*distance)==distance)
	{
		MPI_Recv(uf, nfy, MPI_DOUBLE, my_rank-distance, distance, MPI_COMM_WORLD, &status);
	}
	free(uf_double);
	free(uc_ghost);
}

void interpolation_mpi(double *uc, double *uf, int nfx, int nfy, int my_rank, int prev_rank, int next_rank)
{
	int ic,jc,ncx,ncy,count=0;
	double *uc_ghost;
	int i,j;
	MPI_Status status;
	ncx = nfx>>1;
	ncy = nfy>>1;
	uc_ghost = (double*)malloc((ncx+2)*(ncy+2)*sizeof(double));

	for (i=1; i<ncx+1; i++) 
	{
		for (j=1; j<ncy+1; j++) 
		{
			UIC(i,j) = uc[count];
			count++;
		}
	}
	for (i=1; i<ncx+1; i++) 
	{
		UIC(i,0) = -UIC(i,1);
		UIC(i,ncy+1) = -UIC(i,ncy);
	}

	if(prev_rank!=MPI_PROC_NULL)
	{
		MPI_Recv(&UIC(0,0),ncy+2, MPI_DOUBLE, prev_rank, 0, MPI_COMM_WORLD, &status);
		MPI_Send(&UIC(1,0),ncy+2, MPI_DOUBLE, prev_rank, 0, MPI_COMM_WORLD);
	}
	if(next_rank!=MPI_PROC_NULL)
	{
		MPI_Send(&UIC(ncx  ,0),ncy+2, MPI_DOUBLE, next_rank, 0, MPI_COMM_WORLD);
		MPI_Recv(&UIC(ncx+1,0),ncy+2, MPI_DOUBLE, next_rank, 0, MPI_COMM_WORLD, &status);
	}

	if(prev_rank==MPI_PROC_NULL) 
	{
		for (j=1; j<ncy+1; j++) 
		{
			UIC(0,j) = -UIC(1,j);
		}
		UIC(0,0) = -UIC(1,1);
		UIC(0,ncy+1) = -UIC(1,ncy);
	}
	if(next_rank==MPI_PROC_NULL) 
	{
		for (j=1; j<ncy+1; j++) 
		{
			UIC(ncx+1,j) = -UIC(ncx,j);
		}
		UIC(ncx+1,0) = -UIC(ncx,1);
		UIC(ncx+1,ncy+1) = -UIC(ncx,ncy);
	}


	for (ic=0; ic<ncx; ic++) 
	{
		for (jc=0; jc<ncy; jc++) 
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

	free(uc_ghost);
}
