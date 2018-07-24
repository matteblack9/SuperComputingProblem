#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "multigrid.h"
#include <mpi.h>

#define PI 3.141592653589793

int main(int argc, char **argv) 
{
	int nlevel = 15;
	double length_x = 1.0, length_y = 1.0;
	int maxiteration = 20;

	int ngrid;
	int mlevel,gsize; 
	int numgrid_x[nlevel], numgrid_y[nlevel];
	int matrixDOF[nlevel];  
	int count, ii, jj, kk, i, iter;  
	double gridsize[nlevel];
	double *rhs[nlevel], *solution[nlevel];
    double coord_x, coord_y;
	double time_start, time_end;
	double time_elapsed;

	int mpisize,my_rank,rank;
	int prev_rank, next_rank, distance;

	char myfilename[256];
	FILE *fp;
	MPI_Status status;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&mpisize);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

	time_start = MPI_Wtime();

	ngrid = pow(2,nlevel-1);

	if(my_rank == 0)
	{
		prev_rank=MPI_PROC_NULL;
	}
	else
	{
		prev_rank = my_rank - 1;
	}
	if(my_rank == mpisize-1)
	{
		next_rank=MPI_PROC_NULL;
	}
	else
	{
		next_rank = my_rank + 1;
	}

	numgrid_x[nlevel-1] = ngrid/mpisize;
	numgrid_y[nlevel-1] = ngrid;
	gridsize[nlevel-1] = length_x/(double)ngrid;
	matrixDOF[nlevel-1] = numgrid_x[nlevel-1]*numgrid_y[nlevel-1];

	gsize=mpisize;
	mlevel=0;
	for(i=0;i<nlevel;i++)
	{
		gsize/=2;
		mlevel++;
		if(gsize==1) break;
	}

	for(i=nlevel-1; i>mlevel; i--)
	{
		numgrid_x[i-1] = numgrid_x[i]/2;
		numgrid_y[i-1] = numgrid_y[i]/2;
		gridsize[i-1] = length_y/numgrid_y[i-1];
		matrixDOF[i-1] = numgrid_x[i-1]*numgrid_y[i-1];
	}
	for(i=mlevel; i>0; i--)
	{
		numgrid_x[i-1] = 1;
		numgrid_y[i-1] = numgrid_y[i]/2;
		gridsize[i-1] = length_y/numgrid_y[i-1];
		matrixDOF[i-1] = numgrid_x[i-1]*numgrid_y[i-1];
	}

	if(!my_rank) printf("[Multigrid] Geometry and matrix size initialized.\n");

	for(i=0; i<nlevel; i++) 
	{
		rhs[i] = (double*)malloc(matrixDOF[i]*sizeof(double));
		solution[i] = (double*)malloc(matrixDOF[i]*sizeof(double));
	}

	count = 0;

	for(ii=0; ii<numgrid_x[nlevel-1]; ii++)
	{
		coord_x = length_x/(double)mpisize*(double)my_rank+ii*gridsize[nlevel-1]+0.5*gridsize[nlevel-1];
		for(jj=0; jj<numgrid_y[nlevel-1]; jj++)
		{
			coord_y = jj*gridsize[nlevel-1]+0.5*gridsize[nlevel-1];
			rhs[nlevel-1][count] = sin(coord_x/length_x*PI)* \
				   sin(coord_y/length_y*PI);
			count++;
		} 
	}	

	if(!my_rank) printf("[Multigrid] Geometry and rhs constructed.\n");
	if(!my_rank) printf("[Multigrid] Start solving equations.\n");

	for(iter=0;iter<maxiteration;iter++)
	{
		if(!my_rank) printf("[Multigrid] iteration %d \n",iter+1);
		for(i=nlevel-1; i>mlevel; i--) 
		{
			restriction(rhs[i],rhs[i-1],numgrid_x[i-1],numgrid_y[i-1]);
			if(!my_rank) printf("[Multigrid] Restriction of RHS done (level %d -> level %d).\n", i, i-1);

		}
		distance = 1;
		for(i=mlevel; i>0; i--) 
		{
			restriction_mlevel(rhs[i],rhs[i-1],numgrid_x[i-1],numgrid_y[i-1],my_rank,distance);
			distance = distance * 2;
			if(!my_rank) printf("[Multigrid] Restriction of RHS done (level %d -> level %d).\n", i, i-1);
		}
	
		if(my_rank == 0)
		{
			solution[0][0] = -rhs[0][0]/4.0;
			printf("[Multigrid] Solution at the coarsest level = %f\n",solution[0][0]);
		}
	
		for(i=1; i<=mlevel; i++)
		{
			distance = distance / 2;
			interpolation_mlevel(solution[i-1],solution[i],numgrid_x[i],numgrid_y[i],my_rank,mpisize,distance);
			if(!my_rank) printf("[Multigrid] Interpolation of solution in mlevel done (level %d -> level %d). \n", i-1, i);
		}
	
		for(i=mlevel+1; i<nlevel; i++)
		{
			interpolation_mpi(solution[i-1],solution[i],numgrid_x[i],numgrid_y[i],my_rank,prev_rank,next_rank);
			if(!my_rank) printf("[Multigrid] Interpolation of solution in nlevel done (level %d -> level %d). \n", i-1, i);
		}
		for(ii=0;ii<matrixDOF[nlevel-1];ii++)
		{
			rhs[nlevel-1][ii] = PI*PI*PI*PI*solution[nlevel-1][ii];
		}
	}	
	if(!my_rank) printf("[Multigrid] Final solution obtained.\n");

	MPI_Barrier(MPI_COMM_WORLD);

	i=nlevel-6; 

	sprintf(myfilename,"./solution.dat",i);
	if(my_rank==(mpisize/2-1))
	{
		fp = fopen(myfilename, "a");
		for(ii=-8; ii<0; ii++)
		{
			coord_x = length_x/(double)mpisize*(double)my_rank+(ii+numgrid_x[i])*gridsize[i]+0.5*gridsize[i];
			for(jj = -8; jj < 8; jj++)
			{
				coord_y = (jj+numgrid_y[i]/2)*gridsize[i]+0.5*gridsize[i];

				count = (ii+numgrid_x[i])*numgrid_y[i] + jj+numgrid_y[i]/2;
				fprintf(fp,"%12.6f %12.6f %12.6f\n", coord_x, coord_y, solution[i][count]);
			}
		}
		fclose(fp);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank==(mpisize/2))
	{
		fp = fopen(myfilename, "a");
		for(ii=0; ii<8; ii++)
		{
			coord_x = length_x/(double)mpisize*(double)my_rank+(ii)*gridsize[i]+0.5*gridsize[i];
			for(jj = -8; jj < 8; jj++)
			{
				coord_y = (jj+numgrid_y[i]/2)*gridsize[i]+0.5*gridsize[i];
				count = (ii)*numgrid_y[i] + jj+numgrid_y[i]/2;
				fprintf(fp,"%12.6f %12.6f %12.6f\n", coord_x, coord_y, solution[i][count]);
			}
		}
		fclose(fp);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(!my_rank) printf("[Multigrid] Final solution printed.\n");

	for(i = 0; i < nlevel; i++) 
	{
	        free(solution[i]);
		free(rhs[i]);
	}

	if(!my_rank) printf("[Multigrid] Memory deallocated.\n");

	time_end = MPI_Wtime();
	time_elapsed = time_end - time_start;
	if(!my_rank) printf("[Multigrid] Mission completed in %f (secs).\n", time_elapsed);


	MPI_Finalize();
	return 0;
}


