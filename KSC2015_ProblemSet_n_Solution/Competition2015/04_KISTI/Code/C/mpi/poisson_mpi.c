#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "matrixconstructor.h"
#include "cgsolver_mpi.h"

#define PI 3.141592653589793

int main(int argc, char **argv) {

	int numgrid_x, numgrid_y, firstgrid_x, lastgrid_x;
	int firstrow, lastrow, matrixDOF, globalmatrixDOF;  
	int count, ii, jj, kk, maxiteration, myrank, ncpus;  
	double length_x, length_y, gridsize, tolerance;
	double *poissonmatrix, *rhs, *solution, *coordinate;
	double elapsed_time;
	char myfilename[256];
	FILE *fp;

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        MPI_Comm_size(MPI_COMM_WORLD, &ncpus);

	if(myrank == 0)
		printf("[Poisson] MPI initialized.\n");

	elapsed_time = -MPI_Wtime();

	maxiteration = 20000;
	tolerance = 1.0e-10;
	gridsize = 0.01;
	length_x = 1.0;
	length_y = 1.0;

	numgrid_x = (int)round(length_x/gridsize) + 1;
	numgrid_y = (int)round(length_y/gridsize) + 1;
	globalmatrixDOF = numgrid_x*numgrid_y;

        firstgrid_x = (numgrid_x*myrank/ncpus) + 1;
	lastgrid_x = numgrid_x*(myrank+1)/ncpus;
        firstrow = (firstgrid_x-1)*numgrid_y + 1;
        lastrow = lastgrid_x*numgrid_y;
	matrixDOF = lastrow - firstrow + 1;

	if(myrank == 0)
		printf("[Poisson] Geometry and matrix size initialized.\n");

	coordinate = (double*)malloc(matrixDOF*2*sizeof(double));
	rhs = (double*)malloc(matrixDOF*sizeof(double));

	count = 0;

	for(ii=firstgrid_x; ii<=lastgrid_x; ii++)
	{
		for(jj=1; jj<=numgrid_y; jj++)
		{
			coordinate[2*count] = (ii-1)*gridsize;
			coordinate[2*count+1] = (jj-1)*gridsize;
			rhs[count] = sin(coordinate[2*count]/length_x*PI)*sin(coordinate[2*count+1]/length_y*PI)*gridsize*gridsize;
			count++;
		} 
	}

	if(myrank == 0)	
		printf("[Poisson] Geometry and rhs constructed.\n");

	poissonmatrix = (double*) malloc(matrixDOF*globalmatrixDOF*sizeof(double)); 

	for(ii=0; ii<matrixDOF*globalmatrixDOF; ii++)	
		poissonmatrix[ii] = 0.0;

	construct_poissonmatrix(firstrow, firstgrid_x, lastgrid_x, numgrid_x, numgrid_y, poissonmatrix);

	if(myrank == 0)
		printf("[Poisson] Poisson matrix constructed.\n");

	solution = (double*) malloc(matrixDOF*sizeof(double));

	for(ii=0; ii<matrixDOF; ii++)
		solution[ii] = 1.0;

	if(myrank == 0)
		printf("[Poisson] Start solving equations.\n");

	cgsolver(matrixDOF, myrank, ncpus, poissonmatrix, rhs, solution, maxiteration, tolerance);

	if(myrank == 0)
	{
		printf("[Poisson] Solution obtained.\n");
		system("rm -rf result\n");
		system("mkdir result\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(myfilename,"./result/solution.dat");
	
	for(ii=0; ii<=ncpus; ii++)
	{
		if(myrank == ii)
		{
			fp = fopen(myfilename, "a");
			count = 0;
	
			for(jj=firstgrid_x; jj<=lastgrid_x; jj++)
			{
				for(kk=1; kk<=numgrid_y; kk++)
				{
					fprintf(fp,"%9.5f%10.5f%13.6f\n", coordinate[2*count], coordinate[2*count+1], solution[count]);
					count++;
				}
			}

			fclose(fp);
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	if(myrank == 0)
		printf("[Poisson] Solution printed.\n");

	free(solution);
	free(poissonmatrix);
	free(rhs);
	free(coordinate);

	if(myrank == 0)
		printf("[Poisson] Memory deallocated.\n");

	elapsed_time += MPI_Wtime();	

	if(myrank == 0)
		 printf("[Poisson] Finalizing MPI - Computing finished in %e (secs).\n", elapsed_time);

	MPI_Finalize();
	return 0;
}


