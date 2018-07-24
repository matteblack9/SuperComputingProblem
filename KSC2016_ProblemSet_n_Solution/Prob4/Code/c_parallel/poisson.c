#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "matrixconstructor.h"
#include "cgsolver.h"
#include "multigrid.h"

#define PI 3.141592653589793

int main(int argc, char **argv) 
{
    int nlevel, initial_ngrid;
    int *numgrid, *globalnumgrid, *matrixDOF, *globalmatrixDOF;
    int *firstgrid_x, *lastgrid_x, *firstrow, *lastrow;
    int count, ii, jj, kk, ll, maxiteration, myrank, ncpus;  
    double length, tolerance;
    double *gridsize;
    double **poissonmatrix, **rhs, **solution, **coordinate;
    double time_start, time_end;
    double elapsed_time;
    char myfilename[256];
    FILE *fp;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &ncpus);

    if(myrank==0)
	printf("[Poisson] MPI initialized.\n");

    nlevel = 3;
    initial_ngrid = 128;
    maxiteration = 20000;
    tolerance = 1.0e-10;
    length = 1.0;

    time_start = MPI_Wtime();

    numgrid = (int*)malloc(nlevel*sizeof(int));
    globalnumgrid = (int*)malloc(nlevel*sizeof(int));
    matrixDOF = (int*)malloc(nlevel*sizeof(int));
    globalmatrixDOF = (int*)malloc(nlevel*sizeof(int));
    firstgrid_x = (int*)malloc(nlevel*sizeof(int));
    lastgrid_x = (int*)malloc(nlevel*sizeof(int));
    firstrow = (int*)malloc(nlevel*sizeof(int));
    lastrow = (int*)malloc(nlevel*sizeof(int));
    gridsize = (double*)malloc(nlevel*sizeof(double));

    globalnumgrid[nlevel-1] = initial_ngrid;
    gridsize[nlevel-1] = length/initial_ngrid;
    globalmatrixDOF[nlevel-1] = globalnumgrid[nlevel-1]*globalnumgrid[nlevel-1];
        
    firstgrid_x[nlevel-1] = (globalnumgrid[nlevel-1]*myrank/ncpus)+1;
    lastgrid_x[nlevel-1] = globalnumgrid[nlevel-1]*(myrank+1)/ncpus;
    numgrid[nlevel-1] = lastgrid_x[nlevel-1]-firstgrid_x[nlevel-1]+1;
    matrixDOF[nlevel-1] = numgrid[nlevel-1]*globalnumgrid[nlevel-1];
    firstrow[nlevel-1] = (firstgrid_x[nlevel-1]-1)*globalnumgrid[nlevel-1]+1;
    lastrow[nlevel-1] = lastgrid_x[nlevel-1]*globalnumgrid[nlevel-1];

    for(ii=nlevel-1; ii>=1; ii--)
    {
	globalnumgrid[ii-1] = globalnumgrid[ii]/2;
	gridsize[ii-1] = length/globalnumgrid[ii-1];
	globalmatrixDOF[ii-1] = globalnumgrid[ii-1]*globalnumgrid[ii-1];

    	firstgrid_x[ii-1] = (globalnumgrid[ii-1]*myrank/ncpus)+1;
    	lastgrid_x[ii-1] = globalnumgrid[ii-1]*(myrank+1)/ncpus;
    	numgrid[ii-1] = lastgrid_x[ii-1]-firstgrid_x[ii-1]+1;
    	matrixDOF[ii-1] = numgrid[ii-1]*globalnumgrid[ii-1];
    	firstrow[ii-1] = (firstgrid_x[ii-1]-1)*globalnumgrid[ii-1]+1;
    	lastrow[ii-1] = lastgrid_x[ii-1]*globalnumgrid[ii-1];
    }

    if(globalnumgrid[0]<ncpus || (nlevel>1 && globalnumgrid[0]%ncpus!=0))
    {
 	if(globalnumgrid[0]<ncpus && myrank==0)	
		printf("[Poisson] # of MPI processes must be equal or larger than # of grids in the (smallest) domain!\n"); 
	     
	if(nlevel>1 && globalnumgrid[0]%ncpus!=0 && myrank==0)
		printf("[Poisson/Multigrid] # of grids in the smallest domain must be a multiple of # of MPI processes!\n");

	free(gridsize);
	free(globalnumgrid);
	free(globalmatrixDOF);
	free(matrixDOF);
	free(numgrid);
	free(firstrow);
	free(lastrow);
	free(firstgrid_x);
	free(lastgrid_x);
	
	MPI_Finalize();

       if(myrank==0)
          printf("[Poisson] MPI finalized.\n");

	return 0;
    }

    if(myrank==0)
    	printf("[Poisson] Geometry and matrix size initialized.\n");

    poissonmatrix = (double**)malloc(nlevel*sizeof(double*));
    rhs = (double**)malloc(nlevel*sizeof(double*));
    solution = (double**)malloc(nlevel*sizeof(double*));
    coordinate = (double**)malloc(nlevel*sizeof(double*));

    for(ii=0; ii<nlevel; ii++) 
    {
	poissonmatrix[ii] = (double*)malloc(matrixDOF[ii]*globalmatrixDOF[ii]*sizeof(double));
	rhs[ii] = (double*)malloc(matrixDOF[ii]*sizeof(double));
	solution[ii] = (double*)malloc(matrixDOF[ii]*sizeof(double));
	coordinate[ii] = (double*)malloc(matrixDOF[ii]*2*sizeof(double));
    }

    for(ii=0; ii<nlevel; ii++)
    {
	count = 0;

	for(jj=firstgrid_x[ii]; jj<=lastgrid_x[ii]; jj++)
	{
	    for(kk=1; kk<=globalnumgrid[ii]; kk++)
	    {
		coordinate[ii][2*count] = jj*gridsize[ii]-0.5*gridsize[ii];
		coordinate[ii][2*count+1] = kk*gridsize[ii]-0.5*gridsize[ii];
		count++;
	    } 
	}
    }

    count = 0;

    for(ii=1; ii<=numgrid[nlevel-1]; ii++)
    {
	for(jj=1; jj<=globalnumgrid[nlevel-1]; jj++)
	{
	    rhs[nlevel-1][count] = sin(coordinate[nlevel-1][2*count]/length*PI)* \
				   sin(coordinate[nlevel-1][2*count+1]/length*PI)*gridsize[nlevel-1]*gridsize[nlevel-1];
	    count++;
	} 
    }

    if(myrank==0)
    	printf("[Poisson] Geometry and rhs constructed.\n");

    for(ii=0; ii<nlevel; ii++) 
	for(jj=0; jj<matrixDOF[ii]*globalmatrixDOF[ii]; jj++)	
	    poissonmatrix[ii][jj] = 0.0;
 
    for(ii=0; ii<nlevel; ii++)
	construct_poissonmatrix(firstrow[ii], firstgrid_x[ii], lastgrid_x[ii], globalnumgrid[ii], globalnumgrid[ii], poissonmatrix[ii]);



    if(myrank==0)
    	printf("[Poisson] Poisson matrix constructed.\n");

    for(ii=0; ii<matrixDOF[nlevel-1]; ii++)
	solution[nlevel-1][ii] = 1.0;

    if(myrank==0)
    	printf("[Poisson] Start solving equations.\n");

    for(ii=nlevel-1; ii>0; ii--) 
    {
	restriction(rhs[ii], rhs[ii-1], numgrid[ii-1], globalnumgrid[ii-1]);
	
	if(myrank==0)
		printf("[Multigrid] Restriction of RHS done (level %d -> level %d).\n", ii, ii-1);

	for(jj=0; jj<matrixDOF[ii-1]; jj++)
	    solution[ii-1][jj] = 1.0;

//        for(int aa=0; aa<numgrid[ii-1]*globalnumgrid[ii-1]; aa++)
//                printf("%f %f\n", rhs[ii-1][aa], solution[ii-1][aa]);
    }

    cgsolver(matrixDOF[0], myrank, ncpus, poissonmatrix[0], rhs[0], solution[0], maxiteration, tolerance);

    for(ii=1; ii<nlevel; ii++)
    {
	interpolation(solution[ii-1], solution[ii], numgrid[ii], globalnumgrid[ii], myrank, ncpus);
	if(myrank==0)
		printf("[Multigrid] Interpolation of solution done (level %d -> level %d).\n", ii-1, ii);
	cgsolver(matrixDOF[ii], myrank, ncpus, poissonmatrix[ii], rhs[ii], solution[ii], maxiteration, tolerance);
	if(myrank==0)
		printf("[Poisson] Solution in level %d obtained.\n", ii);
    }

    if(myrank==0)
    	printf("[Poisson] Final solution obtained.\n");

    if(myrank==0)
    {
	system("rm -rf result\n");
    	system("mkdir result\n");
     }

    MPI_Barrier(MPI_COMM_WORLD);

    for(ii=0; ii<nlevel; ii++) 
    {
	sprintf(myfilename,"./result/solution_level%d.dat", ii);

	for(jj=0; jj<ncpus; jj++)
	{
		if(myrank==jj)
		{
			fp = fopen(myfilename, "a");
			count = 0;

			for(kk=firstgrid_x[ii]; kk<=lastgrid_x[ii]; kk++)
			{
				for(ll=1; ll<=globalnumgrid[ii]; ll++)
				{
					fprintf(fp,"%9.5f%10.5f%13.6f\n", coordinate[ii][2*count], coordinate[ii][2*count+1], solution[ii][count]);
					count++;
	    			}
			}

			fclose(fp);
    		}
		
		MPI_Barrier(MPI_COMM_WORLD);
	}
    }
  
    if(myrank==0)
    	printf("[Poisson] Final solution printed.\n");

    for(ii=0; ii<nlevel; ii++) 
    {
	free(coordinate[ii]);
	free(solution[ii]);
	free(rhs[ii]);
	free(poissonmatrix[ii]);
    }

    free(poissonmatrix);
    free(rhs);
    free(solution);
    free(coordinate);

    free(gridsize);
    free(globalnumgrid);
    free(globalmatrixDOF);
    free(matrixDOF);
    free(numgrid);
    free(firstrow);
    free(lastrow);
    free(firstgrid_x);
    free(lastgrid_x);

    if(myrank==0)
    	printf("[Poisson] Memory deallocated.\n");

    time_end = MPI_Wtime();
    elapsed_time = time_end - time_start;

    if(myrank==0)
    	printf("[Poisson] Mission completed in %f (secs).\n", elapsed_time);

    MPI_Finalize();

    if(myrank==0)
	printf("[Poisson] MPI finalized.\n");
    
    return 0;
}

