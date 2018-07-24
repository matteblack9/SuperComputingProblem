#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixconstructor.h"
#include "cgsolver.h"

#define PI 3.141592653589793

int main(int argc, char **argv) {

	int numgrid_x, numgrid_y;
	int matrixDOF;  
	int count, ii, jj, kk, maxiteration;  
	double length_x, length_y, gridsize, tolerance;
	double *poissonmatrix, *rhs, *solution, *coordinate;
	char myfilename[256];
	FILE *fp;

	maxiteration = 20000;
	tolerance = 1.0e-10;
	gridsize = 0.01;
	length_x = 1.0;
	length_y = 1.0;

	numgrid_x = (int)round(length_x/gridsize) + 1;
	numgrid_y = (int)round(length_y/gridsize) + 1;
	matrixDOF = numgrid_x*numgrid_y;

	printf("[Poisson] Geometry and matrix size initialized.\n");

	coordinate = (double*)malloc(matrixDOF*2*sizeof(double));
	rhs = (double*)malloc(matrixDOF*sizeof(double));

	count = 0;

	for(ii=1; ii<=numgrid_x; ii++)
	{
		for(jj=1; jj<=numgrid_y; jj++)
		{
			coordinate[2*count] = (ii-1)*gridsize;
			coordinate[2*count+1] = (jj-1)*gridsize;
			rhs[count] = sin(coordinate[2*count]/length_x*PI)*sin(coordinate[2*count+1]/length_y*PI)*gridsize*gridsize;
			count++;
		} 
	}
	
	printf("[Poisson] Geometry and rhs constructed.\n");

	poissonmatrix = (double*) malloc(matrixDOF*matrixDOF*sizeof(double)); 

	for(ii=0; ii<matrixDOF*matrixDOF; ii++)	
		poissonmatrix[ii] = 0.0;

	construct_poissonmatrix(1, 1, numgrid_x, numgrid_x, numgrid_y, poissonmatrix);

	printf("[Poisson] Poisson matrix constructed.\n");

	solution = (double*) malloc(matrixDOF*sizeof(double));

	for(ii=0; ii<matrixDOF; ii++)
		solution[ii] = 1.0;

	printf("[Poisson] Start solving equations.\n");

	cgsolver(matrixDOF, poissonmatrix, rhs, solution, maxiteration, tolerance);

	printf("[Poisson] Solution obtained.\n");

	system("rm -rf result\n");
	system("mkdir result\n");
	sprintf(myfilename,"./result/solution.dat");
	fp = fopen(myfilename, "wt");
	
	count = 0;
	for(ii=1; ii<=numgrid_x; ii++)
	{
		for(jj=1; jj<=numgrid_y; jj++)
		{
			fprintf(fp,"%9.5f%10.5f%13.6f\n", coordinate[2*count], coordinate[2*count+1], solution[count]);
			count++;
		}
	}

	fclose(fp);

	printf("[Poisson] Solution printed.\n");

	free(solution);
	free(poissonmatrix);
	free(rhs);
	free(coordinate);

	printf("[Poisson] Memory deallocated.\n");
	return 0;
}


