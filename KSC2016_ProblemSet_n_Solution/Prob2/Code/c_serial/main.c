#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "multigrid.h"

#define PI 3.141592653589793

int main(int argc, char **argv) 
{
	int nlevel = 15;
	double length_x = 1.0, length_y = 1.0;
    int maxiteration = 20;
	int ngrid;
	int numgrid_x[nlevel], numgrid_y[nlevel];
	int matrixDOF[nlevel];  
	int count, ii, jj, kk, i, iter;  
	double gridsize[nlevel];
	double *rhs[nlevel], *solution[nlevel];
	double coord_x, coord_y;
	clock_t time_start, time_end;
	double time_elapsed;

	char myfilename[256];
	FILE *fp;

	time_start=clock();

	ngrid = pow(2,nlevel-1);

	numgrid_x[nlevel-1] = ngrid;
	numgrid_y[nlevel-1] = ngrid;
	gridsize[nlevel-1] = length_x/(double)ngrid;
	matrixDOF[nlevel-1] = numgrid_x[nlevel-1]*numgrid_y[nlevel-1];

	for(i=nlevel-1; i>=1; i--)
	{
		numgrid_x[i-1] = numgrid_x[i]/2;
		numgrid_y[i-1] = numgrid_y[i]/2;
		gridsize[i-1] = length_x/numgrid_x[i-1];
		matrixDOF[i-1] = numgrid_x[i-1]*numgrid_y[i-1];
	}

	printf("[Multigrid] Geometry and matrix size initialized.\n");

	for(i=0; i<nlevel; i++) 
	{
		rhs[i] = (double*)malloc(matrixDOF[i]*sizeof(double));
		solution[i] = (double*)malloc(matrixDOF[i]*sizeof(double));
	}

	count = 0;

	for(ii=0; ii<numgrid_x[nlevel-1]; ii++)
	{
		coord_x = ii*gridsize[nlevel-1]+0.5*gridsize[nlevel-1];
		for(jj=0; jj<numgrid_y[nlevel-1]; jj++)
		{
			coord_y = jj*gridsize[nlevel-1]+0.5*gridsize[nlevel-1];
			rhs[nlevel-1][count] = sin(coord_x/length_x*PI)* \
					sin(coord_y/length_y*PI);
			count++;
		} 
	}
	printf("[Multigrid] Geometry and rhs constructed.\n");
	printf("[Multigrid] Start solving equations.\n");

	for(iter=0;iter<maxiteration;iter++) 
	{
		printf("[Multigrid] iteration %d \n",iter+1);
		for(i=nlevel-1; i>0; i--) 
		{
			restriction(rhs[i],rhs[i-1],numgrid_x[i-1]);
			printf("[Multigrid] Restriction of RHS done (level %d -> level %d).\n", i, i-1);
		}

		solution[0][0] = -rhs[0][0]/4.0;
		printf("[Multigrid] Solution at the coarsest level = %f\n",solution[0][0]);
	
		for(i=1; i<nlevel; i++)
		{
			interpolation(solution[i-1],solution[i],numgrid_x[i]);
			printf("[Multigrid] Interpolation of solution done (level %d -> level %d).\n", i-1, i);
		}
		for(ii=0;ii<matrixDOF[nlevel-1];ii++)
		{
			rhs[nlevel-1][ii] = PI*PI*PI*PI*solution[nlevel-1][ii];
		}
	}

	printf("[Multigrid] Final solution obtained.\n");

	i=(nlevel-1)-5;
	sprintf(myfilename,"./solution.dat");
	fp = fopen(myfilename, "wt");

	for(ii=-8; ii<8; ii++)
	{
		coord_x = (ii+numgrid_x[i]/2)*gridsize[i]+0.5*gridsize[i];
		for(jj=-8; jj < 8; jj++)
		{
			coord_y = (jj+numgrid_y[i]/2)*gridsize[i]+0.5*gridsize[i];
			count = (ii+numgrid_x[i]/2)*numgrid_y[i] + jj+numgrid_y[i]/2;
			fprintf(fp,"%12.6f %12.6f %12.6f\n", coord_x, coord_y, solution[i][count]);
		}
	}

	fclose(fp);

	printf("[Multigrid] Final solution printed.\n");

	for(i = 0; i < nlevel; i++) 
	{
		free(solution[i]);
		free(rhs[i]);
	}

	printf("[Multigrid] Memory deallocated.\n");
	time_end=clock();
	time_elapsed = (double)(time_end - time_start)/CLOCKS_PER_SEC;

	printf("[Multigrid] Mission completed in %f (secs).\n", time_elapsed);

	return 0;
}

