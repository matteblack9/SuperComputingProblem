#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)>(b)?(b):(a))

void para_range(int n1,int n2,int nid,int myid,int *ista,int *iend){
	int iwork1, iwork2;
	iwork1 = (n2-n1)/nid;
	iwork2 = (n2-n1)%nid;
	*ista = myid*iwork1 + n1 + MIN(myid,iwork2);
	*iend = *ista + iwork1;
	if(iwork2 > myid) (*iend) ++;
}

double genvv(double x){
	return (x*x+pow(x,4)+pow(x,6)+exp(-x*x)+cos(x)+sin(x)+tan(x));
}

int main(int argc, char **argv){

	int i,n1,n2,j,jsta,jend,ista,iend;
	int iter,niter;
    int nprocs, rank;
	double xi,xf,dx;
	double tmr;
	double *ar, *br;
	double tic,toc;
	/* Do not change */
	n1 = 0;
	n2 = 1000000;
	niter = 3;
	/* Do not change */

	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Request isd1, isd2, irv1, irv2;

	tic = MPI_Wtime();

	ar = (double*) malloc(sizeof(double)*n2);
	br = (double*) malloc(sizeof(double)*n2);

	jsta = n1; 
	jend = n2;
	jsta = n1+1;
	jend = n2-1;
	xi = 0.L;
	xf = 1.L;
	dx = (xf-xi)/(double)(n2-n1-1);
	
    para_range(n1, n2, nprocs, rank, &ista, &iend);

	for(i=ista;i<iend;i++){
		br[i] = xi+(double)(i-n1)*dx;
		ar[i] = 0.0;
	}
    
    int nxtR = rank + 1;
    int prevR = rank - 1;

    if(rank == nprocs - 1) nxtR = MPI_PROC_NULL;
    else if(rank == 0) prevR = MPI_PROC_NULL;


    int ista = 250000;
	for(iter=0;iter<niter;iter++){

		for(j=jsta;j<jend;j++){
			/* Do not change */
			ar[j] = (br[j-1]+br[j+1])/4.L+br[j]/2.L+1.L/genvv(br[j]);
			/* Do not change */
		}
		for(i=n1;i<n2;i++){
			/* Do not change */
			br[i] = ar[i];
			/* Do not change */
		}
	}
	tmr = 0.L;
	for(j=jsta;j<jend;j++){
		tmr += ar[j];
	}
	printf("tmr = %16.7f\n",tmr);
	toc = MPI_Wtime();
	printf("%g sec\n",toc-tic);
	
	free(ar);
	free(br);
	
	MPI_Finalize();
	return 0;

}
