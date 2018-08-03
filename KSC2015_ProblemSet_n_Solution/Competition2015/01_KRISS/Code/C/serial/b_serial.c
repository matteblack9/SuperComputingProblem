#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

double genvv(double x){
	return (x*x+pow(x,4)+pow(x,6)+exp(-x*x)+cos(x)+sin(x)+tan(x));
}

int main(int argc, char **argv){

	int i,n1,n2,j,jsta,jend;
	int iter,niter;
	double xi,xf,dx;
	double tmr;
	double *ar, *br;
	double tic,toc;
	/* Do not change */
	n1 = 0;
<<<<<<< HEAD
	n2 = 8;
=======
	n2 = 100000000; // 병렬화 해야하는 부분.
>>>>>>> a52bcd65976077d8d4ed396d7e5eec526700fb10
	niter = 3;
	/* Do not change */

	MPI_Init(&argc, &argv);
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

	for(i=n1;i<n2;i++){
		br[i] = xi+(double)(i-n1)*dx;
		ar[i] = 0.0;
	}

	int ista = 250000;

	for(iter=0;iter<niter;iter++){

		for(i=n1;i<n2;i++){
			printf("before br[%d] = %lf \n", i, br[i]);
		}

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

		for(i=n1;i<n2;i++){
			printf("After br[%d] = %lf \n", i, br[i]);
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
