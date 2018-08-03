#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)>(b)?(b):(a))

#define TRUE 1
#define FALSE 0



void para_range(int n1,int n2,int nid,int myid,int *ista,int *iend){
	int iwork1, iwork2;
	iwork1 = (n2-n1)/nid;
	iwork2 = (n2-n1)%nid;
	*ista = myid*iwork1 + n1 + MIN(myid,iwork2);
	*iend = *ista + iwork1;
	if(iwork2 > myid) (*iend) ++;
}

double genvv(double x)
{
	double res = x*x + pow(x,4)+ pow(x,6)+exp(-x*x)+cos(x)+sin(x)+tan(x);
	return res;

}

int main (int argc, char **argv)
{
	int i,n1,n2,j,jsta,jend;
	int iter,niter;
	MPI_Status istatus;
	int ierr, myid,nid;
	int iprev, inext, ista, iend;
	MPI_Request isd1,isd2,irv1,irv2;
	int itag, iroot;
	double xi,xf,dx;
	double tmr;
	double *ar, *br;
	double ptmr, tic,toc;

	/* do not change ------ */
	n1 = 0;
	n2 = 8;
	niter = 3;
	/* do not change ------ */

	ar = (double*) malloc(n2*sizeof(double));
	br = (double*) malloc(n2*sizeof(double));
	
	xi = 0.L;
	xf = 1.;
	dx = (xf-xi)/(double)(n2-n1-1);

	MPI_Init(&argc, &argv);
	tic = MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &nid);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	para_range(n1,n2,nid,myid,&ista,&iend);
	printf("rank:%10d ista=%15d iend=%15d\n", myid, ista, iend);

	jsta = ista;
	jend = iend;
	if(myid==0) jsta = n1+1;
	if(myid == nid-1) jend = n2-1;

	inext = myid + 1;
	iprev = myid - 1;
	if(myid == nid-1) inext = MPI_PROC_NULL;
	if(myid == 0) iprev = MPI_PROC_NULL;

	
	for(i=ista;i<iend;i++){
		br[i] = xi+ (double)(i-n1)*dx;
		ar[i] = 0.0;
	}

	for(iter=0;iter<niter;iter++){
		itag = 101;
		MPI_Isend(br+iend-1,   1, MPI_DOUBLE, inext, itag, MPI_COMM_WORLD, &isd1);
		MPI_Isend(br+ista,   1, MPI_DOUBLE, iprev, itag, MPI_COMM_WORLD, &isd2);
		MPI_Irecv(br+ista-1,   1, MPI_DOUBLE, iprev, itag, MPI_COMM_WORLD, &irv1);
		MPI_Irecv(br+iend,   1, MPI_DOUBLE, inext, itag, MPI_COMM_WORLD, &irv2);
		MPI_Wait(&isd1,&istatus);
		MPI_Wait(&isd2,&istatus);
		MPI_Wait(&irv1,&istatus);
		MPI_Wait(&irv2,&istatus);

		for(i=ista;i<iend;i++){
			printf("before iter = %d rank = %d br[%d] = %lf \n", iter, myid, i, br[i]);
		}

		for(j=jsta;j<jend;j++){
			/*  not change -----{ */
			ar[j] = (br[j-1]+br[j+1])/4.L + br[j]/2.L + 1.L/genvv(br[j]);
			/*  not change -----} */
		}
		for(i=ista;i<iend;i++){
			/*  not change -----{ */
			br[i] = ar[i];
			/*  not change -----} */
		}

		for(i=ista;i<iend;i++){
			printf("after iter = %d rank = %d br[%d] = %lf \n", iter, myid, i, br[i]);
		}

	}
	ptmr = 0.L;
	for(j=jsta;j<jend;j++){
		ptmr += ar[j];
	}

	printf("rank = %d ptmr = %f \n", myid, ptmr);

	iroot = 0;
	MPI_Reduce(&ptmr, &tmr, 1, MPI_DOUBLE, MPI_SUM, iroot, MPI_COMM_WORLD);
	if(myid==0) printf("tmr = %16.6f\n",tmr);
	toc = MPI_Wtime();
	if(myid==0) printf("%g sec\n",toc-tic);

	free(ar);
	free(br); 

	MPI_Finalize();

}
