#include<stdio.h>
#include<stdlib.h>
#include<math.h>
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
	n2 = 100000000;
	niter = 3;
	/* do not change ------ */

	ar = (double*) malloc(n2*sizeof(double));
	br = (double*) malloc(n2*sizeof(double));
	
	xi = 0.L;
	xf = 1.;
	dx = (xf-xi)/(double)(n2-n1-1);

	for(i=n1;i<n2;i++){
		br[i] = xi+(double)(i-n1)*dx;
	}

	MPI_Init(&argc, &argv);
	tic = MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &nid);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	para_range(n1,n2,nid,myid,&ista,&iend); // para range를 통해 작업 범위 나눔.
	printf("rank:%10d ista=%15d iend=%15d\n", myid, ista, iend);

	jsta = ista;
	jend = iend;
	if(myid==0) jsta = n1+1;
	if(myid == nid-1) jend = n2-1;

	// send/recv할때 보낼 rank에 사용
	inext = myid + 1; 
	iprev = myid - 1;
	if(myid == nid-1) inext = MPI_PROC_NULL;
	if(myid == 0) iprev = MPI_PROC_NULL;
	for(i=ista;i<iend;i++){
		br[i] = xi+ (double)(i-n1)*dx;
	}

	for(iter=0;iter<niter;iter++){
		itag = 101;
		/**
		 * 각 부분에서 idx-1과 idx+1부분이 필요하기 때문에 
		 * Isend/Irecv로 비동기적으로 보낸다.
		 * Wait를 통해 통신이 동작하는지 확인한다.
		 */
		MPI_Isend(br+iend-1,   1, MPI_DOUBLE, inext, itag, MPI_COMM_WORLD, &isd1); // inext에 b[j-1]전달
		MPI_Isend(br+ista,   1, MPI_DOUBLE, iprev, itag, MPI_COMM_WORLD, &isd2); // iprev b[j+1] 전달
		MPI_Irecv(br+ista-1,   1, MPI_DOUBLE, iprev, itag, MPI_COMM_WORLD, &irv1); // b[j-1] 받음
		MPI_Irecv(br+iend,   1, MPI_DOUBLE, inext, itag, MPI_COMM_WORLD, &irv2); // b[j+1] 받음
		MPI_Wait(&isd1,&istatus);
		MPI_Wait(&isd2,&istatus);
		MPI_Wait(&irv1,&istatus);
		MPI_Wait(&irv2,&istatus);

		for(j=jsta;j<jend;j++) {
			/*  not change -----{ */
			ar[j] = (br[j-1]+br[j+1])/4.L + br[j]/2.L + 1.L/genvv(br[j]);
			/*  not change -----} */
		}
		for(i=ista;i<iend;i++) {
			/*  not change -----{ */
			br[i] = ar[i];
			/*  not change -----} */
		}
	}
	ptmr = 0.L;
	for(j=jsta;j<jend;j++){
		ptmr += ar[j];
	}
	iroot = 0;
	MPI_Reduce(&ptmr, &tmr, 1, MPI_DOUBLE, MPI_SUM, iroot, MPI_COMM_WORLD); // MPI_Reduce로 ptmr 합침
	if(myid==0) printf("tmr = %16.6f\n",tmr);
	toc = MPI_Wtime();
	if(myid==0) printf("%g sec\n",toc-tic);

	free(ar);
	free(br); 

	MPI_Finalize();

}
