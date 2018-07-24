#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include<mpi.h>


#define min(a,b) ((a)<(b)?(a):(b))


int ncall;

double simpson2(double (*func1)(double),int , double , double , int , int );

double func(double x){
	ncall ++;
	return 4.L/(x*x+1.L);
}

int main(int argc, char **argv){
	double rslt,aa,bb,rslt0;
	int n,ncall0;
	int ierr, kount, iroot;
	int myid,nproc;
	double time1,time2,time_start, time_end;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);

	if(myid ==0 && nproc > 1) printf("%d processes are alive\n",nproc);
	if(myid ==0 && nproc == 1) printf("%d processes is alive\n",nproc);


	time_start = MPI_Wtime();
	n=20;
	n =6; 
	n  = 200000;
	n = 2000000000;

	aa = 0.L;
	bb = 1.L;
	ncall = 0;
	ncall = 0;

	rslt = simpson2(func,n,aa,bb,myid,nproc);

	iroot = 0; kount = 1;
	MPI_Reduce(&rslt,&rslt0, kount, MPI_DOUBLE, MPI_SUM,iroot, MPI_COMM_WORLD);
	if(myid==0) printf("%d %18.14g n, rslt0\n",n,rslt0);
	iroot = 0; kount = 1;
	MPI_Reduce(&ncall, &ncall0,kount,MPI_INT,MPI_SUM,iroot, MPI_COMM_WORLD);
	if(myid==0) printf("%d ncall0\n",ncall0);

	time_end = MPI_Wtime();
	if(myid==0) {
		double timer = time_end - time_start;
		double timerm = timer/60.L;
		double timerh = timer/3600.L;
		double timerd = timer/3600.L/24.L;
		printf("%14.5g s %14.5g m %14.5g h %14.5g d\n",timer, timerm,timerh,timerd);
	}
	MPI_Finalize();
	return 0;
}

void equal_load(int n1, int n2, int nproc, int myid, int *istart, int *ifinish){
	int iw1,iw2;
	iw1 = (n2-n1+1)/nproc;
	iw2 = (n2-n1+1)%nproc;
	*istart = myid*iw1+n1+min(myid,iw2);
	*ifinish= *istart +iw1 -1;
	if(iw2>myid) *ifinish = *ifinish + 1;
	if(n2<*istart) *ifinish = *istart - 1;
}

double simpson2(double (*func1)(double),int n, double aa, double bb, int myid, int nproc){
	double rslt;
	double h,xx;
	int j,n1,n2;
	int istart,ifinish;

	rslt = 0.; 
	if((n%2) != 0) {
		fprintf(stderr,"input error, n must be even number %d\n",n);
		exit(9);
	}

	n1 = 1; n2 = n-1;
	equal_load(n1,n2,nproc,myid,&istart,&ifinish);


	h = (bb-aa)/(double)n;
	if(myid==0) rslt = (func1(aa)+func1(bb));
	for(j=istart;j<=ifinish;j++){
		xx = aa + h*(double)j;
		/*
		if((j%2) == 1) {
		*/
		if((j&1)) {
			rslt += 4.L*func1(xx);
		}
		else {
			rslt += 2.L*func1(xx);
		}
	}
	rslt = rslt * h/3.L;
	return rslt;
}


