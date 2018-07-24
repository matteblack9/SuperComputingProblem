#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>

#include<omp.h>


#define min(a,b) ((a)<(b)?(a):(b))


int ncall0;

double simpson2(double (*func1)(double),int , double , double);

double func(double x){
	return 4.L/(x*x+1.L);
}

float gettime();
int main(int argc, char **argv){
	double rslt,aa,bb,rslt0;
	int n;
	int ierr, kount, iroot;
	double time1,time2,time_start, time_end;


	time_start = gettime();
	time_start = omp_get_wtime();
	n=20;
	n =6; 
	n  = 200000;
	n = 2000000000;

	aa = 0.L;
	bb = 1.L;

	rslt = simpson2(func,n,aa,bb);

	printf("%d %18.14g n, rslt0\n",n,rslt);
	printf("%d ncall0\n",ncall0);

	time_end = gettime();
	time_end = omp_get_wtime();
	{
		double timer = time_end - time_start;
		double timerm = timer/60.L;
		double timerh = timer/3600.L;
		double timerd = timer/3600.L/24.L;
		printf("%14.5g s %14.5g m %14.5g h %14.5g d\n",timer, timerm,timerh,timerd);
	}
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

double simpson2(double (*func1)(double),int n, double aa, double bb){
	int myid,nproc;
	double rslt,trslt=0;
	double h,xx;
	int j,n1,n2;
	int istart,ifinish;

	if((n%2) != 0) {
		fprintf(stderr,"input error, n must be even number %d\n",n);
		exit(9);
	}

	n1 = 1; n2 = n-1;
	h = (bb-aa)/(double)n;

#pragma omp parallel
	{
		nproc = omp_get_num_threads();
		myid = omp_get_thread_num();
		if(myid==0) printf("Total number of threads is %d\n",nproc);
	}


	int ncall ;
	ncall0 = 0;



#ifdef OLD
#pragma omp parallel private(nproc,myid,istart,ifinish,rslt,xx,j,ncall)
	{
		rslt = 0.; 
		nproc = omp_get_num_threads();
		myid = omp_get_thread_num();
		equal_load(n1,n2,nproc,myid,&istart,&ifinish);
		ncall = 0;
		if(myid==0) {
			rslt = (func1(aa)+func1(bb));
			ncall += 2;
		}
		for(j=istart;j<=ifinish;j++){
			xx = aa + h*(double)j;
			if((j%2) == 1) {
				rslt += 4.L*func1(xx);
			}
			else {
				rslt += 2.L*func1(xx);
			}
			ncall ++;
		}
		rslt = rslt * h/3.L;
#pragma omp critical
		{
			trslt += rslt;
			ncall0 += ncall;
		}
	}
#else
#pragma omp parallel private(rslt, ncall)
	{
		rslt = 0;
		ncall = 0;
#pragma omp for private(j,xx)
		for(j=1;j<=n-1;j++){
			xx = aa + h*(double)j;
			if((j%2)==1) rslt += 4.L*func1(xx);
			else rslt += 2.L*func1(xx);
			ncall ++;
		}
#pragma omp critical
		{
			trslt += rslt;
			ncall0 += ncall;
		}
	}
	trslt += func1(aa)+func1(bb);
	ncall0 += 2;
	trslt = trslt* h/3.L;
#endif
	return trslt;

}

struct timeval tv;

float gettime()
{
	static int startflag=1;
	static double tsecs0, tsecs1;
	
	if( startflag ) {
		(void) gettimeofday(&tv,0);
		tsecs0 = tv.tv_sec + tv.tv_usec*1.0e-6;
		startflag = 0;
	}
	(void) gettimeofday(&tv,0);
	tsecs1 = tv.tv_sec + tv.tv_usec*1.0e-6;

	return ((float) (tsecs1-tsecs0));

}
