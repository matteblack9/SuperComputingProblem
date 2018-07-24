[문제1번 C 기본 코드]

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>

int ncall;
typedef enum { FALSE=0, TRUE=1} boolean;
float gettime();

double func(double x){
	ncall ++;
	return 4.L/(x*x+1.L);
}

double simpson0(double (*func1)(double),int n,double aa,double bb){
	double h,xx,rslt;
	boolean lodd;
	int j;

	if((n%2) !=0){
		fprintf(stderr,"input error, n must be even number %d\n",n);
		exit(9);
	}
	rslt = 0.L;
	h = (bb-aa)/(double)n;
	rslt = (func1(aa) + func1(bb));
	for(lodd = TRUE,j = 1; j<n;j++){
		xx = aa + h * (double)j;
		if(lodd) {
			rslt +=  4.L*func1(xx);
		}
		else {
			rslt += 2.L*func1(xx);
		}
		lodd = ((lodd)? (FALSE):(TRUE));
	}
	rslt = rslt*h/3.L;
	return rslt;
}

int main(int argc, char **argv){
	double rslt,aa,bb;
	int n;
	double time_start,time_end;

	time_start = gettime();
	n = 2000000000;
	aa = 0.L;
	bb = 1.L;
	ncall = 0;

	rslt = simpson0(func,n,aa,bb);
	time_end = gettime();
	printf("%d %18.16g n, rslt,%d\n",n,rslt,ncall);
	printf("Wallclock time %g\n",time_end-time_start);
	return 0;
}
struct timeval tv;
float gettime()
{
	static int startflag=1;
	static double tsecs0, tsecs1;	
	if( startflag ) {
		(void ) gettimeofday(&tv,0);
		tsecs0 = tv.tv_sec + tv.tv_usec*1.0e-6;
		startflag = 0;
	}
	(void) gettimeofday(&tv,0);
	tsecs1 = tv.tv_sec + tv.tv_usec*1.0e-6;
	return ((float) (tsecs1-tsecs0));
}
