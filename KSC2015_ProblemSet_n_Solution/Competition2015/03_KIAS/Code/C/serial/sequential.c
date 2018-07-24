#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<sys/time.h>

float ran2(long *);
long iseed=-9;

struct timeval tv;
float gettime(){
	static int startflag = 1;
	static double tsecs0, tsecs1;
	if(startflag) {
		(void ) gettimeofday(&tv, NULL);
		tsecs0 = tv.tv_sec + tv.tv_usec*1.0E-6;
		startflag = 0;
	}
	(void) gettimeofday(&tv, NULL);
	tsecs1 = tv.tv_sec + tv.tv_usec*1.0e-6;
	return (float) (tsecs1 - tsecs0);
}



int main(int argc, char **argv){
	long i, j;
	long np,niter,maxnp=100000000;
	float *val,**lval;

	niter = atoi(argv[1]);

	float time1, time2;
	time1 = gettime();


	val = (float*)malloc(sizeof(float)*maxnp*niter);
	lval = (float**)malloc(sizeof(float*)*niter);
	long nlval[niter];

	float step = (1.L-0.L)/(float)niter;

	for(i=0;i<niter;i++){
		lval[i] = (float*)malloc(sizeof(float)*maxnp*2);
		nlval[i] = 0;
	}
	np = 0;
	for(i=0;i<niter;i++){
		iseed = -1*(i+1284L);
		for(j=0;j<maxnp;j++){
			val[np++] = ran2(&iseed);
		}
	}
	for(i=0;i<np;i++){
		j = val[i]/step;
		*(lval[j]+nlval[j]) = val[i];
		nlval[j]++;
	}
	for(i=0;i<niter;i++){
		printf("p%d has %ld members ::: %g %g\n", (int)i,nlval[i], step*i, step*(+1));
	}
	time2 = gettime();
	printf("Wallclock time = %g second\n", (time2-time1));
}
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) {
				iv[j] = *idum;
			}
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */
