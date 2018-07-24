#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<sys/time.h>




typedef struct Pos{
	float x,y,z;
}Pos;

double potential(Pos *r, int np){
	int i,j;
	double potent = 0;
	for(j=0;j<np;j++){
		for(i=j+1;i<np;i++){
			float x= r[i].x-r[j].x;
			float y= r[i].y-r[j].y;
			float z= r[i].z-r[j].z;
			float dist = sqrtf(x*x+y*y+z*z);
			if(dist > 0.) potent += 1./dist;
		}
	}
	return potent;
}

int getrandomnp(int istep, int niter){
	if(istep%3 ==0) return (int)(50000*(rand()/(RAND_MAX+1.)));
	else  return (int)(500*(rand()/(RAND_MAX+1.)));
}

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
	int i, j;
	int np,niter,maxnp=5000000;
	Pos *r;
	double totpotent=0;
	unsigned int iseed = 100;

	niter = 50;
	r = (Pos*)malloc(sizeof(Pos)*maxnp);
	srand( iseed);

	float time1, time2;
	time1 = gettime();



	for(i=0;i<niter;i++){
		np =  getrandomnp(i,niter);
		for(j=0;j<np;j++){
			r[j].x = 2.*(rand()/(RAND_MAX+1.))-1.;
			r[j].y = 2.*(rand()/(RAND_MAX+1.))-1.;
			r[j].z = 2.*(rand()/(RAND_MAX+1.))-1.;
		}
		totpotent += potential(r, np);
	}
	time2 = gettime();
	printf("Total potential is %20.10g in wallclock time = %g second\n",totpotent, (time2-time1));
}
