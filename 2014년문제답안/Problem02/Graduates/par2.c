#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<sys/time.h>
#include<mpi.h>


#define READY 0
#define WRITING -1
#define FINISH -991
#define NP_TAG 99
#define R_TAG 990
#define POT_TAG 123


#define max(a,b) (a)>(b)? (a):(b)
#define min(a,b) (a)<(b)? (a):(b)

long iseed = -9;
float ran2(long *);

typedef struct Pos{
	float x,y,z;
}Pos;

double potential(Pos *r, int np, int jstart, int jfinal){
	int i,j;
	double potent = 0;
	for(j=jstart;j<jfinal;j++){
		for(i=j+1;i<np;i++){
			float x= r[i].x-r[j].x;
			float y= r[i].y-r[j].y;
			float z= r[i].z-r[j].z;
			float dist = sqrtf(x*x+y*y+z*z);
			if(dist > 0.) potent += 1.L/dist;
		}
	}
	int jjstart = min(np-jfinal,np/2);
	int jjfinal = min(np/2,np - jstart);
	for(j=jjstart;j<jjfinal;j++){
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
	if(istep %3 ==0) {
		return (int)(100000*(ran2(&iseed)));
	}
	else {
		return (int)(500*(ran2(&iseed)));
	}
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


int myid,nid;

int main(int argc, char **argv){
	int i, j;
	int np,niter,maxnp=5000000;
	Pos *r;
	double totpotent=0;
	double potent;
	int finish = FINISH;
	MPI_Status mstatus;
	MPI_Request request;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nid);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	niter = 50*3;
	r = (Pos*)malloc(sizeof(Pos)*maxnp);
	ran2(&iseed);

	float time1, time2;
	time1 = gettime();

	long nlocalsize = 10000*1000;
	int jstart, jfinal;


	if(myid==0){
		int nsend,nrecv;
		int ready,src,dest;
		nsend = nrecv = 0;
		for(i=0;i<niter;i++){
			np =  getrandomnp(i,niter);
			for(j=0;j<np;j++){
				r[j].x = 2.*(ran2(&iseed))-1.;
				r[j].y = 2.*(ran2(&iseed))-1.;
				r[j].z = 2.*(ran2(&iseed))-1.;
			}
			long nwork = (long)np*((long)np-1)/2;
			long nsplit,njump;
			nsplit = (nwork-1)/nlocalsize + 1;
			if(nsplit <=1) {
				njump = np/2+1;
			}
			else {
				njump = (np/2-1)/nsplit+1;
			}
			/*
			if(np == 90583) {
				int kkkk = 1;
				while(kkkk) {
					kkkk = 1;
				}
			}
			*/
			for(j=np;j>np/2;j-=njump)
			{
				jstart = max(np/2,j-njump);
				jfinal = j;
				do {
					MPI_Probe(MPI_ANY_SOURCE, READY,MPI_COMM_WORLD, &mstatus);
					MPI_Recv(&ready,1, MPI_INT, mstatus.MPI_SOURCE, READY,MPI_COMM_WORLD,&mstatus);
					if(ready == READY){
						nsend ++;
						MPI_Send(&np,1, MPI_INT,mstatus.MPI_SOURCE, NP_TAG,MPI_COMM_WORLD);
						MPI_Send(r,np*sizeof(Pos), MPI_BYTE,mstatus.MPI_SOURCE, R_TAG,MPI_COMM_WORLD);
						MPI_Send(&jstart,sizeof(int), MPI_INT,mstatus.MPI_SOURCE, R_TAG,MPI_COMM_WORLD);
						MPI_Send(&jfinal,sizeof(int), MPI_INT,mstatus.MPI_SOURCE, R_TAG,MPI_COMM_WORLD);
/*
						printf("P%d is sending %d data to pid= %d with np= %d :    [%d , %d) && [%d , %d)\n",myid,i,mstatus.MPI_SOURCE,
								np,jstart,jfinal, min(np/2,np-jfinal), min(np/2,np-jstart));
*/
					}
					else {
						nrecv ++;
						MPI_Recv(&potent,1,MPI_DOUBLE,mstatus.MPI_SOURCE, POT_TAG,MPI_COMM_WORLD,&mstatus);
						totpotent += potent;
					}
				}while(ready != READY);
			}
		}
		j = 0;
		for(i=nrecv;i<nsend;){
			MPI_Probe(MPI_ANY_SOURCE, READY, MPI_COMM_WORLD, &mstatus);
			MPI_Recv(&ready, 1, MPI_INT, mstatus.MPI_SOURCE, READY, MPI_COMM_WORLD, &mstatus);
			if(ready == WRITING){
				MPI_Recv(&potent,1,MPI_DOUBLE,mstatus.MPI_SOURCE, POT_TAG,MPI_COMM_WORLD, &mstatus);
				totpotent += potent;
				MPI_Send(&finish,1, MPI_INT, mstatus.MPI_SOURCE, NP_TAG, MPI_COMM_WORLD);
				i++;
			}
			else {
				MPI_Send(&finish,1, MPI_INT, mstatus.MPI_SOURCE, NP_TAG, MPI_COMM_WORLD);
				j ++;
			}
		}
		for(i=1;i<nid-j;i++){
			MPI_Recv(&ready, 1, MPI_INT, MPI_ANY_SOURCE, READY, MPI_COMM_WORLD, &mstatus);
			MPI_Send(&finish, 1, MPI_INT, mstatus.MPI_SOURCE, NP_TAG, MPI_COMM_WORLD);
		}
	}
	else {
		int ready = READY;
		MPI_Send(&ready,1, MPI_INT, 0, READY, MPI_COMM_WORLD);
		MPI_Recv(&np, 1, MPI_INT, 0, NP_TAG, MPI_COMM_WORLD, &mstatus);
		while(np != FINISH){
			MPI_Recv(r, sizeof(Pos)*np, MPI_BYTE, 0, R_TAG, MPI_COMM_WORLD, &mstatus);
			MPI_Recv(&jstart, sizeof(int), MPI_INT, 0, R_TAG, MPI_COMM_WORLD, &mstatus);
			MPI_Recv(&jfinal, sizeof(int), MPI_INT, 0, R_TAG, MPI_COMM_WORLD, &mstatus);
			potent = potential(r, np,jstart,jfinal);
			ready = WRITING;
			MPI_Issend(&ready, 1, MPI_INT, 0, READY, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &mstatus);
			MPI_Send(&potent, 1, MPI_DOUBLE, 0, POT_TAG, MPI_COMM_WORLD);

			ready = READY;
			MPI_Send(&ready, 1, MPI_INT,0,READY, MPI_COMM_WORLD);
			MPI_Recv(&np, 1, MPI_INT, 0, NP_TAG, MPI_COMM_WORLD, &mstatus);
		}
	}
	time2 = gettime();
	if(myid==0) printf("[C] np=%d TPOT=%20.10f in Wtimes = %g (sec)\n",nid,totpotent, (time2-time1));
	MPI_Finalize();
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
			if (j < NTAB) iv[j] = *idum;
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
