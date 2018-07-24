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
	if(istep < niter/10) {
		return (int)(3000+5000*(rand()/(RAND_MAX+1.)));
	}
	else {
		return (int)(100*(rand()/(RAND_MAX+1.)));
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
	unsigned int iseed = 100;
	MPI_Status mstatus;
	MPI_Request request;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nid);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	niter = 4000;
	r = (Pos*)malloc(sizeof(Pos)*maxnp);
	srand( iseed);

	float time1, time2;
	time1 = gettime();



	if(myid==0){
		int nsend,nrecv;
		int ready,src,dest;
		nsend = nrecv = 0;
		for(i=0;i<niter;i++){
			np =  getrandomnp(i,niter);
			for(j=0;j<np;j++){
				r[j].x = 2.*(rand()/(RAND_MAX+1.))-1.;
				r[j].y = 2.*(rand()/(RAND_MAX+1.))-1.;
				r[j].z = 2.*(rand()/(RAND_MAX+1.))-1.;
			}
			do {
				MPI_Probe(MPI_ANY_SOURCE, READY,MPI_COMM_WORLD, &mstatus);
				MPI_Recv(&ready,1, MPI_INT, mstatus.MPI_SOURCE, READY,MPI_COMM_WORLD,&mstatus);
				if(ready == READY){
					nsend ++;
					MPI_Send(&np,1, MPI_INT,mstatus.MPI_SOURCE, NP_TAG,MPI_COMM_WORLD);
					MPI_Send(r,np*sizeof(Pos), MPI_BYTE,mstatus.MPI_SOURCE, R_TAG,MPI_COMM_WORLD);
				}
				else {
					nrecv ++;
					MPI_Recv(&potent,1,MPI_DOUBLE,mstatus.MPI_SOURCE, POT_TAG,MPI_COMM_WORLD,&mstatus);
					totpotent += potent;
				}
			}while(ready != READY);
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
			potent = potential(r, np);
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
	if(myid==0) printf("Total potential is %20.10g in wallclock time = %g second\n",totpotent, (time2-time1));
	MPI_Finalize();
}
