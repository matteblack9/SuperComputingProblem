#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "kmeans.h"

/**
* kindex와 pt가 N_PT만큼 할당됨.
*/
void assignment_step(const PPOINT kmeans, const PPOINT pt,
		const int istart, const int iend, int *kindex)
{
	double d1, d2;

	for(int i=istart; i<iend; i++) {
		kindex[i] = 0;
		d1 = distance(pt[i], kmeans[0]);
		for(int j=1; j<N_K; j++) {
			d2 = distance(pt[i], kmeans[j]);
			if( d1 > d2 ) {
				d1 = d2;
				kindex[i] = j;
			}
		}
	}
}

/**
* 가장 크기가 큰 N_PT를 분할하여 작업한다.
* N_K와 N_PT의 작업을 구분하여 구현해야 한다.
* kmeans와 num_pt를 모두 합쳐서(ALLREDUCE) 마지막 for-loop를 수행하도록 해여한다.
*/
void update_step(PPOINT kmeans, const PPOINT pt,
		const int istart, const int iend, const int *kindex)
{
	int i, idx;
	int num_pt[N_K], local_num_pt[N_K];
	POINT local_kmeans[N_K];

	for(i=0; i<N_K; i++) {
		kmeans[i].x = kmeans[i].y = 0.0;
		local_kmeans[i].x = local_kmeans[i].y = 0.0;
		num_pt[i] = 0;
		local_num_pt[i] = 0;
	}
	for(i=istart; i<iend; i++) {
		idx = kindex[i];
		local_kmeans[idx].x += pt[i].x;
		local_kmeans[idx].y += pt[i].y;
		local_num_pt[idx]++;
	}

	/**
	* kmeans와 num_pt를 모두 합쳐서(ALLREDUCE) 마지막 for-loop를 수행하도록 해여한다.
	* local_kmeans를 따로 선언하여 작업(send buffer가 필요)
	*/
	MPI_Allreduce(local_kmeans, kmeans, N_K*2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // x,y좌표가 있기 때문에 N_K*2임.
	MPI_Allreduce(local_num_pt, num_pt, N_K, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	for(i=0; i<N_K; i++) {
		kmeans[i].x /= num_pt[i];
		kmeans[i].y /= num_pt[i];
	}
}
