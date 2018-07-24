#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "kmeans.h"

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

	MPI_Allreduce(local_kmeans, kmeans, N_K*2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(local_num_pt, num_pt, N_K, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	for(i=0; i<N_K; i++) {
		kmeans[i].x /= num_pt[i];
		kmeans[i].y /= num_pt[i];
	}
}
