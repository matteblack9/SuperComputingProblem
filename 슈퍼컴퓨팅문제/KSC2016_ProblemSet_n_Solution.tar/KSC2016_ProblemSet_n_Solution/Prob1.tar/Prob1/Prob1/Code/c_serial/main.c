#include <stdio.h>
#include <mpi.h>
#include "kmeans.h"

int main(int argc, char *argv[])
{
	PPOINT kmeans, pt, kmeans_old;
	int *kindex;
	int iter, nprocs, rank;
	double time1, time2;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

///////////////////////////////////////////////////////////////
//  1. read point data and centers
///////////////////////////////////////////////////////////////
	initialize_data(&kmeans, &pt, &kmeans_old, &kindex);

///////////////////////////////////////////////////////////////
//  2. execute k-means clustering
///////////////////////////////////////////////////////////////
	time1 = MPI_Wtime();

	for(iter=0; iter<100; iter++) {
		assignment_step(kmeans, pt, kindex);
		for(int i=0; i<N_K; i++) {
			kmeans_old[i].x = kmeans[i].x;
			kmeans_old[i].y = kmeans[i].y;
		}
		update_step(kmeans, pt, kindex);
		if( check_diff(kmeans, kmeans_old) == 0 )
			break;
	}

	time2 = MPI_Wtime();

///////////////////////////////////////////////////////////////
//  3. report and verify results
///////////////////////////////////////////////////////////////
	if( rank == 0 ) {
		printf("Iteration: %d\tExecution Time: %lf\n", iter, time2-time1);
		check_result(kmeans, kmeans_old);
	}
	release_data(kmeans, pt, kmeans_old, kindex);

	MPI_Finalize();

	return 0;
}
