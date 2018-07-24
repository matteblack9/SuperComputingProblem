#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"
#define GRANULARITY 5		// it should be less than nx
#define MASTER 0



int nprocs, myrank;


struct job {
	bool working;			// true:working, false:quit
	int site[GRANULARITY];
};


struct job_msg {
	int solutions_found;
	int origin;
};



double time_diff(struct timeval t1, struct timeval t2) {
	double sec1, sec2;

	sec1 = t1.tv_sec + t1.tv_usec*1e-6;
	sec2 = t2.tv_sec + t2.tv_usec*1e-6;

	return sec2-sec1;
}



int master_build_house(int nx, int column, int site[]) {
   	int num_sol = 0;
	bool is_sol;
	int i, j;

	for (i=0; i<nx; i++) {
   		// Try to build a house in each line of column
	   	site[column] = i;

	   	// Check if this placement is still a solution
		is_sol = true;
		for (j=column-1; j>=0; j--) {
		   	if ((site[column] == site[j])               ||
				(site[column] == site[j] - (column-j))  ||
				(site[column] == site[j] + (column-j))) {
			   	is_sol = false;
			   	break;
		   	}
	   	}

	   	if (is_sol) {
		  	if (column == GRANULARITY-1) {
				// If this is the last level (granularity of the job),
				// send a next job to a worker
				num_sol += send_job_worker(nx, site);
		   	} 
			else {
			   	// The placement is not complete.
				// Try to place the house on the next column
				num_sol += master_build_house(nx, column+1, site);
		   	}
	   	}
   	}

   	return num_sol;
}



int worker_build_house(int nx, int column, int sub_site[]) {
   	int num_sol = 0;
	bool is_sol;
	int i, j;

   	// Try to build a house in each line of column
	for (i=0; i<nx; i++) {
	   	sub_site[column] = i;

	   	// Check if this placement is still a solution
		is_sol = true;
		for (j=column-1; j>=0; j--) {
		   	if ((sub_site[column] == sub_site[j])               ||
				(sub_site[column] == sub_site[j] - (column-j))  ||
				(sub_site[column] == sub_site[j] + (column-j))) {
			   	is_sol = false;
			   	break;
		   	}
	   	}

	   	if (is_sol) { 
		  	if (column == nx-1) {
			   	// If this is the last column, printout the solution
				num_sol += 1;

		   	} else {
			   	// The placement is not complete.
				// Try to place the house on the next column
				num_sol += worker_build_house(nx, column+1, sub_site);
		   	}
	   	}
   	}

   	return num_sol;
}



int send_job_worker(int nx, int site[]) {
	int num_sol = 0;		// The number of solutions found meanwhile
	int i;
	struct job todo;
	struct job_msg msg;
	
	// Set the job
	todo.working = true;

	for (i=0; i<GRANULARITY; i++)
		todo.site[i] = site[i];

	// Recieve the last result from a worker
	MPI_Recv(&msg, sizeof(msg), MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, 
			 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	num_sol = msg.solutions_found;

	// Send the new job to the worker
	MPI_Send(&todo, sizeof(todo), MPI_BYTE, msg.origin, 0, MPI_COMM_WORLD);

	return num_sol;
}



int wait_remaining_results() {
	// Wait for remaining results, sending a quit whenever a new result arrives
	
	int num_sol = 0;
	int n_workers = nprocs-1;
	struct job todo;
	struct job_msg msg;

	todo.working = false;

	while (n_workers > 0) {
		// Receive a message from a worker
		MPI_Recv(&msg, sizeof(msg), MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, 
				 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		num_sol += msg.solutions_found;

		MPI_Send(&todo, sizeof(todo), MPI_BYTE, msg.origin, 0, MPI_COMM_WORLD);

		n_workers -= 1;
	}

	return num_sol;
	
}



void worker(int nx) {
	// There is a default message which lets a worker request a 
	// job reporting the number of solutions found in the last iteration
	
	int num_sol;
	struct job_msg msg;

    msg.origin = myrank;
	msg.solutions_found = 0;

	// Request initial job
	MPI_Send(&msg, sizeof(msg), MPI_BYTE, MASTER, 0, MPI_COMM_WORLD);

	while (true) {
		// Wait for a job or a quit message
	    struct job todo;
		MPI_Recv(&todo, sizeof(todo), MPI_BYTE, MPI_ANY_SOURCE, 
				 MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if (todo.working == false)
		    break;

		num_sol = worker_build_house(nx, GRANULARITY, todo.site);

		// Ask for more work
		msg.origin = myrank;
		msg.solutions_found = num_sol;
		MPI_Send(&msg, sizeof(msg), MPI_BYTE, MASTER, 0, MPI_COMM_WORLD);
	}
}



void main(int argc, char* argv[]) {
	int num_sol;
	int nx=14;
	int *site;
	struct timeval start, end;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


	if (myrank == MASTER) {
		printf("nx=%d\n", nx);
		printf("nprocs=%d\n", nprocs);
		site = (int*)malloc(nx*sizeof(int));

		gettimeofday(&start, NULL);
		num_sol = master_build_house(nx, 0, site);
		num_sol += wait_remaining_results();
		gettimeofday(&end, NULL);

		printf("num_sol=%d\n", num_sol);
		printf("elapsed time=%.3lf sec\n", time_diff(start,end));
	}
	else {
		worker(nx);
	}

	MPI_Finalize();
}
