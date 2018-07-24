#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <sys/time.h>



double time_diff(struct timeval t1, struct timeval t2) {
	double sec1, sec2;

	sec1 = t1.tv_sec + t1.tv_usec*1e-6;
	sec2 = t2.tv_sec + t2.tv_usec*1e-6;

	return sec2-sec1;
}




int build_house(int nx, int column, int site[]) {
   	int num_sol = 0;
	bool is_sol;
	int i, j;

   	// Try to build a house in each line of column
	for (i=0; i<nx; i++) {
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
		  	if (column == nx-1) {
			   	// If this is the last column, printout the solution
				num_sol += 1;

		   	} else {
			   	// The placement is not complete. 
				// Try to place the house on the next column
				num_sol += build_house(nx, column+1, site);
		   	}
	   	}
   	}

   	return num_sol;
}



void main() {
	int num_sol;
	int nx=14;
	int *site;
	struct timeval start, end;

	printf("nx=%d\n", nx);
	site = (int*)malloc(nx*sizeof(int));

	gettimeofday(&start, NULL);
   	num_sol = build_house(nx, 0, site);
	gettimeofday(&end, NULL);

	printf("num_sol=%d\n", num_sol);
	printf("elapsed time=%.6lf sec\n", time_diff(start,end));
}
