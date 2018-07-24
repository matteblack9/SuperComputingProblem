#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* File write */
int file_write(char* f_name, int m_size, int d_size, double* matrix, int myrank, int nprocs)
{
	FILE* f_wrt = fopen(f_name, "w");
	int i, j;

	for (i=0; i< nprocs; i++)
	{
		if(myrank == i)
		{
			if(i == 0)
				f_wrt = fopen(f_name, "w");
			else 
				f_wrt = fopen(f_name, "a");

			for(j=0; j < m_size * d_size; j++)
				fprintf(f_wrt, "%f\n", matrix[j]);

			fclose(f_wrt);
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}			      

	return 0;
}

/* Covariance matrix generate */
double *cov_gen(int start_index, int my_n, int n, int myrank, int nprocs)
{
	int i, j, row_offset;
	double *L = (double*)malloc(my_n*n*sizeof(double));
	if (L == NULL)
		exit(EXIT_FAILURE); 

	for (i = 0; i < my_n*n; i++)
		L[i] = 0.0;

	row_offset = start_index - 1;

	for (i = 0; i < my_n; i++)
	{
		L[i * n + i + row_offset] = log(n) * log(n) - cos((i + row_offset) * n);

		if(!(i == 0 && myrank == 0))
			L[i * n + i + row_offset - 1] = sin((i + row_offset) * n);

		if(!(i == my_n - 1 && myrank == nprocs - 1))
			L[i * n + i + row_offset + 1] = sin((i + row_offset + 1) * n);
	}

	return L;
}

/* Multivariate number generate */
double *z_gen(int start_index, int my_n, int d)
{
	int i, j, row_offset;
	double *L = (double*)malloc(my_n*d*sizeof(double));
	if (L == NULL)
		exit(EXIT_FAILURE);

	for (i = 0; i < my_n*d; i++)
		L[i] = 0.0;

	row_offset = start_index - 1;

	for (i = 0; i < my_n; i++)
		for (j = 0; j < d; j++)
			L[i * d + j] = exp( -cos((row_offset+i+1) * j) * sin((row_offset+i+1) * (j+1)) );

	return L;
}

/* Cholesky Decomposition */
double *cholesky(double *A, int start_index, int last_index, int n, int myrank, double *previous_result) 
{
	int i, j, k, row_offset;
	int my_size = last_index;
	double *L = (double*)malloc(last_index*n*sizeof(double));
	double *result = (double*)malloc((last_index-start_index+1)*n*sizeof(double));
	double s;

	if (L == NULL || result == NULL)
		exit(EXIT_FAILURE);

	for (i = 0; i < last_index*n; i++)
		L[i] = 0.0;

	for (i = 0; i < (last_index-start_index+1)*n-1; i++)
		result[i] = 0.0;

	if(myrank != 0)
		memcpy(L, previous_result, (start_index-1)*n*sizeof(double));	

	row_offset = start_index - 1;

	for (i = row_offset; i < last_index; i++)
		for (j = 0; j < i+1; j++) 
		{
			s = 0.0;
			for (k = 0; k < j; k++)
				s += L[i * n + k] * L[j * n + k];
			L[i * n + j] = (i == j) ?
				sqrt(A[(i-row_offset) * n + i] - s) :
				(1.0 / L[j * n + j] * (A[(i-row_offset) * n + j] - s));
		}

	memcpy(result, &L[row_offset*n], (last_index-start_index+1)*n*sizeof(double)); 
	free(L);

	return result;
}

/* Matrix multiplication */
double *mat_mul(double *A, double *B, int n, int m, int myrank, int nprocs, int mysize)
{
	int i, j, k, leftrank, rightrank;
	MPI_Request req[2];
	MPI_Status status[2];
	double *L = (double*)malloc(mysize*m*sizeof(double));
	double *Bl = (double*)malloc(mysize*m*sizeof(double));
	double *Br = (double*)malloc(mysize*m*sizeof(double));

	if (L == NULL || Bl == NULL || Br == NULL)
		exit(EXIT_FAILURE);

	for (i = 0; i < mysize*m; i++)
	{
		L[i] = 0.0;
		Bl[i] = 0.0;
		Br[i] = 0.0;
	}

	leftrank = (myrank-1+nprocs)%nprocs;
	rightrank = (myrank+1)%nprocs;

	if(nprocs > 1)
	{
		MPI_Irecv(Bl, mysize*m, MPI_DOUBLE, rightrank, 1002, MPI_COMM_WORLD, &req[0]);
		MPI_Isend(B, mysize*m, MPI_DOUBLE, leftrank, 1002, MPI_COMM_WORLD, &req[1]);
	}

	for (i = 0; i < mysize; i++)
		for (j = 0; j < m; j++)
			for (k = 0; k < mysize; k++)
				L[i * m + j] = L[i * m +j] + A[i * n + k + mysize*myrank] * B[k * m + j];    

	if(nprocs > 1)
		MPI_Waitall(2, req, status);

	if(nprocs > 2)
	{
		MPI_Irecv(Br, mysize*m, MPI_DOUBLE, leftrank, 1002, MPI_COMM_WORLD, &req[0]);
		MPI_Isend(B, mysize*m, MPI_DOUBLE, rightrank, 1002, MPI_COMM_WORLD, &req[1]);
	}

	if(nprocs > 1)
	{
		for (i = 0; i < mysize; i++)
			for (j = 0; j < m; j++)
				for (k = 0; k < mysize; k++)
					L[i * m + j] = L[i * m +j] + A[i * n + k + mysize*rightrank] * Bl[k * m + j];
	}

	if(nprocs > 2)
		MPI_Waitall(2, req, status);

	if(nprocs > 2)
	{
		for (i = 0; i < mysize; i++)
			for (j = 0; j < m; j++)
				for (k = 0; k < mysize; k++)
					L[i * m + j] = L[i * m +j] + A[i * n + k + mysize*leftrank] * Br[k * m + j];
	}

	free(Bl); 
	free(Br);
	return L;
}

void show_matrix(double *A, int n, int m) 
{
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++)
			printf("%8f ", A[i * m + j]);
		printf("\n");
	}
}

int main(int argc, char **argv) 
{
	int nprocs, myrank;
	int m_size = 3200;
	int d_size = 30;
	int debug_step = 0;
	int start_index, last_index, my_size;
	int i,j,k;
	double *cov_numbers, *z_numbers, *c2, *c3, *temp;
	double elapsed_time = 0.0;
	MPI_Status status;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if(m_size % nprocs != 0)
	{
		if(myrank == 0)
			printf("m_size must be an integer multiple of nprocs!\n");
		MPI_Finalize();
		return 0;
	}

	elapsed_time = -MPI_Wtime();

	start_index = m_size*myrank/nprocs + 1;
	last_index = m_size*(myrank+1)/nprocs;
	my_size = last_index - start_index + 1;

	cov_numbers = cov_gen(start_index, my_size, m_size, myrank, nprocs);

	if(myrank == 0)   
		printf("Step %d completed - generation of a covariance matrix\n", debug_step);

	z_numbers = z_gen(start_index, my_size, d_size);

	if(myrank == 0)
		printf("Step %d completed - generation of a Z matrix\n", debug_step+1);

	if(myrank != 0)
		temp = (double*)malloc(sizeof(double)*(start_index-1)*m_size);  

	if(nprocs > 1)
	{
		if(myrank == 0)
		{               
			c2 = cholesky(cov_numbers, start_index, last_index, m_size, myrank, temp);
			for(int ii = 1; ii < nprocs; ii++)
				MPI_Send(c2, my_size*m_size, MPI_DOUBLE, ii, 1002, MPI_COMM_WORLD);
		}
		else
			MPI_Recv(&temp[0], my_size*m_size, MPI_DOUBLE, 0, 1002, MPI_COMM_WORLD, &status);

		for(int jj = 1; jj < nprocs-1; jj++)
		{
			if(myrank == jj)
			{
				c2 = cholesky(cov_numbers, start_index, last_index, m_size, myrank, temp);
				for(int ii = myrank+1; ii < nprocs; ii++)
					MPI_Send(c2, my_size*m_size, MPI_DOUBLE, ii, 1002, MPI_COMM_WORLD);
			}
			else if(myrank > jj)
				MPI_Recv(&temp[my_size*jj*m_size], my_size*m_size, MPI_DOUBLE, jj, 1002, MPI_COMM_WORLD, &status);

			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

	if(myrank == nprocs-1)
		c2 = cholesky(cov_numbers, start_index, last_index, m_size, myrank, temp);

	if(myrank != 0)
		free(temp);

	if(myrank == 0)
		printf("Step %d completed - Cholesky decomposition\n", debug_step+2);

	/* Step 3 : Correlated Random Number Generator: Matrix Multiplication */
	c3 = mat_mul(c2, z_numbers, m_size, d_size, myrank, nprocs, my_size);

	if(myrank == 0)
		printf("Step %d completed - matrix multiplication\n", debug_step+3);

	/* Step 4 : Write result files*/
	if(myrank == 0)
	{
		system("rm -rf result");
		system("mkdir result");
	}

	file_write("./result/result.txt", my_size, d_size, c3, myrank, nprocs);

	if(myrank == 0)
		printf("Step %d completed - generation of results\n", debug_step+4);

	free(cov_numbers);
	free(z_numbers);
	free(c2);
	free(c3);

	if(myrank == 0)
		printf("Step %d completed - memory release\n", debug_step+5);

	elapsed_time += MPI_Wtime();

	if(myrank == 0)
		printf("Total elapsed time = %g (secs)\n", elapsed_time);

	MPI_Finalize();
	return 0;
}





