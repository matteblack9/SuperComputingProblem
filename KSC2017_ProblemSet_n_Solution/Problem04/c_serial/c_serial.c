#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* File write */
int file_write(char* f_name, int m_size, int d_size, double* matrix)
{
	FILE* f_wrt = fopen(f_name, "w");
	int i = 0;

	for(i=0; i < m_size * d_size; i++)
		fprintf(f_wrt, "%f\n", matrix[i]);

	fclose(f_wrt);
	return 0;
}

/* Covariance matrix generate */
double *cov_gen(int n)
{
	int i, j;
	double *L = (double*)malloc(n*n*sizeof(double));
	if (L == NULL)
		exit(EXIT_FAILURE); 

	for (i = 0; i < n*n; i++)
		L[i] = 0.0;

	for (i = 0; i < n; i++)
		L[i * n + i ] = log(n) * log(n) - cos(i * n);

	for (j = 1; j < n; j++)
	{
		L[(j - 1) * n + j ] = sin(j * n);
		L[ j * n + j - 1] = L[(j - 1) * n + j ];
	}

	return L;
}


/* Multivariate number generate */
double *z_gen(int n, int d)
{
	int i, j;
	double *L = (double*)malloc(n*d*sizeof(double));
	if (L == NULL)
		exit(EXIT_FAILURE); 

	for (i = 0; i< n*d; i++)
		L[i] = 0.0;

	for (i = 0; i < n; i++)
		for (j = 0; j < d; j++)
			L[i * d + j] = exp( - cos((i+1) * j) * sin((i+1)* (j+1)) );

	return L;
}


/* Cholesky Decomposition */
double *cholesky(double *A, int n) 
{
	int i, j, k;
	double s;
	double *L = (double*)malloc(n*n*sizeof(double));
	if (L == NULL)
		exit(EXIT_FAILURE);

	for (i = 0; i < n*n; i++)
		L[i] = 0.0;

	for (i = 0; i < n; i++)
		for (j = 0; j < (i+1); j++) 
		{
			s = 0.0;
			for (k = 0; k < j; k++)
				s += L[i * n + k] * L[j * n + k];
			L[i * n + j] = (i == j) ?
				sqrt(A[i * n + i] - s) :
				(1.0 / L[j * n + j] * (A[i * n + j] - s));
		}

	return L;
}

/* Matrix multiplication */
double *mat_mul(double *A, double *B, int n, int m)
{
	int i, j, k;
	double *L = (double*)malloc(n*m*sizeof(double));
	if (L == NULL)
		exit(EXIT_FAILURE);

	for (i = 0; i < n*m; i++)
		L[i] = 0.0;

	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			for (k = 0; k < n; k++)
				L[i * m + j] = L[i * m +j] + A[i * n + k] * B[k * m + j];

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

int main() 
{
	/* Step 1 : Generate covariance matrix & standard normal random number files */
	int m_size = 3200;
	int d_size = 30;
	int debug_step = 0;

	double *cov_numbers = cov_gen(m_size);
	printf("Step %d completed - generation of a covariance matrix\n", debug_step);

	double *z_numbers = z_gen(m_size, d_size);

	printf("Step %d completed - generation of a Z matrix\n", debug_step+1);

	/* Step 2 : Cholesky Decomposition of Sigma-matrix */
	double *c2 = cholesky(cov_numbers, m_size);

	printf("Step %d completed - Cholesky decomposition\n", debug_step+2);

	/* Step 3 : Correlated Random Number Generator: Matrix Multiplication */
	double *c3 = mat_mul(c2, z_numbers, m_size, d_size);

	printf("Step %d completed - matrix multiplication\n", debug_step+3);

	/* Step 4 : Write result files*/
	system("rm -rf result");
	system("mkdir result");

	file_write("./result/result.txt", m_size, d_size, c3);
	printf("Step %d completed - generation of results\n", debug_step+4);

	free(cov_numbers);
	free(z_numbers);
	free(c2);
	free(c3);

	printf("Step %d completed - memory release\n", debug_step+5);

	return 0;
}





