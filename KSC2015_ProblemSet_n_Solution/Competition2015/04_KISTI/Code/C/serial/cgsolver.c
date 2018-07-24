#include "cgsolver.h"

void cgsolver(int size, double *matrix, double *rhs, double *solution, int maxiteration, double tolerance)
{
	int ii, jj, kk;
	double alpha=0.0, beta=0.0, temp1, temp2, res0tol=0.0;
	double *res, *p, *Ax, *Ap;

	res = (double*) malloc(size*sizeof(double));
	p   = (double*) malloc(size*sizeof(double));
	Ax  = (double*) malloc(size*sizeof(double));
	Ap  = (double*) malloc(size*sizeof(double));

	multiply(size, matrix, solution, Ax);

	for (ii=0; ii<size; ii++)
	{
		res[ii] = rhs[ii]-Ax[ii];
		p[ii]   = res[ii];
	}

	res0tol = innerproduct(res, res, size);

	printf("[CG] Conjugate gradient is started.\n");

	for (ii=0; ii<maxiteration; ii++)
	{
		if ((ii%20==0)&&(ii!=0))
			printf("[CG] mse %e with a tolerance criteria of %e at %5d iterations.\n", sqrt(temp2/res0tol), tolerance, ii);

		temp1 = innerproduct(res, res, size);
		multiply(size, matrix, p, Ap);
		temp2 = innerproduct(Ap, p, size);

		alpha=temp1/temp2;

		for (jj=0; jj<size; jj++)
		{
			solution[jj] = solution[jj] + alpha*p[jj];
			res[jj] = res[jj] - alpha*Ap[jj];
		}

		temp2 = innerproduct(res, res, size);
		
		if (sqrt(temp2/res0tol) < tolerance)
			break;

		beta = temp2/temp1;

		for (jj=0; jj<size; jj++)
			p[jj]= res[jj] + beta*p[jj];

	}

	printf("[CG] Finished with total iteration = %d, mse = %e.\n", ii, sqrt(temp2/res0tol));

	free(res);
	free(p);
	free(Ax);
	free(Ap);
}

double innerproduct(double *x, double *y, int size)
{
	int ii;
	double result;

	result = 0.0;

	for(ii=0; ii<size; ii++)
		result += x[ii]*y[ii];

	return result;
}

void multiply(int size, double *matrix, double *x, double *y)
{
	int ii, jj;

        for (ii=0; ii<size; ii++)       // initialize y
                y[ii]=0.0;

        for (ii=0; ii<size; ii++)
                for (jj=0; jj<size; jj++)
                        y[ii] += matrix[ii*size+jj] * x[jj];
}

