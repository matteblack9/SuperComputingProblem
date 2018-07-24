#include "cgsolver.h"

void cgsolver(int size, int myrank, int ncpus, double *matrix, double *rhs, double *solution, int maxiteration, double tolerance)
{
	int ii, jj, kk, rank[3];
	int *load, *sindex, maxsize;
	double alpha=0.0, beta=0.0, temp1, temp2, res0tol=0.0;
	double *res, *p, *Ax, *Ap, *xtl, *xtr;

	load = (int*) malloc(ncpus*sizeof(int));
	sindex = (int*) malloc(ncpus*sizeof(int));

	MPI_Allgather(&size, 1, MPI_INT, load, 1, MPI_INT, MPI_COMM_WORLD);

	res = (double*) malloc(size*sizeof(double));
	p   = (double*) malloc(size*sizeof(double));
	Ax  = (double*) malloc(size*sizeof(double));
	Ap  = (double*) malloc(size*sizeof(double));

	rank[0] = (myrank-1+ncpus)%ncpus;
	rank[1] = myrank;
	rank[2] = (myrank+1)%ncpus;

	maxsize = 0;
	sindex[0] = 0;
	for(ii=0; ii<ncpus; ii++)
	{
		if(maxsize<load[ii])
			maxsize=load[ii];
                if(ii!=0)
                        sindex[ii]=sindex[ii-1]+load[ii-1];
	}

	xtl = (double*) malloc(maxsize*sizeof(double));
	xtr = (double*) malloc(maxsize*sizeof(double));

	multiply(rank, load, sindex, ncpus, matrix, solution, Ax, xtl, xtr);

	for (ii=0; ii<size; ii++)
	{
		res[ii] = rhs[ii]-Ax[ii];
		p[ii]   = res[ii];
	}

	res0tol = innerproduct(res, res, size);


	if(myrank == 0)
		printf("[CG] Conjugate gradient is started.\n");

	for (ii=0; ii<maxiteration; ii++)
	{
		if ((myrank==0)&&(ii%20==0)&&(ii!=0))
			printf("mse = %e with criteria of %e at %5d iterations.\n", sqrt(temp2/res0tol), tolerance, ii);

		temp1 = innerproduct(res, res, size);
		multiply(rank, load, sindex, ncpus, matrix, p, Ap, xtl, xtr);
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

	if(myrank==0)
		printf("[CG] Finished with total iteration = %d, mse = %e.\n", ii+1, sqrt(temp2/res0tol));

	free(xtl);
	free(xtr);
	free(res);
	free(p);
	free(Ax);
	free(Ap);
	free(sindex);
	free(load);
}

double innerproduct(double *x, double *y, int size)
{
	int ii;
	double result, globalresult;

	result = 0.0;
	for(ii=0; ii<size; ii++)
		result += x[ii]*y[ii];

	MPI_Allreduce(&result, &globalresult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return globalresult;
}

void multiply(int *rank, int *load, int *sindex, int ncpus, double *matrix, double *x, double *y, double *xtl, double *xtr)
{
        int ii, jj, leftrank, myrank, rightrank, globalmatrixDOF;
        MPI_Request req[2];
        MPI_Status status[2];

	leftrank  = rank[0];
	myrank    = rank[1];
	rightrank = rank[2];
	globalmatrixDOF = sindex[ncpus-1]+load[ncpus-1];

        for(ii=0; ii<load[myrank]; ii++)
                y[ii] = 0.0;

	if(ncpus > 1)
	{
		MPI_Irecv(xtl, load[rightrank], MPI_DOUBLE, rightrank, 1002, MPI_COMM_WORLD, &req[0]);
		MPI_Isend(x, load[myrank], MPI_DOUBLE, leftrank, 1002, MPI_COMM_WORLD, &req[1]);
	}

	for (ii=0; ii<load[myrank]; ii++)
		for (jj=0; jj<load[myrank]; jj++)
			y[ii] += matrix[ii*globalmatrixDOF+sindex[myrank]+jj]*x[jj];
	if(ncpus > 1)
		MPI_Waitall(2, req, status);

        if(ncpus > 2)
        {
                MPI_Irecv(xtr, load[leftrank], MPI_DOUBLE, leftrank, 1002, MPI_COMM_WORLD, &req[0]);
                MPI_Isend(x, load[myrank], MPI_DOUBLE, rightrank, 1002, MPI_COMM_WORLD, &req[1]);
        }

	if(ncpus > 1)
	{
        	for (ii=0; ii<load[myrank]; ii++)
                	for (jj=0; jj<load[rightrank]; jj++)
                        	y[ii] += matrix[ii*globalmatrixDOF+sindex[rightrank]+jj]*xtl[jj];
	}

	if(ncpus > 2)
		MPI_Waitall(2, req, status);	

	if(ncpus > 2)
	{
        	for (ii=0; ii<load[myrank]; ii++)
                	for (jj=0; jj<load[leftrank]; jj++)
				y[ii] += matrix[ii*globalmatrixDOF+sindex[leftrank]+jj]*xtr[jj];
	}	
}

