#include <mpi.h>
#define m 6
#define n 9

void para_range(int n1, int n2, int nprocs, int myrank, int *ista, int *iend){
    int iwork1, iwork2;
    iwork1 = (n2-n1+1)/nprocs;
    iwork2 = (n2-n1+1)%nprocs;
    *ista = myrank*iwork1 + n1 + min(myrank, iwork2);
    *iend = *ista + iwork1 - 1;
    if(iwork2 > myrank) *iend = *iend + 1;
}

int min(int a, int b){
	if(a >= b) return b;
	else return a;
}

int main(int argc, char *argv[])
{
	int i, j , nprocs, myrank;
	double a[m][n], b[m][n];
	double works1[m], workr1[m], works2[m], workr2[m];
	int ista, iend, ista2, iend1, jnext, jprev;
	MPI_Request jsend1, jsend2, jrecv1, jrecv2;
	MPI_Status jstatus;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	para_range(0, n-1, nprocs, myrank, &ista, &iend);
	ista2 = ista; iend1 = iend;
	if(myrank == 0) ista2 = 1;
	if(myrank == nprocs - 1) iend1 = n - 2;
	jnext = myrank + 1;
	jprev = myrank - 1;
	for(j = 0; j < n; j++)
		for(i = ista; i <= iend; i++) a[i][j] = (i + 1) + 10.0 * j;
	if(myrank != nprocs - 1)
		for(j = 0; j < n; j++) works2[j] = a[ista][j];
	if(myrank != 0)
		for(j = 0; j < n; j++) works1[j] = a[iend][j];
	MPI_Isend(works1, n, MPI_DOUBLE, iprev,  1, MPI_COMM_WORLD, &jsend1);
	MPI_Isend(works2, n, MPI_DOUBLE, inext,  1, MPI_COMM_WORLD, &jsend2);
	MPI_Irecv(works1, n, MPI_DOUBLE, inext,  1, MPI_COMM_WORLD, &jsend3);
	MPI_Irecv(works2, n, MPI_DOUBLE, iprev,  1, MPI_COMM_WORLD, &jsend4);

	MPI_Wait(&jsend1, &jstatus)
	MPI_Wait(&jsend2, &jstatus)
	MPI_Wait(&jsend3, &jstatus)
	MPI_Wait(&jsend4, &jstatus)

	for(i = ista2; i < iend1; i++)
		for(j = 1; j < n-1; j++)
			b[i][j] = a[i-1][j]
	MPI_Finalize();
}