#include "motion_estimation.h"

#define WIDTH		6400
#define HEIGHT		3600
#define SCAN		320
#define FRAMENUM	100

void para_range(int n1, int n2, int nprocs, int myrank, int *ista, int *iend, int *pload, int *pdisp)
{
	int iwork1, iwork2, i;
	iwork1 = (n2-n1+1)/nprocs;
	iwork2 = (n2-n1+1) % nprocs;
	*ista= myrank*iwork1 + n1 + MIN(myrank, iwork2);
	*iend = *ista + iwork1 -1;
	if(iwork2>myrank) *iend = *iend +1;
	for ( i = 0 ; i < nprocs ; i++)
	{
		pload[i] = iwork1;
		if( i < iwork2) pload[i]++;
	}
	pdisp[0] = 0;
	for ( i = 1 ; i < nprocs ; i++)
	{
		pdisp[i] = pdisp[i-1] + pload[i-1];
	}
}

int main(int argc, char *argv[])
{
	unsigned char *pImages;
	int nprocs,myrank,ista,iend,isize,iendp;
	int *pload,*pdisp;
	clock_t t_sta,t_end;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

	pload = (int*)malloc(sizeof(int) * nprocs);
	pdisp = (int*)malloc(sizeof(int) * nprocs);
	para_range(0,FRAMENUM-2,nprocs,myrank,&ista,&iend,pload,pdisp);

	isize = iend - ista + 1;

	printf("nprocs = %d, myrank = %d, ista = %d, iend = %d, isize = %d, pload = %d, pdisp =  %d\n",nprocs,myrank,ista,iend,isize,pload[myrank],pdisp[myrank]);

	pImages = (unsigned char *)malloc(sizeof(unsigned char)*WIDTH*HEIGHT*(isize+1));
	if (pImages == NULL)
	{
		printf("memory alloc error\n");
		return 0;
	}

	int *nMV, *nMV_sum;
	nMV     = (int*)malloc(sizeof(int) * (2 * isize));
	nMV_sum = (int*)malloc(sizeof(int) * (2 * FRAMENUM - 2));
	if (nMV == NULL)
	{
		free(pImages);
		printf("memory alloc error\n");
		return 0;
	}

	t_sta = clock();
	FILE *fp;

	if ((fp = fopen("../video_6400_3600.yuv", "rb")) == NULL)
	{
		printf("file open error\n");
		return 0;
	}

	int i;
	for (i = 0; i < isize + 1; i++)
	{
		fseek(fp, (long int)(ista+i)*WIDTH*HEIGHT * 3 / 2, SEEK_SET);
		fread(pImages + (long int)WIDTH*HEIGHT * i, 1, WIDTH*HEIGHT, fp);
	}
	fclose(fp);

	t_end = clock();

	printf("Read   time = %f msec\n", (double)(t_end - t_sta)/1000.0);

	t_sta = clock();
	MotionEstimation(pImages, WIDTH, HEIGHT,SCAN, nMV, 0, isize+1, ista);
	t_end = clock();
	printf("Motion time = %f msec\n", (double)(t_end - t_sta)/1000.0);

	for (i = 0; i < 2 * FRAMENUM - 2; i++)
	{
		nMV_sum[i] = 0;
	}

	for (i = 0; i < nprocs; i++)
	{
		pload[i] = pload[i] * 2;
		pdisp[i] = pdisp[i] * 2;
	}

	MPI_Gatherv(nMV, 2*isize, MPI_INT, nMV_sum, pload, pdisp, MPI_INT, 0, MPI_COMM_WORLD);

	if(myrank == 0)
	{
		if ((fp = fopen("./result.txt", "w")) != NULL)
		{
			for (i = 0; i < FRAMENUM-1; i++)
			{
				fprintf(fp, "%5d %5d %5d\n",i,nMV_sum[2 * i + 0], nMV_sum[2 * i + 1]);
			}
			fclose(fp);
		}
	}
	
	free(pImages);
	free(nMV);
	free(nMV_sum);
	free(pload);
	free(pdisp);
	MPI_Finalize();

    return 0;
}

