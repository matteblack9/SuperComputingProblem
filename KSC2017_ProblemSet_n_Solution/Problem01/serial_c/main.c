#include "motion_estimation.h"

#define WIDTH		6400
#define HEIGHT		3600
#define SCAN		320
#define FRAMENUM	100

int main()
{
	unsigned char *pImages;
	clock_t t_sta,t_end;
	pImages = (unsigned char *)malloc(sizeof(unsigned char)*WIDTH*HEIGHT*FRAMENUM);
	if (pImages == NULL)
	{
		printf("memory alloc error\n");
		return 0;
	}

	int *nMV;
	nMV = (int*)malloc(sizeof(int) * (2 * FRAMENUM - 2));
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
	for (i = 0; i < FRAMENUM; i++)
	{
		fseek(fp, (long int)i*WIDTH*HEIGHT * 3 / 2, SEEK_SET);
		fread(pImages + (long int)WIDTH*HEIGHT * i, 1, WIDTH*HEIGHT, fp);
	}
	fclose(fp);

	t_end = clock();

	printf("Read   time = %f msec\n", (double)(t_end - t_sta)/1000.0);

	t_sta = clock();
	MotionEstimation(pImages, WIDTH, HEIGHT,SCAN, nMV, FRAMENUM);
	t_end = clock();
	printf("Motion time = %f msec\n", (double)(t_end - t_sta)/1000.0);

	if ((fp = fopen("./result.txt", "w")) != NULL)
	{
		for (i = 0; i < FRAMENUM - 1; i++)
		{
			fprintf(fp, "%5d %5d %5d\n",i,nMV[2 * i + 0], nMV[2 * i + 1]);
		}
		fclose(fp);
	}
//	if ((fp = fopen("./result.txt", "rb")) != NULL)
//	{
//		int nErrorCount = 0;
//		int nMVX, nMVY;
//		for (i = 0; i < FRAMENUM - 1; i++)
//		{
//			fscanf(fp, "%d %d\n", &nMVX, &nMVY);
//			if (nMV[2 * i + 0] != nMVX || nMV[2 * i + 1] != nMVY)
//				nErrorCount++;
//		}
//		printf("num of error: %d\n", nErrorCount);
//		fclose(fp);
//	}
	
	free(pImages);
	free(nMV);

    return 0;
}

