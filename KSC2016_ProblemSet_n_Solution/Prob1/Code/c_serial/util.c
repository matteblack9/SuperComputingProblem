#include "kmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double distance(const POINT pt1, const POINT pt2)
{   
    return( sqrt((pt2.x-pt1.x)*(pt2.x-pt1.x)+(pt2.y-pt1.y)*(pt2.y-pt1.y)) );
}   

int check_diff(const PPOINT kmeans, const PPOINT kmeans_old)
{   
    for(int i=0; i<N_K; i++)
        if( distance(kmeans[i], kmeans_old[i]) > 0 ) return 1;
    return 0;
}   

void initialize_data(PPOINT *kmeans, PPOINT *pt,
		PPOINT *kmeans_old, int **kindex)
{
	*kmeans = (PPOINT) malloc(sizeof(POINT) * N_K);
	*kmeans_old = (PPOINT) malloc(sizeof(POINT) * N_K);
	*pt = (PPOINT) malloc(sizeof(POINT) * N_PT);
	*kindex = (int *) malloc(sizeof(int) * N_PT);

	FILE *fp = fopen("input.dat", "rb");
	fread(*pt, sizeof(POINT), N_PT, fp);
	fread(*kmeans, sizeof(POINT), N_K, fp);
	fclose(fp);
}

void check_result(const PPOINT kmeans, PPOINT kmeans_old)
{
	FILE *fp = fopen("result.dat", "rb");
	fread(kmeans_old, sizeof(POINT), N_K, fp);
	fclose(fp);	

	double err = 0.0;
	for(int i=0; i<N_K; i++)
		err += distance(kmeans[i], kmeans_old[i]);

	if( err < 1.0E-8 )
//	if( check_diff(kmeans, kmeans_old) )
		printf("Result: Pass!\n");
	else
		printf("Result: Fail!\n");
}

void release_data(PPOINT kmeans, PPOINT pt, PPOINT kmeans_old, int *kindex)
{
	free(kmeans);
	free(pt);
	free(kmeans_old);
	free(kindex);
}
