#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kmeans.h"

void assignment_step(const PPOINT kmeans, const PPOINT pt, int *kindex)
{
	double d1, d2;

	for(int i=0; i<N_PT; i++) {
		kindex[i] = 0;
		d1 = distance(pt[i], kmeans[0]);
		for(int j=1; j<N_K; j++) {
			d2 = distance(pt[i], kmeans[j]);
			if( d1 > d2 ) {
				d1 = d2;
				kindex[i] = j;
			}
		}
	}
}

void update_step(PPOINT kmeans, const PPOINT pt, const int *kindex)
{
	int i, idx;
	int num_pt[N_K];

	for(i=0; i<N_K; i++) {
		kmeans[i].x = kmeans[i].y = 0.0;
		num_pt[i] = 0;
	}
	for(i=0; i<N_PT; i++) {
		idx = kindex[i];
		kmeans[idx].x += pt[i].x;
		kmeans[idx].y += pt[i].y;
		num_pt[idx]++;
	}

	for(i=0; i<N_K; i++) {
		kmeans[i].x /= num_pt[i];
		kmeans[i].y /= num_pt[i];
	}
}
