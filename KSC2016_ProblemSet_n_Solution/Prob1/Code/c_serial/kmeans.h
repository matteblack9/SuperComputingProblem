#ifndef _KMEANS_H_
#define _KMEANS_H_

#define N_K 16
#define N_PT 20000000

typedef struct {
	double x, y;
} POINT, *PPOINT;

double distance(const POINT p1, const POINT p2);
int check_diff(const PPOINT kmeans, const PPOINT kmeans_old);

void initialize_data(PPOINT *kmeans, PPOINT *pt, PPOINT *kmeans_old, int **kindex);
void check_result(const PPOINT kmeans, PPOINT kmeans_old);
void release_data(PPOINT kmeans, PPOINT pt, PPOINT kmeans_old, int *kindex);

void assignment_step(const PPOINT kmeans, const PPOINT pt, int *kindex);
void update_step(PPOINT kmeans, const PPOINT pt, const int *kindex);

#endif  //  _KMEANS_H_
