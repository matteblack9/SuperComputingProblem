#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#define num_steps 1000000000
#define NUM_THREADS 4

int main()
{
	double sum[4],a, step, x, pi;
	double t1, t2;
	clock_t t;
	int i,id;

	step = 1./(double)num_steps;

	t = clock();
#pragma omp parallel private(id,x)
{
	id = omp_get_thread_num();
	sum[id] = 0.0;
	for(i = omp_get_thread_num(); i < num_steps; i = i + NUM_THREADS)
	{
		x = (i + 0.5)*step;
		sum[id] += 4.0/(1.0 + x*x);
	}
}
	t = clock() - t;

	printf("%d ms\n", t);
	printf("\n");
	for(i = 0, pi = 0.0; i < NUM_THREADS; i++)
		pi += step * sum[i];
	printf(" Sum = %d.15\n", sum[0] + sum[1] + sum[2] + sum[3]);
	printf(" numerical pi = %.15f \n", pi);
	printf(" analytical pi = %.15f\n", acos(-1.0));
	printf(" Error = %E \n", fabs(acos(-1.0)-pi));

	return 0;
}