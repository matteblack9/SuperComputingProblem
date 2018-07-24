#include <stdio.h>
#include <math.h>
#include <time.h>
#define num_steps 1000000000

int main()
{
	double sum, step, x, pi;
	double t1, t2;
	clock_t t;
	int i;

	sum = 0.0;
	step = 1./(double)num_steps;

	t = clock();
	for(i = 1; i < num_steps; i++)
	{
		x = (i - 0.5)*step;
		sum = sum + 4.0/(1.0 + x*x);
	}
	t = clock() - t;

	printf("%d ms", t);
	printf("\n");
	pi = step * sum;
	printf(" Sum = %.15f", sum);
	printf(" numerical pi = %.15f \n", pi);
	printf(" analytical pi = %.15f\n", acos(-1.0));
	printf(" Error = %E \n", fabs(acos(-1.0)-pi));

	return 0;
}