#include <stdio.h>
#include <omp.h>

int main()
{
	int a[10], tid, i;
	omp_set_num_threads(4);

#pragma omp parallel shared(a) private(tid)
{
	tid = omp_get_thread_num();
	a[0] = tid + 1;
}

	for(i = 0; i < 4; i++)
	printf("a[%d] = %d\n", i ,a[i]);
	return 0;

}