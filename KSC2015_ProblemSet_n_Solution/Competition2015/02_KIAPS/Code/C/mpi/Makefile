CC=mpicc
CFLAGS=-O3 

.PHONY: all

all: cparallel

cparallel: 1d_parallel.c
	$(CC) $(CFLAGS) 1d_parallel.c -std=c99 -o c_parallel.ex
clean:
	rm -rf *.o *.mod ?_*l.ex
