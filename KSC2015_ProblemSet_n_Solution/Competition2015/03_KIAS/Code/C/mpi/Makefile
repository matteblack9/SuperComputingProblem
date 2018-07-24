CC=mpicc
CFLAGS=-O3 

.PHONY: all

all: cparallel

cparallel: parallel.c
	$(CC) $(CFLAGS) parallel.c -lm -std=c99 -o c_parallel.ex
clean:
	rm -rf *.o *.mod ?_*l.ex
