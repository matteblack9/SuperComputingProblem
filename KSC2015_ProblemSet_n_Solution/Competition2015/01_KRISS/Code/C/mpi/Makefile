CC=mpicc
CFLAGS=-O3 

.PHONY: all

all: cparallel
cparallel: b.c
	$(CC) $(CFLAGS) b.c -lm -std=c99 -o c_parallel.ex
clean:
	rm -rf *.o *.mod ?_*l.ex
