FC=mpif90
FFLAGS=-O3

.PHONY: all

all: fparallel

fparallel: parallel.f90
	$(FC) $(FFLAGS) parallel.f90 -o f_parallel.ex
clean:
	rm -rf *.o *.mod ?_*l.ex
