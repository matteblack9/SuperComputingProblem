FC=mpif90
FFLAGS=-O3

.PHONY: all

all: fparallel

fparallel: b.f90
		$(FC) $(FFLAGS) b.f90 -o f_parallel.ex
clean:
		rm -rf *.o *.mod ?_*l.ex
