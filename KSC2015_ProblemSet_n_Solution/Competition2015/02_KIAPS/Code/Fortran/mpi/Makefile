FC=mpif90
FFLAGS=-O3

.PHONY: all

all: fparallel

fparallel: sub_mod1.o sub_mod2par.o 1d_parallel.f90
	$(FC) $(FFLAGS) 1d_parallel.f90 sub_mod1.o sub_mod2par.o -o f_parallel.ex
sub_mod2par.o: sub_mod2par.f90
	$(FC) -c $(FFLAGS) sub_mod2par.f90
sub_mod1.o: sub_mod1.f90
	$(FC) -c $(FFLAGS) sub_mod1.f90
clean:
	rm -rf *.o *.mod ?_*l.ex
