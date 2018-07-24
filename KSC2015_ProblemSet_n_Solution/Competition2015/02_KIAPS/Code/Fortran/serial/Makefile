FC=gfortran
FFLAGS=-O3

.PHONY: all

all: fserial

fserial: sub_mod1.o sub_mod2.o 1d_serial.f90
	$(FC) $(FFLAGS) 1d_serial.f90 sub_mod1.o sub_mod2.o -o f_serial.ex
sub_mod2.o: sub_mod2.f90
	$(FC) -c $(FFLAGS) sub_mod2.f90
sub_mod1.o: sub_mod1.f90
	$(FC) -c $(FFLAGS) sub_mod1.f90
clean:
	rm -rf *.o *.mod ?_*l.ex
