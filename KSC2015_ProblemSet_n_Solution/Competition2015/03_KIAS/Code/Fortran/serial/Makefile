FC=gfortran
FFLAGS=-O3

.PHONY: all

all: fserial

fserial: sequential.f90
	$(FC) $(FFLAGS) sequential.f90 -o f_serial.ex
clean:
	rm -rf *.o *.mod ?_*l.ex
