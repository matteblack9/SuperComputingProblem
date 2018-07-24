FC=ifort
FFLAGS=-O3

.PHONY: all

all: fserial

fserial: b_serial.f90
		$(FC) $(FFLAGS) b_serial.f90 -o f_serial.ex
clean:
		rm -rf *.o *.mod ?_*l.ex
