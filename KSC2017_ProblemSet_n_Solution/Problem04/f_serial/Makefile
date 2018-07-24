MPIFC = mpif90 
FFLAGS = -O3
SERIAL_SRCS = f_serial.f90
MPI_SRCS = f_mpi.f90

all :
	$(MPIFC) $(FFLAGS) $(SERIAL_SRCS) -o f_serial.ex

clean :
	rm -rf f_serial.ex result f_serial.o *.out


