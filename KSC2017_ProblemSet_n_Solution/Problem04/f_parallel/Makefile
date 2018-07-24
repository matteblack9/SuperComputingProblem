MPIFC = mpif90 
FFLAGS = -cpp -O3  
MPI_SRCS = f_parallel.f90

all: mpi

mpi:
	$(MPIFC) $(FFLAGS) $(MPI_SRCS) -o f_parallel.ex

clean :
	rm -rf f_parallel.ex f_parallel.o result *.out

