CC = mpicc
CFLAGS = -O0 -W -Wall
SERIAL_SRCS = Boyer-Moore.c Main_serial.c
MPI_SRCS = Boyer-Moore.c Main_mpi.c

all: serial mpi

serial :
	$(CC) $(CFLAGS) $(SERIAL_SRCS) -o BM_serial.ex
mpi:
	$(CC) $(CFLAGS) $(MPI_SRCS) -o BM_mpi.ex

clean :
	rm *.ex
