C    	= mpicc
CFLAGS  = -O3 -std=c99 
OBJS 	= cgsolver_mpi.o matrixconstructor.o poisson_mpi.o
LIBS	:= -lm

all	: mytest

mytest	: $(OBJS)
	$(C) -o c_parallel.ex $(CCFLAGS) $(OBJS) $(LIBS)

clean :
	rm -rf *.o c_parallel.ex result

.SUFFIXES : .c .o

.c.o:   ; $(C) -c $(CFLAGS) $(INC) $?


