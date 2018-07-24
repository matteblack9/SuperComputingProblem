F    	= mpif90
FFLAGS = -O3 -w 
OBJS 	= poisson_mpi.o cgsolver_mpi.o matrixconstructor.o
LIBS	:= -lm

all	: mytest

mytest	: $(OBJS)
	$(F) -o f_parallel.ex $(CCFLAGS) $(OBJS) $(LIBS)

clean :
	rm -rf *.o f_parallel.ex result

.SUFFIXES : .f90 .o

.f90.o:   ; $(F) -c $(FFLAGS) $(INC) $?


