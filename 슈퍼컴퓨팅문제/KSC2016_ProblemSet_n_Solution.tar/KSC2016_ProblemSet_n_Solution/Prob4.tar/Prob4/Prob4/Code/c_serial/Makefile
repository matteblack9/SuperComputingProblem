C    	= mpicc
CFLAGS  = -O3 -std=c99 -fopenmp
OBJS 	= multigrid.o cgsolver.o matrixconstructor.o poisson.o
LIBS	:= -lm

all	: mytest

mytest	: $(OBJS)
	$(C) -o poisson.ex $(CFLAGS) $(OBJS) $(LIBS)

clean :
	rm -rf *.o poisson.ex result

.SUFFIXES : .c .o

.c.o:   ; $(C) -c $(CFLAGS) $(INC) $?


