C    	= mpicc
CFLAGS  = -O3 -std=c99 -fopenmp 
OBJS 	= multigrid.o main.o
LIBS	:= -lm

all	: mytest

mytest	: $(OBJS)
	$(C) -o multigrid.ex $(CFLAGS) $(OBJS) $(LIBS)

clean :
	rm -rf *.o multigrid.ex solution*

.SUFFIXES : .c .o

.c.o:   ; $(C) -c $(CFLAGS) $(INC) $?


