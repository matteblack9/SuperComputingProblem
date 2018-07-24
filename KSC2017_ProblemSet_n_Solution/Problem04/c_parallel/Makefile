C    	= mpicc
CFLAGS  = -O3 -std=c99 -fopenmp
OBJS 	= c_parallel.o
LIBS	:= -lm

all	: mytest

mytest	: $(OBJS)
	$(C) -o c_parallel.ex $(CFLAGS) $(OBJS) $(LIBS)

clean :
	rm -rf *.o c_parallel.ex result *.out

.SUFFIXES : .c .o

.c.o:   ; $(C) -c $(CFLAGS) $(INC) $?


