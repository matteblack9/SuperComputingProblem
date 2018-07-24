C    	= gcc
CFLAGS  = -O3 -std=c99 
OBJS 	= cgsolver.o matrixconstructor.o poisson.o
LIBS	:= -lm

all	: mytest

mytest	: $(OBJS)
	$(C) -o c_serial.ex $(CCFLAGS) $(OBJS) $(LIBS)

clean :
	rm -rf *.o c_serial.ex result

.SUFFIXES : .c .o

.c.o:   ; $(C) -c $(CFLAGS) $(INC) $?


