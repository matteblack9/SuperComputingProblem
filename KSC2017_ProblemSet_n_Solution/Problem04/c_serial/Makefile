C    	= mpicc
CFLAGS  = -O3 -std=c99 -fopenmp
OBJS 	= c_serial.o
LIBS	:= -lm

all	: mytest

mytest	: $(OBJS)
	$(C) -o c_serial.ex $(CFLAGS) $(OBJS) $(LIBS)

clean :
	rm -rf *.o c_serial.ex result *.out

.SUFFIXES : .c .o

.c.o:   ; $(C) -c $(CFLAGS) $(INC) $?


