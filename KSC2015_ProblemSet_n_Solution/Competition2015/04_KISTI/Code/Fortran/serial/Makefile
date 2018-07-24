F    	= ifort
FFLAGS  = -O3 -w
OBJS 	= poisson.o cgsolver.o matrixconstructor.o
LIBS	:= -lm

all	: mytest

mytest	: $(OBJS)
	$(F) -o f_serial.ex $(FFLAGS) $(OBJS) $(LIBS)

clean :
	rm -rf *.o f_serial.ex result

.SUFFIXES : .f90 .o

.f90.o:   ; $(F) -c $(FFLAGS) $(INC) $?


