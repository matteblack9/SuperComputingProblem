.SUFFIXES : .f90 .f .o

FC   	= mpif90
FCFLAGS	= -O3 -fdefault-real-8 -fopenmp
OBJS 	= multigrid.o main.o
LIBS	:= -lm

all	: mytest

mytest	: $(OBJS)
	$(FC) -o multigrid.ex $(FCFLAGS) $(OBJS) $(LIBS)

clean :
	rm -rf *.o *.mod multigrid.ex.ex solution*


.f90.o:
	$(FC) $(FCFLAGS) -c $<

