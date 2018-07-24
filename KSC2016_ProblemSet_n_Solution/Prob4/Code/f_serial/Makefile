.SUFFIXES : .f90 .f .o

FC   	= mpif90
FCFLAGS	= -O3 -fdefault-real-8 -fopenmp
OBJS 	= multigrid.o cgsolver.o matrixconstructor.o poisson.o
LIBS	:= -lm

all	: mytest

mytest	: $(OBJS)
	$(FC) -o poisson.ex $(FCFLAGS) $(OBJS) $(LIBS)

clean :
	rm -rf *.o *.mod poisson.ex result


.f90.o:
	$(FC) $(FCFLAGS) -c $<

