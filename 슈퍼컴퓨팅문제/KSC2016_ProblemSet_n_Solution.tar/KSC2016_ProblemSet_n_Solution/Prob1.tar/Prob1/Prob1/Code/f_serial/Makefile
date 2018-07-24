CC		= mpif90
CFLAGS	= -O3 -fopenmp
LDFLAGS	= $(CFLAGS)

TARGET	= kmeans.x
OBJS	= util.o kmeans.o main.o

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

%.o : %.f90
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f *.o *.mod $(TARGET)
