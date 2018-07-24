CC		= mpicc
CFLAGS	= -O3 -std=gnu99 -fopenmp
LDFLAGS	= $(CFLAGS)

TARGET	= kmeans.x
OBJS	= kmeans.o util.o main.o

all: $(TARGET)

$(TARGET) : $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

%.o : %.c kmeans.h
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f *.o $(TARGET)
