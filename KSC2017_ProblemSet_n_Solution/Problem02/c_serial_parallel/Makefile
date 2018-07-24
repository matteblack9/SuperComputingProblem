CC = mpicc
TARGET = bhlab bhlab_parallel bhlab_parallel_
OBJECTS = blackhole_lab.o SBLspacetime.o
CFLAGS = -std=c99 -fopenmp
LINK = -lm

all : $(TARGET)

bhlab : main.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LINK)

bhlab_parallel : main_parallel.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LINK)

bhlab_parallel_ : main_parallel_.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LINK)

clean :
	rm $(TARGET) $(OBJECTS) main.o main_parallel.o main_parallel_.o
