.SUFFIXES : .cpp .o

CC    	= mpicc
CFLAGS  = -O2 -m64
OBJS 	= motion_estimation.o main.o
LIBS	= -lm
TARGET = c_parallel.ex

$(TARGET) : $(OBJS)
	$(CC) -o $(TARGET) $(CCFLAGS) $(OBJS) $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o $(TARGET) result

