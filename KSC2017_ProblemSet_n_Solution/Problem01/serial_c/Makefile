.SUFFIXES : .cpp .o

CC    	= gcc
CFLAGS  = -O2 -m64
OBJS 	= motion_estimation.o main.o
LIBS	= -lm
TARGET = c_serial.ex

$(TARGET) : $(OBJS)
	$(CC) -o $(TARGET) $(CCFLAGS) $(OBJS) $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o $(TARGET) result

