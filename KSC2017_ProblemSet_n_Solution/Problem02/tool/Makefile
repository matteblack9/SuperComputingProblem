CC = cc
TARGET = tool
OBJECTS = main.o
#IDIR = -I/usr/local/include
CFLAGS = -std=c99 $(IDIR)
#LDIR = -L/usr/local/lib
LINK = -lpng -lm

all : $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDIR) $(LINK)

clean :
	rm $(TARGET) $(OBJECTS)
