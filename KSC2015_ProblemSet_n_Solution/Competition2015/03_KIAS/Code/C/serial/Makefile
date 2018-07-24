CC=gcc
CFLAGS=-O3 

.PHONY: all

all: cserial

cserial: sequential.c
	$(CC) $(CFLAGS) sequential.c -lm -std=c99 -o c_serial.ex
clean:
	rm -rf *.o *.mod ?_*l.ex
