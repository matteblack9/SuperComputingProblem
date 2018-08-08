CC=gcc
MCC = mpicc
CFLAGS=-O3 

.PHONY: all

all: cserial

cserial: sequential.c
	$(CC) $(CFLAGS) sequential.c -lm -std=c99 -o c_serial.ex
	$(MCC) $(CFLAGS) test.c -lm -std=c99 -o test.ex
clean:
	rm -rf *.o *.mod ?_*l.ex
