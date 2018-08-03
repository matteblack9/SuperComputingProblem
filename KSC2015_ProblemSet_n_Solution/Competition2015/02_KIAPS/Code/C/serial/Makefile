CC=gcc
MCC = mpicc
CFLAGS=-O3 

.PHONY: all

all: cserial

cserial: 1d_serial.c
	$(CC) $(CFLAGS) 1d_serial.c -lm -std=c99 -o c_serial.ex
	$(MCC) $(CFLAGS) test.c -lm -std=c99 -o test.ex
clean:
	rm -rf *.o *.mod ?_*l.ex
