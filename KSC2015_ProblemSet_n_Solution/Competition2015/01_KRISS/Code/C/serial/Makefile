CC=mpicc
CPP = mpic++
CFLAGS=-O3 

.PHONY: all

all: cserial
cserial: b_serial.c
	$(CC) $(CFLAGS) b_serial.c -lm -std=c99 -o c_serial.ex
clean:
	rm -rf *.o *.mod ?_*l.ex
