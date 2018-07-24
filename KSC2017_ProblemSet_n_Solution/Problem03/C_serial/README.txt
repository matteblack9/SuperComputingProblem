1. BUILD
	$ make serial

2. RUN
	$ ./BM_serial.ex [target file path] [pattern string]
	ex) ./BM_serial.ex ./PPAP.txt apple
	or
	$ mpirun -np 1 ./BM_serial.ex [target file path] [pattern string]
	ex) mpirun -np 1 ./BM_serial.ex ./PPAP.txt apple
