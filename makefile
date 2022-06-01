CC=g++
CFLAGS=-std=c++11
MPICC=mpic++
OBJ = main.o

main: main.o
	$(MPICC) -o main main.o

main.o: main.cpp
	$(MPICC) $(CFLAGS) -c main.cpp

.PHONY: clean

clean:
	rm main main.o helpers.h.gch
