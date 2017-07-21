#
# Hopper - NERSC 
#
# Portland Group Compilers PGI are loaded by default; for other compilers please check the module list
#
CC = g++
#CC 
MPCC = CC
OPENMP = -mp
CFLAGS = -O1 
# here is the O3 should be
LIBS =


TARGETS = fast_serial autograder

all:	$(TARGETS)


autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) autograder.o common.o grid.o
fast_serial: fast_serial.o common.o grid.o
	$(CC) -o $@ $(LIBS) fast_serial.o common.o grid.o

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp
common.o: common.cpp common.h grid.cpp
	$(CC) -c $(CFLAGS) common.cpp grid.cpp
fast_serial.o: fast_serial.cpp common.h 
	$(CC) -c $(CFLAGS) fast_serial.cpp 
grid.o: grid.cpp common.h
	$(CC) -c $(CFLAGS) grid.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
