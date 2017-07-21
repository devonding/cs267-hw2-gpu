CC = nvcc
MPCC = nvcc
OPENMP = 
CFLAGS = -O3 -arch=sm_35
NVCCFLAGS = -O3 -arch=sm_35
LIBS = -lm

TARGETS = gpu

all:	$(TARGETS)

gpu: gpu.o common.o 
	$(CC) -o $@ $(NVCCLIBS) gpu.o common.o 

gpu.o: gpu.cu common.h
	$(CC) -c $(NVCCFLAGS) gpu.cu
common.o: common.cu common.h
	$(CC) -c $(CFLAGS) common.cu

clean:
	rm -f *.o $(TARGETS)
