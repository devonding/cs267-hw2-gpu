# Load CUDA using the following command
# module load cuda

#
# Stampede
#
CC = nvcc
MPCC = nvcc
OPENMP = 
CFLAGS = -O3 -arch=sm_35
NVCCFLAGS = -O3 -arch=sm_35
LIBS = -lm

TARGETS = serial gpu_serial gpu

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o
gpu_serial: gpu_serial.o common.o
	$(CC) -o $@ $(NVCCLIBS) gpu_serial.o common.o
gpu: gpu.o common.o 
	$(CC) -o $@ $(NVCCLIBS) gpu.o common.o 

serial.o: serial.cu common.h
	$(CC) -c $(CFLAGS) serial.cu
gpu_serial.o: gpu_serial.cu common.h
	$(CC) -c $(NVCCFLAGS) gpu_serial.cu
gpu.o: gpu.cu common.h
	$(CC) -c $(NVCCFLAGS) gpu.cu
common.o: common.cu common.h
	$(CC) -c $(CFLAGS) common.cu

clean:
	rm -f *.o $(TARGETS)
