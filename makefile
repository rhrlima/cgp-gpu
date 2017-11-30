CC = gcc
CUC = nvcc

FLAGS = -lm
GPUFLAGS = -Wno-deprecated-gpu-targets

all: src/cgp.c src/cgp.h src/main.c
	$(CC) src/cgp.c src/main.c -o main $(FLAGS)

run: src/cgp.c src/cgp.h src/main.c
	$(CC) src/cgp.c src/main.c -o main $(FLAGS)
	./main

test: src/cgp.c src/cgp.h tests/test1.c tests/test2.c 
	$(CC) src/cgp.c tests/test1.c -o test1
	$(CC) src/cgp.c tests/test2.c -o test2
	./test1
	./test2

gpu: cuda/cgp.cu cuda/cgp.cuh cuda/main.cu
	$(CUC) cuda/cgp.cu cuda/main.cu -o gmain $(GPUFLAGS)
	./gmain

clean:
	rm -f main