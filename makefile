CC = gcc
CUC = nvcc

FLAGS = -lm
GPUFLAGS = -Wno-deprecated-gpu-targets

all:
	$(CC) src/cgp.c src/main.c -o main $(FLAGS)

run:
	$(CC) src/cgp.c src/main.c -o main $(FLAGS)
	./main

gpu:
	$(CUC) cuda/cgp.cu cuda/main.cu -o gmain $(GPUFLAGS)
	./gmain

test:
	$(CC)  cpu/cgp.c tests/cpu_teste.c -o cpu_test $(FLAGS)
	$(CUC) cuda/cgp.cu tests/gpu_teste.cu -o gpu_test $(GPUFLAGS)

clean:
	rm ./*.exe