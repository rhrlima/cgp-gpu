
CPUC = gcc
GPUC = nvcc

CPUFLAGS = -lm
GPUFLAGS = -Wno-deprecated-gpu-targets

all:
	$(info Compilando TUDO (exemplos e testes) )
	$(CPUC) cpu/cgp.c cpu/main.c -o main.o $(CPUFLAGS)
	$(CPUC) cpu/cgp.c tests/cpu_test.c -o cpu_test.o $(CPUFLAGS)

	$(GPUC) cuda/cgp.cu cuda/main.cu -o gmain.o $(GPUFLAGS)
	$(GPUC) cuda/cgp.cu tests/gpu_test.cu -o gpu_test.o $(GPUFLAGS)

cpuonly:
	$(info Compilando CPU (exemplos e testes) )
	$(CPUC) cpu/cgp.c cpu/main.c -o main.o $(CPUFLAGS)
	$(CPUC) cpu/cgp.c tests/cpu_test.c -o cpu_test.o $(CPUFLAGS)

gpuonly:
	$(info Compilando GPU (exemplos e testes) )
	$(GPUC) cuda/cgp.cu cuda/main.cu -o gmain.o $(GPUFLAGS)
	$(GPUC) cuda/cgp.cu tests/gpu_test.cu -o gpu_test.o $(GPUFLAGS)

clean:
	$(info Limpando executáveis)
	rm ./*.o