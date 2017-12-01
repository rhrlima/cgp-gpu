echo "Compilando versão CPU:"

echo "Exemplo de execução (main)"
gcc cpu/cgp.c cpu/main.c -o main -lm

echo "Experimentos (cpu_test)"
gcc cpu/cgp.c tests/cpu_test.c -o cpu_test -lm

echo ""
echo "Compilando versão GPU:"

echo "Exemplo de execução (gmain)"
nvcc cuda/cgp.cu cuda/main.cu -o gmain -Wno-deprecated-gpu-targets

echo "Experimentos (gpu_test)"
nvcc cuda/cgp.cu tests/gpu_test.cu -o gpu_test -Wno-deprecated-gpu-targets

echo "Pronto"