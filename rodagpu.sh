echo -e "Experimentos para a GPU, o algoritmo é executado 10x e o tempo em segundos é obtido.\n"

echo -e "DATASET symbolic:\n"

echo -e "\nTeste 1000 instancias"
./gpu_test datasets/symbolic_1000.data

echo -e "\nTeste 10000 instancias"
./gpu_test datasets/symbolic_10000.data

echo -e "\nTeste 100000 instancias"
./gpu_test datasets/symbolic_100000.data

echo -e "\nTeste 1000000 instancias"
./gpu_test datasets/symbolic_1000000.data



echo -e "\n\nDATASET symbolic2:\n"

echo -e "\nTeste 1000 instancias"
./gpu_test datasets/symbolic2_1000.data

echo -e "\nTeste 10000 instancias"
./gpu_test datasets/symbolic2_10000.data

echo -e "\nTeste 100000 instancias"
./gpu_test datasets/symbolic2_100000.data

echo -e "\nTeste 1000000 instancias"
./gpu_test datasets/symbolic2_1000000.data

echo -e "\nPronto\n"