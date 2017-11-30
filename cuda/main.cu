#include "cgp.cuh"

#define DATASETBUFFER 100

#define NUMNODES 9
#define MAXARITY 2
#define NUMFUNCTIONS 4

#define POPSIZE 5
#define MAXGENS 100

int main() {

	unsigned int seed = time(NULL);
	srand(seed);

	char dataset_file[DATASETBUFFER];

	struct dataset *data;
	struct parameters *params;

	int *result;
	struct chromosome *chromo;

	strcpy(dataset_file, "datasets/symbolic2_1024.data");

	printf("Dataset: '%s'\n", dataset_file);
	data = loadDataset(dataset_file);
	
	params = initialiseParameters(NUMNODES, MAXARITY, NUMFUNCTIONS, data);
	printParameters(params);

	printf("CUDA CGP:\n");

	result = CUDAexecuteCGP(params, data, POPSIZE, MAXGENS);

	chromo = createChromosomeFromArray(params, result);
	calculateFitness(chromo, data);
	printChromosome(chromo);

	free(result);
	freeChromosome(chromo);
	freeDataset(data);
	free(params);

	return 0;
}