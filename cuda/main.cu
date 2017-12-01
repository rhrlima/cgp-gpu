#include "cgp.cuh"

#define DATASETBUFFER 100

#define NUMNODES 9
#define MAXARITY 2
#define NUMFUNCTIONS 4

#define POPSIZE 5
#define MAXGENS 100

int main(int argc, char *argv[]) {

	unsigned int seed = time(NULL);
	srand(seed);

	char dataset_file[DATASETBUFFER];

	struct dataset *data;
	struct parameters *params;

	int *result;
	struct chromosome *chromo;

	if (argc > 1) strcpy(dataset_file, argv[1]);
	else {
		printf("No dataset was specified.\nExiting.\n");
		exit(0);
	}

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