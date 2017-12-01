#include <windows.h>
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
	else exit(0);

	printf("Dataset: '%s'\n", dataset_file);
	data = loadDataset(dataset_file);
	
	params = initialiseParameters(NUMNODES, MAXARITY, NUMFUNCTIONS, data);
	printParameters(params);

	printf("CUDA CGP:\n");

	LARGE_INTEGER frequency;
    LARGE_INTEGER start;
    LARGE_INTEGER end;
    double interval;

    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&start);

	result = CUDAexecuteCGP(params, data, POPSIZE, MAXGENS);

	QueryPerformanceCounter(&end);
    interval = (double) (end.QuadPart - start.QuadPart) / frequency.QuadPart;

    printf("TIME: %f seconds\n", interval);

	chromo = createChromosomeFromArray(params, result);
	calculateFitness(chromo, data);
	printChromosome(chromo);

	free(result);
	freeChromosome(chromo);
	freeDataset(data);
	free(params);

	return 0;
}