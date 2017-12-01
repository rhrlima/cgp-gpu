#include <windows.h>
#include "../cpu/cgp.h"

#define DATASETBUFFER 100

#define NUMNODES 9
#define MAXARITY 2
#define NUMFUNCTIONS 4

#define POPSIZE 5
#define MAXGENS 100

#define RUNS 10

int main(int argc, char *argv[]) {

	unsigned int seed = time(NULL);
	srand(seed);

	char dataset_file[DATASETBUFFER];

	struct dataset *data;
	struct parameters *params;
	struct chromosome *chromo;

	if (argc > 1) strcpy(dataset_file, argv[1]);
	else exit(0);

	data = loadDataset(dataset_file);

	params = initialiseParameters(NUMNODES, MAXARITY, NUMFUNCTIONS, data);

	LARGE_INTEGER frequency;
	LARGE_INTEGER start;
	LARGE_INTEGER end;
	double interval;

	for(int i=0; i<RUNS; i++) {
		QueryPerformanceFrequency(&frequency);
		QueryPerformanceCounter(&start);
		chromo = executeCGP(params, data, POPSIZE, MAXGENS);
		QueryPerformanceCounter(&end);
		interval = (double) (end.QuadPart - start.QuadPart) / frequency.QuadPart;
		printf("%f\n", interval);
	}

	freeChromosome(chromo);
	freeDataset(data);
	free(params);

	return 0;
}