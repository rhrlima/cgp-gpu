#include <time.h>
#include "../cuda/cgp.cuh"

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

	int *result;

	if (argc > 1) strcpy(dataset_file, argv[1]);
	else exit(0);

	data = loadDataset(dataset_file);
	
	params = initialiseParameters(NUMNODES, MAXARITY, NUMFUNCTIONS, data);

	clock_t start, end;
	double interval;

	int i;
	for(i = 0; i < RUNS; i++) {
		start = clock();
		result = CUDAexecuteCGP(params, data, POPSIZE, MAXGENS);
		end = clock();
		interval = (double) (end - start) / CLOCKS_PER_SEC;
		printf("%f\n", interval);
	}

	free(result);
	freeDataset(data);
	free(params);

	return 0;
}