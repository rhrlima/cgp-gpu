#include <time.h>
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

	time_t start, end;
	double interval;

	int i;
	for(i = 0; i < RUNS; i++) {
		start = clock();
		chromo = executeCGP(params, data, POPSIZE, MAXGENS);
		end = clock();
		interval = (double) (end - start) / CLOCKS_PER_SEC;
		printf("%f\n", interval);
	}

	freeChromosome(chromo);
	freeDataset(data);
	free(params);

	return 0;
}