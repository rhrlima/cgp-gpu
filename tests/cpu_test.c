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
	else {
		printf("No dataset specified.\nExiting.\n");
		exit(0);
	}

	data = loadDataset(dataset_file);

	params = initialiseParameters(NUMNODES, MAXARITY, NUMFUNCTIONS, data);

	// chromo = createChromosome(params);
	// calculateFitness(chromo, data);
	// printChromosome(chromo);

	chromo = executeCGP(params, data, 10, 1000);
	printChromosome(chromo);

	int sol[] = {2, 0, 0, 3, 1, 0, 1, 0, 0, 0, 0, 2, 0, 1, 4, 0, 4, 1, 2, 3, 1, 3, 5, 6, 2, 5, 4, 6};
	struct chromosome *chromo2 = createChromosomeFromArray(params, sol);
	calculateFitness(chromo2, data);
	printChromosome(chromo2);

	freeChromosome(chromo);
	freeDataset(data);
	free(params);

	return 0;
}