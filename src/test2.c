#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cgp.h"

int main() {

	srand(time(NULL));

	char dataset_file[100];

	struct dataset *data;
	struct parameters *params;
	struct chromosome *chromo;

	strcpy(dataset_file, "datasets/symbolic3.data"); //(x0+x1) + (x0*x1) + (-x0)*(x1*x1) | (8, 2, 4)

	data = loadDataset(dataset_file);

	params = initialiseParameters(8, 2, 4, data);

	printf("Dataset: '%s'\n", dataset_file);

	printParameters(params);

	printf("Running CGP\n");
	chromo = executeCGP(params, data, 10000);

	printf("Best solution found\n");
	printChromosome(chromo);

	int array[] = {0, 0, 1, 2, 0, 1, 1, 0, 0, 1, 4, 0, 2, 1, 1, 2, 5, 6, 0, 2, 3, 0, 8, 7, 9};
	chromo = createChromosomeFromArray(params, array);
	calculateFitness(chromo, data);

	printf("Best hardcoded\n");
	printChromosome(chromo);

	return 0;
}