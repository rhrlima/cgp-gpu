#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cgp.h"

int main() {

	srand(time(NULL));

	char dataset_file[100];

	struct dataset *dset;
	struct parameters *params;
	struct chromosome *chromo, *chromo2;

	//strcpy(dataset_file, "datasets/symbolic.data");
	strcpy(dataset_file, "datasets/symbolic2.data");
	//strcpy(dataset_file, "datasets/symbolic3.data");

	printf("Dataset: %s\n", dataset_file);
	dset = loadDataset(dataset_file);

	params = (struct parameters*)malloc(sizeof(struct parameters));

	params->numInputs 	 = 1;
	params->numNodes 	 = 6;
	params->numOutputs 	 = 1;
	params->arity 		 = 2;
	params->numFunctions = 4;
	
	// chromo = createChromosome(params);
	// chromo2 = copyChromosome(chromo);

	// calculateFitness(chromo, dset);
	// calculateFitness(chromo2, dset);

	// printChromosome(chromo);
	// printChromosome(chromo2);

	// singleMutation(params, chromo2);
	// calculateFitness(chromo2, dset);

	// printChromosome(chromo);
	// printChromosome(chromo2);

	//From array chromosome
	int array[] = {0, 0, 0, 1, 0, 0, 1, 1, 1, 2, 0, 3, 0, 0, 0, 3, 4, 5, 5};
	chromo = createChromosomeFromArray(params, array);
	calculateFitness(chromo, dset);
	printChromosome(chromo);

	printf("Running CGP\n");
	chromo = executeCGP(params, dset, 1000);

	printf("Best solution found\n");
	printChromosome(chromo);

	return 0;
}