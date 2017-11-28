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
	struct chromosome *chromo, *best;

	strcpy(dataset_file, "datasets/symbolic2.data");	//x*x + x+x | (9, 2, 4)
	//strcpy(dataset_file, "datasets/symbolic3.data");	//(x0+x1) + (x0*x1) + (-x0)*(x1*x1) | (8, 2, 4)
	//strcpy(dataset_file, "datasets/symbolic4.data"); //x0*x1*x1 + x2*x1 + x*3

	data = loadDataset(dataset_file);

	params = initialiseParameters(9, 2, 4, data);//numNodes, maxArity, numFunctions

	printf("Dataset: '%s'\n", dataset_file);

	printParameters(params);

	printf("Running CGP\n");
	chromo = executeCGP(params, data, 10000);

	printf("Best solution found\n");
	printChromosome(chromo);

	int array[] = {1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 2, 4, 1, 4, 1, 2, 2, 5, 3, 2, 3, 6}; //best for dataset 2
	// int array[] = {0, 0, 1, 2, 0, 1, 1, 0, 0, 1, 4, 0, 2, 1, 1, 2, 5, 6, 0, 2, 3, 0, 8, 7, 9}; //best for dataset 3
	// int array[] = {2, , 1, 2};
	best = createChromosomeFromArray(params, array);
	calculateFitness(best, data);

	printf("Best hardcoded\n");
	printChromosome(best);

	freeChromosome(chromo);
	freeChromosome(best);
	
	freeDataset(data);
	free(params);

	return 0;
}