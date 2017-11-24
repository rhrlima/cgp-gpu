#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "cgp.h"

int main() {

	srand(time(NULL)); //No seed
	//srand(2);

	GENES *pop;

	GENES G1, G2, offspring;

	CONFIG params;
	params.numInputs = 2;
	params.numOutputs = 4;
	params.numFunctions = 4;
	params.functions;
	params.numCols = 3;
	params.numRows = 2;
	params.levelsBack = 5;
	params.numGenes = 6;

	//Y = X0+X1 + X0*X1 + -X0*(X1*X1) + 0
	DATASET dset = loadDataset("oi");

	pop = createPopulation(5, params);
	double best = pop[0].fitness;
	for(int i=0; i<5; i++) {
		calculateFitness(&pop[i], params.numInputs, params.numOutputs, params.numGenes, dset);
		printSolution(pop[i], 6, 2, 4);
		if (pop[i].fitness < best) best = pop[i].fitness;
	}
	for(int i=0; i<10; i++) {
		printf("Gen: %d\n", i);
		for(int j=0; j<5; j++) {
			if (randfloat(0, 1) < 0.1) {
				printf("Mutate: %d\n", j);
				GENES offspring = mutation(pop[j], params);
				calculateFitness(&offspring, params.numInputs, params.numOutputs, params.numGenes, dset);
				if (offspring.fitness < pop[i].fitness) {
					printf("Improved\n");
					pop[j] = offspring;
				}
			}
		}
		for(int j=0; j<5; j++) {
			if (pop[j].fitness < best) best = pop[j].fitness;
		}
		printf("Best: %f\n", best);
	}

	printf("\n");
	hardCodedSolution(&G1, params.numGenes, params.numOutputs);
	calculateFitness(&G1, params.numInputs, params.numOutputs, params.numGenes, dset);
	printSolution(G1, 6, 2, 4);

	return 0;
}