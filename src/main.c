#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "cgp.h"

//add modular functions
//add support to different datasets

int main() {

	srand(time(NULL));

	//DATASET dset;
	struct parameters *params;
	struct chromosome *chromo;

	//dset = loadDataset("datasets/symbolic2.data"); //(x*x) + (x+x) + (x-x)

	params = (struct parameters*)malloc(sizeof(struct parameters));

	params->numInputs = 1;
	params->numNodes = 6;
	params->numOutputs = 1;
	params->arity = 2;
	params->numFunctions = 4;
	//params.numCols = 3;
	//params.numRows = 3;
	//params.levelsBack = 9;

	double inputs[] = {1, 2, 3};
	
	chromo = createChromosome(params);
	executeChromosome(chromo, inputs);
	printChromosome(chromo);

	// calculateFitness(&G, params.numInputs, params.numOutputs, params.numGenes, dset);
	// printSolution(G, params.numGenes, params.arity, params.numOutputs);
	// freeSolution(G, params.numGenes);
	// freeDataset(dset);

	// dset = loadDataset("datasets/symbolic3.data"); //(x*x) + (x+x) + (x-x)
	// params.numInputs = dset.numInputs; 
	// params.numOutputs = dset.numOutputs;
	// params.numFunctions = 4;
	// params.arity = 2;
	// params.numCols = 3;
	// params.numRows = 2;
	// params.levelsBack = 6;
	// params.numGenes = params.numCols * params.numRows;
	
	// G = createSolution(params);
	// calculateFitness(&G, params.numInputs, params.numOutputs, params.numGenes, dset);
	// printSolution(G, params.numGenes, params.arity, params.numOutputs);
	// freeSolution(G, params.numGenes);
	// freeDataset(dset);

	//execute(5, 1000, params, dset);

	// int v1[] = {2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 5, 4, 0, 0, 0, 0, 0, 0, 8};
	// GENES G1 = hardCodedSolution(v1, c * r, i, o);
	// calculateFitness(&G1, i, o, c * r, dset);

	// int v2[] = {0, 1, 0, 2, 0, 0, 1, 0, 1, 0, 3, 4, 0, 2, 4, 0, 2, 3, 0, 2, 5, 0, 4, 7, 3, 5, 5, 9};
	// GENES G2 = hardCodedSolution(v2, c * r, i, o);
	// calculateFitness(&G2, i, o, c * r, dset);

	return 0;
}