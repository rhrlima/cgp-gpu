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
	
	pop = createPopulation(5, params);
	for(int i=0; i<5; i++)
		printSolution(pop[i], 6, 2, 4);

	//Y = X0+X1 + X0*X1 + -X0*(X1*X1) + 0
	int fitCases = 2;
	double in[2][2] = {{2.0, 2.0}, {1.0, 3.0}};
	double out[] = {0, -2};
	double fit;

	hardCodedSolution(&G1, params.numGenes, params.numOutputs);
	fit = calculateFitness(&G1, params.numInputs, params.numOutputs, params.numGenes, in, out, fitCases);
	printSolution(G1, 6, 2, 4);
	// printf("Fitness: %2.2f\n", fit);

	// createSolution(&G2, params);
	// fit = calculateFitness(&G2, params.numInputs, params.numOutputs, params.numGenes, in, out, fitCases);
	// printSolution(G2, params.numGenes, params.numOutputs);
	// printf("Fitness: %2.2f\n", fit);

	// mutation(G1, &offspring, params);
	// fit = calculateFitness(&offspring, params.numInputs, params.numOutputs, params.numGenes, in, out, fitCases);
	// printSolution(offspring, 6, 2, 4);
	// printf("Fitness: %2.2f\n", fit);

	//freeSolution(G1, params.numGenes);
	//freeSolution(G2, params.numGenes);
	//freeSolution(offspring, params.numGenes);

	return 0;
}