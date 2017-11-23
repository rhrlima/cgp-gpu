#ifndef HEADER_CGP_
#define HEADER_CGP_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "cgp.h"

#define ADD 0
#define SUB 1
#define MUL 2
#define DIV 3

/* Random int between [min, max] exclusive */
int randint(int min, int max) {

	return rand() % (max-min) + min;
}

/* Random float between [min, max] exclusive */
float randfloat(float min, float max) {

	return (rand() / (float) RAND_MAX) * (max - min) + min;
}

/* Prints the data of a solution [nodes] and [active nodes] */
void printSolution(GENES genes, int numGenes, int numInputs, int numOutputs) {
	for(int i=0; i<numGenes; i++) {
		printf("%d ", genes.nodes[i].func);
		for(int j=0; j<numInputs; j++) printf("%d ", genes.nodes[i].links[j]);
	}
	for(int i=0; i<numOutputs; i++) {
		printf("%d ", genes.nodes[numGenes+i].func);
	}
	printf("| ");
	for(int i=0; i<numGenes; i++) {
		printf("%d ", genes.toEvaluate[i]);
	}
	printf("\n");
}

// void createNode(NODE *node, int index, int numRows, int levelsBack, int numFunctions, int numInputs) {

// 	node -> func = randint(0, numFunctions);
// 	node -> links = (int *) malloc(numInputs * sizeof(int));

// 	int col = index / numRows;
// 	int min = numInputs + (col - levelsBack) * numRows;
// 	int max = numInputs + col * numRows;
// 	for(int l=0; l<numInputs; l++) {
// 		int value = (col >= levelsBack) ? randint(min, max) : randint(0, max);
// 		node -> links[l] = value;
// 	}
// }

NODE createNode(NODE *node, int index, int numRows, int levelsBack, int numFunctions, int numInputs) {

	node -> func = randint(0, numFunctions);
	node -> links = (int *) malloc(numInputs * sizeof(int));

	int col = index / numRows;
	int min = numInputs + (col - levelsBack) * numRows;
	int max = numInputs + col * numRows;
	for(int l=0; l<numInputs; l++) {
		int value = (col >= levelsBack) ? randint(min, max) : randint(0, max);
		node -> links[l] = value;
	}
}

void createOutputNode(NODE *node, int maxValue) {
	node -> func = randint(0, maxValue);
	node -> links = NULL;
}

// void createSolution(GENES *genes, CONFIG params) {

// 	genes -> nodes = (NODE *) malloc((params.numGenes + params.numOutputs) * sizeof(NODE));
// 	genes -> toEvaluate = (int *) malloc( params.numGenes * sizeof(int));
// 	genes -> nodesOutput = (double *) malloc(params.numOutputs * sizeof(double));

// 	for(int i=0; i<params.numGenes; i++) {
// 		NODE n;
// 		createNode(&n, i, params.numRows, params.levelsBack, params.numFunctions, params.numInputs);
// 		genes -> nodes[i] = n;
// 		genes -> toEvaluate[i] = 0;
// 	}
// 	for(int i=0; i<params.numOutputs; i++) {
// 		NODE n;
// 		createOutputNode(&n, params.numInputs + params.numGenes );
// 		genes -> nodes[params.numGenes + i] = n;
// 		genes -> nodesOutput[i] = 0;
// 	}
// }

GENES createSolution(CONFIG params) {
	GENES solution;
	int numGenes = params.numGenes;
	int numInputs = params.numInputs;
	int numOutputs = params.numOutputs;

	solution.nodes = (NODE *) malloc( (numGenes + numOutputs) * sizeof(NODE) );
	solution.toEvaluate = (int *) malloc(numGenes * sizeof(int));
	solution.nodesOutput = (double *) malloc(numOutputs * sizeof(double));

	for(int i=0; i<numGenes; i++) {
		NODE n;
		createNode(&n, i, params.numRows, params.levelsBack, params.numFunctions, params.numInputs);
		solution.nodes[i] = n;
		solution.toEvaluate[i] = 0;
	}
	for(int i=0; i<numOutputs; i++) {
		NODE n;
		createOutputNode(&n, numInputs + numGenes );
		solution.nodes[numGenes + i] = n;
		solution.nodesOutput[i] = 0;
	}
	return solution;
}

void hardCodedSolution(GENES *genes, int numGenes, int numOutputs) {
	
	genes -> nodes = (NODE *) malloc( (numGenes + numOutputs) * sizeof(NODE) );
	genes -> toEvaluate = (int *) malloc( numGenes * sizeof(int));
	genes -> nodesOutput = (double *) malloc(numOutputs * sizeof(double));

	int vet[22] = {0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 2, 5, 7, 3};

	for(int i=0; i<numGenes; i++) {
		NODE aux;
		aux.links = (int *) malloc(2 * sizeof(int));
		aux.func = vet[i*3];
		aux.links[0] = vet[i*3+1];
		aux.links[1] = vet[i*3+2];

		genes -> nodes[i] = aux;
		genes -> toEvaluate[i] = 0;
	}
	
	for(int i=0; i<numOutputs; i++) {
		NODE aux;
		aux.func = vet[numGenes*3+i];
		aux.links = NULL;
		genes -> nodes[numGenes+i] = aux;
		genes -> nodesOutput[i] = 0;
	}
}

void getActiveNodes(GENES *genes, int numInputs, int numOutputs, int numGenes) {

	for(int i=0; i<numOutputs; i++) {
		if (genes -> nodes[ numGenes + i ].func >= numInputs)
			genes -> toEvaluate[ genes -> nodes[ numGenes + i ].func - numInputs ] = 1;
	}

	for(int i=numGenes-1; i>=0; i--) {
		if (genes -> toEvaluate[i]) {
			for (int j=0; j<numInputs; j++) {
				if (genes -> nodes[i].links[j] >= numInputs)
					genes -> toEvaluate[ genes -> nodes[i].links[j] - numInputs ] = 1;
			}
		}
	}
}

double executeNode(int func, double *inputs, int numInputs) {
	double result = inputs[0];
	switch(func) {
		case ADD:
			for(int i=1; i<numInputs; i++) {
				result += inputs[i];
			}
			break;
		case SUB:
			for(int i=1; i<numInputs; i++) {
				result -= inputs[i];
			}
			break;
		case MUL:
			for(int i=1; i<numInputs; i++) {
				result *= inputs[i];
			}
			break;
		case DIV:
			for(int i=1; i<numInputs; i++) {
				result /= (inputs[i] > 0.0) ? inputs[i] : 1.0;
			}
			break;
		default:
			printf("Eu não devia ter sido printado!");
			break;
	}
	return result;
}

void executeNodes(GENES *genes, int numGenes, int numInputs, int numOutputs, double *inputData) {

	double outputs[numInputs + numGenes];
	double inputs[numInputs];

	for (int i=0; i<numInputs; i++) {
		outputs[i] = inputData[i];
	}
	for (int i=0; i<numGenes; i++) {
		if (genes -> toEvaluate[i]) {
			for(int j=0; j<numInputs; j++) {
				inputs[j] = outputs[ genes -> nodes[i].links[j] ];
			}
			outputs[numInputs + i] = executeNode(genes -> nodes[i].func, inputs, numInputs);
		}
	}
	for(int i=0; i<numOutputs; i++) {
		genes -> nodesOutput[i] = outputs[ genes -> nodes[ numGenes + i ].func ];
	}
}

double calculateFitness(GENES *genes, int numInputs, int numOutputs, int numGenes, double inputData[][2], double *outputData, int numCases) {
	double rss = 0.0;

	getActiveNodes(genes, numInputs, numOutputs, numGenes);

	for(int i=0; i<numCases; i++) {
		executeNodes(genes, numGenes, numInputs, numOutputs, inputData[i]);
		double temp = 0.0;
		for(int j=0; j<numOutputs; j++) {
			temp += genes -> nodesOutput[j];
		}
		rss += pow(temp - outputData[i], 2);
	}
	return rss;
}

void freeSolution(GENES genes, int numGenes) {
	for(int i=0; i<numGenes; i++) {
		free(genes.nodes[i].links);
	}
	free(genes.nodes);
	free(genes.toEvaluate);
	free(genes.nodesOutput);
}

GENES * createPopulation(int popSize, CONFIG params) {
	GENES *population = (GENES *) malloc(popSize * sizeof(GENES));

	for(int i=0; i<popSize; i++)
		population[i] = createSolution(params);

	return population;
}

void mutation(GENES parent, GENES *offspring, CONFIG params) {

	offspring -> nodes = (NODE *) malloc( (params.numGenes + params.numOutputs) * sizeof(NODE) );
	offspring -> toEvaluate = (int *) malloc( params.numGenes * sizeof(int));
	offspring -> nodesOutput = (double *) malloc(params.numOutputs * sizeof(double));

	int rand1 = randint(0, params.numGenes + params.numOutputs);

	for(int i=0; i<params.numGenes+params.numOutputs; i++) {
		offspring -> nodes[i] = parent.nodes[i];
		if(i < params.numGenes)
			offspring -> toEvaluate[i] = 0;
		else
			offspring -> nodesOutput[i] = 0;
	}

	if (rand1 < params.numGenes) {
		int rand2 = randint(0, params.numInputs + 1);
		if (rand2 == 0) {
			int old = offspring -> nodes[rand1].func;
			int new = randint(0, params.numFunctions);
			offspring -> nodes[rand1].func = new;
			printf("Old: %d New: %d\n", old, new);
		} else {
			int col = rand1 / params.numRows;
			int min = params.numInputs + (col - params.levelsBack) * params.numRows;
			int max = params.numInputs + col * params.numRows;

			int old = offspring -> nodes[rand1].links[rand2];
			int new = (col >= params.levelsBack) ? randint(min, max) : randint(0, max);
			offspring -> nodes[rand1].links[rand2] = new;
			printf("Old: %d New: %d\n", old, new);
		}
	} else {
		int old = offspring -> nodes[rand1].func;
		int new = randint(0, params.numInputs + params.numGenes);
		offspring -> nodes[rand1].func = new;
		printf("Old: %d New: %d\n", old, new);
	}
}

void execute(CONFIG params) {

	GENES *population;

}

#endif