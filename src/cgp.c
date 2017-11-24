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

DATASET loadDataset(char *fileName) {
	DATASET dset;

	dset.numInputs = 2;
	dset.numOutputs = 4;
	dset.numCases = 2;

	dset.inputs = (double **) malloc(dset.numCases * sizeof(double));
	for(int i=0; i<dset.numInputs; i++) dset.inputs[i] = (double *) malloc(dset.numInputs * sizeof(double));
	dset.outputs = (double *) malloc(dset.numCases * sizeof(double));

	dset.inputs[0][0] = 2;
	dset.inputs[0][1] = 2;
	dset.inputs[1][0] = 1;
	dset.inputs[1][1] = 3;

	dset.outputs[0] = 0;
	dset.outputs[1] = -2;

	return dset;
}

/* Prints the data of a solution [nodes] and [active nodes] */
void printSolution(GENES solution, int numGenes, int numInputs, int numOutputs) {
	for(int i=0; i<numGenes; i++) {
		printf("%d ", solution.nodes[i].func);
		for(int j=0; j<numInputs; j++) printf("%d ", solution.nodes[i].links[j]);
	}
	for(int i=0; i<numOutputs; i++) {
		printf("%d ", solution.nodes[numGenes+i].func);
	}
	printf("| ");
	for(int i=0; i<numGenes; i++) {
		printf("%d ", solution.toEvaluate[i]);
	}
	printf("| % .2f\n", solution.fitness);
}

NODE createNode(int index, int numRows, int levelsBack, int numFunctions, int numInputs) {
	NODE node;
	node.func = randint(0, numFunctions);
	node.links = (int *) malloc(numInputs * sizeof(int));
	int col = index / numRows;
	int min = numInputs + (col - levelsBack) * numRows;
	int max = numInputs + col * numRows;
	for(int i=0; i<numInputs; i++) {
		node.links[i] = (col >= levelsBack) ? randint(min, max) : randint(0, max);
	}
	return node;
}

NODE createOutputNode(int maxValue) {
	NODE node;
	node.func = randint(0, maxValue);
	node.links = NULL;
	return node;
}

GENES createSolution(CONFIG params) {
	GENES solution;
	int numGenes = params.numGenes;
	int numInputs = params.numInputs;
	int numOutputs = params.numOutputs;

	solution.nodes = (NODE *) malloc( (numGenes + numOutputs) * sizeof(NODE) );
	solution.toEvaluate = (int *) malloc(numGenes * sizeof(int));
	solution.fitness = -1;

	for(int i=0; i<numGenes; i++) {
		solution.nodes[i] = createNode(i, params.numRows, params.levelsBack, params.numFunctions, numInputs);
		solution.toEvaluate[i] = 0;
	}
	for(int i=0; i<numOutputs; i++) {
		solution.nodes[numGenes + i] = createOutputNode(numInputs + numGenes);
	}
	return solution;
}

void hardCodedSolution(GENES *genes, int numGenes, int numOutputs) {
	
	genes -> nodes = (NODE *) malloc( (numGenes + numOutputs) * sizeof(NODE) );
	genes -> toEvaluate = (int *) malloc( numGenes * sizeof(int));

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
	}
}

void getActiveNodes(GENES *solution, int numInputs, int numOutputs, int numGenes) {

	for(int i=0; i<numOutputs; i++) {
		if (solution -> nodes[ numGenes + i ].func >= numInputs)
			solution -> toEvaluate[ solution -> nodes[ numGenes + i ].func - numInputs ] = 1;
	}

	for(int i=numGenes-1; i>=0; i--) {
		if (solution -> toEvaluate[i]) {
			for (int j=0; j<numInputs; j++) {
				if (solution -> nodes[i].links[j] >= numInputs)
					solution -> toEvaluate[ solution -> nodes[i].links[j] - numInputs ] = 1;
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

double * executeNodes(GENES *solution, int numGenes, int numInputs, int numOutputs, double *inputData) {

	double *outputs = (double *) malloc( (numInputs + numGenes) * sizeof(double) );
	double inputs[numInputs];

	for (int i=0; i<numInputs; i++) {
		outputs[i] = inputData[i];
	}
	for (int i=0; i<numGenes; i++) {
		if (solution -> toEvaluate[i]) {
			for(int j=0; j<numInputs; j++) {
				inputs[j] = outputs[ solution -> nodes[i].links[j] ];
			}
			outputs[numInputs + i] = executeNode(solution -> nodes[i].func, inputs, numInputs);
		}
	}
	return outputs;
}

void calculateFitness(GENES *solution, int numInputs, int numOutputs, int numGenes, DATASET dataset) {
	double fitness = 0.0;

	int numCases = dataset.numCases;
	double **inputData = dataset.inputs;
	double *outputData = dataset.outputs;

	getActiveNodes(solution, numInputs, numOutputs, numGenes);

	for(int i=0; i<numCases; i++) {
		double *outputs = executeNodes(solution, numGenes, numInputs, numOutputs, inputData[i]);
		double temp = 0.0;
		for(int j=0; j<numOutputs; j++) {
			temp += outputs[ solution -> nodes[numGenes + j].func ];
		}
		fitness += pow(temp - outputData[i], 2);
	}
	solution -> fitness = fitness;
}

void freeSolution(GENES genes, int numGenes) {
	for(int i=0; i<numGenes; i++) {
		free(genes.nodes[i].links);
	}
	free(genes.nodes);
	free(genes.toEvaluate);
}

GENES * createPopulation(int popSize, CONFIG params) {
	GENES *population = (GENES *) malloc(popSize * sizeof(GENES));
	for(int i=0; i<popSize; i++) {
		population[i] = createSolution(params);
	}
	return population;
}

GENES mutation(GENES parent, CONFIG params) {

	GENES offspring;
	offspring.nodes = (NODE *) malloc( (params.numGenes + params.numOutputs) * sizeof(NODE) );
	offspring.toEvaluate = (int *) malloc( params.numGenes * sizeof(int));

	int rand1 = randint(0, params.numGenes + params.numOutputs);

	for(int i=0; i<params.numGenes+params.numOutputs; i++) {
		offspring.nodes[i] = parent.nodes[i];
		if(i < params.numGenes)
			offspring.toEvaluate[i] = 0;
	}

	if (rand1 < params.numGenes) {
		int rand2 = randint(0, params.numInputs + 1);
		if (rand2 == 0) {
			int old = offspring.nodes[rand1].func;
			int new = randint(0, params.numFunctions);
			offspring.nodes[rand1].func = new;
			printf("Old: %d New: %d\n", old, new);
		} else {
			int col = rand1 / params.numRows;
			int min = params.numInputs + (col - params.levelsBack) * params.numRows;
			int max = params.numInputs + col * params.numRows;

			int old = offspring.nodes[rand1].links[rand2];
			int new = (col >= params.levelsBack) ? randint(min, max) : randint(0, max);
			offspring.nodes[rand1].links[rand2] = new;
			printf("Old: %d New: %d\n", old, new);
		}
	} else {
		int old = offspring.nodes[rand1].func;
		int new = randint(0, params.numInputs + params.numGenes);
		offspring.nodes[rand1].func = new;
		printf("Old: %d New: %d\n", old, new);
	}
	return offspring;
}

void execute(CONFIG params) {

	GENES *population;

}

#endif