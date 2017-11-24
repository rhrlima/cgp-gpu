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

void freeDataset(DATASET dataset) {
	for(int i=0; i<dataset.numCases; i++)
		free(dataset.inputs[i]);
	free(dataset.inputs);
	free(dataset.outputs);
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

GENES hardCodedSolution() {
	GENES solution;

	solution.nodes = (NODE *) malloc( 10 * sizeof(NODE) );
	solution.toEvaluate = (int *) malloc( 6 * sizeof(int));
	solution.fitness = -1;

	int vet[22] = {0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 2, 5, 7, 3};
	for(int i=0; i<6; i++) {
		NODE aux;
		aux.links = (int *) malloc(2 * sizeof(int));
		aux.func = vet[i*3];
		aux.links[0] = vet[i*3+1];
		aux.links[1] = vet[i*3+2];

		solution.nodes[i] = aux;
		solution.toEvaluate[i] = 0;
	}

	for(int i=0; i<4; i++) {
		NODE aux;
		aux.func = vet[6*3+i];
		aux.links = NULL;
		solution.nodes[6+i] = aux;
	}

	return solution;
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

NODE copyNode(NODE source, int numInputs) {
	NODE copy;
	copy.func = source.func;
	copy.links = NULL;
	if(source.links != NULL) {
		copy.links = (int *) malloc(numInputs * sizeof(int));
		for(int i=0; i<numInputs; i++) copy.links[i] = source.links[i];
	}
	return copy;
}

GENES copySolution(GENES source, int numGenes, int numInputs, int numOutputs) {
	GENES copy;
	copy.nodes = (NODE *) malloc( (numGenes + numOutputs) * sizeof(NODE) );
	copy.toEvaluate = (int *) malloc(numGenes * sizeof(int));
	copy.fitness = source.fitness;
	for(int i=0; i<numGenes+numOutputs; i++) {
		copy.nodes[i] = copyNode(source.nodes[i], numInputs);
		if(i < numGenes) copy.toEvaluate[i] = source.toEvaluate[i];
	}
	return copy;
}

void freeSolution(GENES solution, int numGenes) {
	for(int i=0; i<numGenes; i++) {
		free(solution.nodes[i].links);
	}
	free(solution.nodes);
	free(solution.toEvaluate);
}

GENES * createPopulation(int popSize, CONFIG params) {
	GENES *population = (GENES *) malloc(popSize * sizeof(GENES));
	for(int i=0; i<popSize; i++) {
		population[i] = createSolution(params);
	}
	return population;
}

GENES mutation(GENES parent, CONFIG params) {

	int numInputs 	= params.numInputs;
	int numOutputs 	= params.numOutputs;
	int numGenes 	= params.numGenes;

	GENES offspring = copySolution(parent, numGenes, numInputs, numOutputs);

	int rand1 = randint(0, numGenes + numOutputs);

	if (rand1 < numGenes) {
		int rand2 = randint(0, numInputs + 1);
		if (rand2 == 0) {
			int old = offspring.nodes[rand1].func;
			int new = randint(0, params.numFunctions);
			offspring.nodes[rand1].func = new;
			//printf("Old: %d New: %d\n", old, new);

		} else {
			int col = rand1 / params.numRows;
			int min = numInputs + (col - params.levelsBack) * params.numRows;
			int max = numInputs + col * params.numRows;

			int old = offspring.nodes[rand1].links[rand2];
			int new = (col >= params.levelsBack) ? randint(min, max) : randint(0, max);
			offspring.nodes[rand1].links[rand2] = new;
			//printf("Old: %d New: %d\n", old, new);

		}
	} else {
		int old = offspring.nodes[rand1].func;
		int new = randint(0, numInputs + numGenes);
		offspring.nodes[rand1].func = new;
		//printf("Old: %d New: %d\n", old, new);

	}
	return offspring;
}

void execute(int popSize, int numGen, CONFIG params, DATASET dataset) {

	GENES *population, best;

	double mtr = 0.5;

	int numInputs 	= params.numInputs;
	int numOutputs 	= params.numOutputs;
	int numGenes 	= params.numGenes;

	population = createPopulation(popSize, params);

	best = copySolution(population[0], numGenes, numInputs, numOutputs);
	calculateFitness(&best, numInputs, numOutputs, numGenes, dataset);

	for(int i=0; i<popSize; i++) {
		calculateFitness(&population[i], numInputs, numOutputs, numGenes, dataset);
		printSolution(population[i], 6, 2, 4);
		if (population[i].fitness < best.fitness)
			best = copySolution(population[i], numGenes, numInputs, numOutputs);
	}
	for(int i=0; i<numGen; i++) {
		for(int j=0; j<popSize; j++) {
			if (randfloat(0, 1) < mtr) {
				//printf("Mutated\n");
				GENES offspring = mutation(population[j], params);
				calculateFitness(&offspring, numInputs, numOutputs, numGenes, dataset);
				if (offspring.fitness < population[j].fitness)
					population[j] = copySolution(offspring, numGenes, numInputs, numOutputs);
			}
		}
		for(int j=0; j<popSize; j++) {
			if (population[j].fitness < best.fitness)
				best = copySolution(population[j], numGenes, numInputs, numOutputs);
		}
		//printf("Gen: % 4d | Best: % 4.2f\n", i, best.fitness);
	
	}

	for(int i=0; i<popSize; i++) {
		freeSolution(population[i], numGenes);
	}

	printSolution(best, 6, 2, 4);
}

#endif