#ifndef HEADER_CGP_
#define HEADER_CGP_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "cgp.h"

#define ADD 0
#define SUB 1
#define MUL 2
#define DIV 3

float sum(float a, float b) {
	return a + b;
}

/* Random int between [min, max] exclusive */
int randint(int min, int max) {

	return rand() % (max-min) + min;
}

/* Random float between [min, max] exclusive */
float randfloat(float min, float max) {

	return (rand() / (float) RAND_MAX) * (max - min) + min;
}

DATASET loadDataset(char *fileName) {
	FILE *file;
	DATASET dset;
	char line[100]; //line buffer

	file = fopen(fileName, "r");
	if (!file) {
		printf("Erro ao tentar abrir arquivo: %s\n", fileName);
	}

	fgets(line, 100, file);
	
	//First line has always 3 values [numInputs, numOutputs, numCases]
	dset.numInputs 	= atoi(strtok(line, ","));
	dset.numOutputs = atoi(strtok(NULL, ","));
	dset.numCases 	= atoi(strtok(NULL, ","));

	dset.inputs = (double **) malloc(dset.numCases * sizeof(double));
	for(int i=0; i<dset.numCases; i++) dset.inputs[i] = (double *) malloc(dset.numInputs * sizeof(double));

	dset.outputs = (double **) malloc(dset.numCases * sizeof(double));
	for(int i=0; i<dset.numCases; i++) dset.outputs[i] = (double *) malloc(dset.numOutputs * sizeof(double));

	for(int i=0; i<dset.numCases; i++) {

		char *r = fgets(line, 256, file);

		//FIXING NUM INPUTS AND OUTPUTS AS 1
		r = strtok(line, ",");
		dset.inputs[i][0] = atof(r);

		r = strtok(NULL, ",");
		dset.outputs[i][0] = atof(r);

		//----

		// dset.inputs[i][0] = atof(strtok(line, ","));
		// printf("%d ", dset.inputs[i][0]);
		// for(int j=1; j<dset.numInputs; j++) {
		// 	dset.inputs[i][j] = atof(strtok(NULL, ","));
		// 	printf("%d ", dset.inputs[i][j]);
		// }
		// printf("\n");
		// for(int j=0; j<dset.numOutputs; j++) {
		// 	dset.outputs[i][j] = atof(strtok(NULL, ","));
		// 	printf("%d ", dset.outputs[i][j]);
		// }
		// printf("\n");
	}

	return dset;
}

void freeDataset(DATASET dataset) {
	for(int i=0; i<dataset.numCases; i++) {
		free(dataset.inputs[i]);
		free(dataset.outputs[i]);
	}
	free(dataset.inputs);
	free(dataset.outputs);
}

/* Prints the data of a solution [nodes] and [active nodes] */
void printSolution(struct chromosome *chromo) {
	int i, j;

	for(i = 0; i < chromo->numNodes; i++) {

		printf("%d ", chromo->nodes[i]->function);
		for(j = 0; j < chromo->arity; j++) printf("%d ", chromo->nodes[i]->links[j]);

	}
	// for(int i=0; i<numOutputs; i++) {
	// 	printf("%d ", solution.nodes[numGenes+i].func);
	// }
	// printf("| ");
	// for(int i=0; i<numGenes; i++) {
	// 	printf("%d ", solution.toEvaluate[i]);
	// }
	// printf("| % .2f\n", solution.fitness);
}

static int getRandomNodeInput(int numChromoInputs, int numNodes, int nodePosition) {
	/* pick any previous node including inputs */
	return randint(0, numChromoInputs + nodePosition);
}

static struct node *createNode(int numInputs, int numNodes, int arity, int numFunctions, int nodePosition) {

	struct node *n;
	int i;

	/* allocate memory for node */
	n = (struct node*)malloc(sizeof(struct node));

	/* allocate memory for the node's inputs and connection weights */
	n->inputs = (int*)malloc(arity * sizeof(int));

	/* set the node's function */
	n->function = randint(0, numFunctions);

	/* set as active by default */
	n->active = 1;

	/* set the nodes inputs and connection weights */
	for (i = 0; i < arity; i++) {
		n->inputs[i] = getRandomNodeInput(numInputs, numNodes, nodePosition, levelsBack);
	}

	/* set the output of the node to zero*/
	n->output = 0;

	/* set the arity of the node */
	n->maxArity = arity;

	return n;
}

static int getRandomChromosomeOutput(int numInputs, int numNodes) {

	/* returns any previous node */
	return randint(0, numInputs + numNodes);
}

/* Creates a solution (GENES) */
struct chromosome *createChromosome(struct parameters *params) {

	struct chromosome *chromo;
	int i;

	chromo = (struct chromosome*)malloc(sizeof(struct chromosome));

	chromo->nodes = (struct node**)malloc(params->numNodes * sizeof(struct node*));

	chromo->outputNodes = (int*)malloc(params->numOutputs * sizeof(int));

	chromo->activeNodes = (int*)malloc(params->numNodes * sizeof(int));

	chromo->outputValues = (double*)malloc(params->numOutputs * sizeof(double));

	for (i = 0; i < params->numNodes; i++) {
		chromo->nodes[i] = createNode(paras->numInputs, paras->numNodes, paras->arity, paras->numFunctions, i);
	}

	for (i = 0; i < params->numOutputs; i++) {
		chromo->outputNodes[i] = createOutputNode(paras->numInputs, paras->numNodes);
	}

	chromo->numInputs = params->numInputs;
	chromo->numNodes = params->numNodes;
	chromo->numOutputs = params->numOutputs;
	chromo->arity = params->arity;

	/* set the number of active node to the number of nodes (all active) */
	chromo->numActiveNodes = params->numNodes;

	/* set the fitness to initial value */
	chromo->fitness = -1;

	/* set the active nodes in the newly generated chromosome */
	setChromosomeActiveNodes(chromo);

	/* used interally when executing chromosome */
	chromo->nodeInputsHold = (double*)malloc(params->arity * sizeof(double));

	return chromo;
}

struct chromosome * createSolutionFromArray(int *vet, int numGenes, int arity, int numOutputs) {
	GENES solution;

	solution.nodes = (NODE *) malloc( (numGenes + numOutputs) * sizeof(NODE) );
	solution.toEvaluate = (int *) malloc( numGenes * sizeof(int));
	solution.fitness = -1;	

	for(int i=0; i<numGenes; i++) {
		NODE aux;
		aux.links = (int *) malloc(arity * sizeof(int));
		aux.func = vet[i * (arity+1)];
		for(int j=0; j<arity; j++){
			aux.links[j] = vet[i*(arity+1)+j+1];
		}
		solution.nodes[i] = aux;
		solution.toEvaluate[i] = 0;
	}

	for(int i=0; i<numOutputs; i++) {
		NODE aux;
		aux.func = vet[numGenes*(arity+1)+i];
		aux.links = NULL;
		solution.nodes[numGenes+i] = aux;
	}

	return solution;
}

void getActiveNodes(struct chromosome *chromo) {

	for(int i=0; i<numOutputs; i++) {
		printf("%d ", solution -> nodes[ numGenes + i ].func);
		if (solution -> nodes[ numGenes + i ].func >= arity)
			solution -> toEvaluate[ solution -> nodes[ numGenes + i ].func - arity ] = 1;
	}
	printf("\n");

	for(int i=numGenes-1; i>=0; i--) {
		if (solution -> toEvaluate[i]) {
			for (int j=0; j<arity; j++) {
				if (solution -> nodes[i].links[j] >= arity)
					solution -> toEvaluate[ solution -> nodes[i].links[j] - arity ] = 1;
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
		//printf("%.2f\n", inputData[i]);
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

void calculateFitness(GENES *solution, int arity, int numInputs, int numOutputs, int numGenes, DATASET dataset) {
	double fitness = 0.0;

	int numCases = dataset.numCases;
	double **inputData = dataset.inputs;
	double **outputData = dataset.outputs;

	getActiveNodes(solution, arity, numOutputs, numGenes);

	for(int i=0; i<numCases; i++) {
		double *outputs = executeNodes(solution, numGenes, numInputs, numOutputs, inputData[i]);
		double solutionTemp = 0.0;
		double outputTemp = 0.0;
		for(int j=0; j<numOutputs; j++) {
			solutionTemp += outputs[ solution -> nodes[numGenes + j].func ];
			outputTemp += outputData[i][j];
		}
		fitness += pow(solutionTemp - outputTemp, 2);
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

/* Creates a population of solutions (GENES) */
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

	GENES *population, best, aux;

	int gens = 0;
	int numInputs 	= params.numInputs;
	int numOutputs 	= params.numOutputs;
	int numGenes 	= params.numGenes;

	population = createPopulation(popSize, params);

	best = copySolution(population[0], numGenes, numInputs, numOutputs);
	calculateFitness(&best, numInputs, numOutputs, numGenes, dataset);

	for(int i=0; i<popSize; i++) {
		calculateFitness(&population[i], numInputs, numOutputs, numGenes, dataset);
		//printSolution(population[i], numGenes, numInputs, numOutputs);
		if (population[i].fitness <= best.fitness)
			best = copySolution(population[i], numGenes, numInputs, numOutputs);
		gens += popSize;
	}

	for(int i=0; i<numGen; i++) {
		for(int j=0; j<popSize-1; j++) {
			printf("A");
			population[j] = mutation(best, params);
			printf("B");
			printf("\n");
			printSolution(population[j], numGenes, numInputs, numOutputs);
			calculateFitness(&population[j], numInputs, numOutputs, numGenes, dataset);
			printf("C\n");
			if(population[j].fitness <= best.fitness) {
				aux = copySolution(population[j], numGenes, numInputs, numOutputs);
			}
		}
		best = copySolution(aux, numGenes, numInputs, numOutputs);
		//printSolution(best, numGenes, numInputs, numOutputs);
		if(best.fitness == 0.0) {break;}
		gens += popSize-1;
	}

	for(int i=0; i<popSize; i++) {
		freeSolution(population[i], numGenes);
	}

	printSolution(best, numGenes, numInputs, numOutputs);
}

#endif