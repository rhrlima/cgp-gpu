#ifndef HEADER_CGP_
#define HEADER_CGP_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include "cgp.h"

#define ADD 0
#define SUB 1
#define MUL 2
#define DIV 3


struct node *createNode(int numInputs, int numNodes, int arity, int numFunctions, int nodePosition) {

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
		n->inputs[i] = getRandomNodeInput(numInputs, numNodes, nodePosition);
	}

	/* set the output of the node to zero*/
	n->output = 0;

	/* set the arity of the node */
	n->maxArity = arity;
	n->actArity = arity;

	return n;
}


void copyNode(struct node *dst, struct node *src) {

	int i;

	dst->function = src->function;

	for(i = 0; i < src->maxArity; i++) {
		dst->inputs[i] = src->inputs[i];
	}

	dst->active = src->active;
	dst->output = src->output;
	dst->maxArity = src->maxArity;
	dst->actArity = src->actArity;
}


void freeNode(struct node *node) {

	if (node == NULL) {
		printf("Warning: freeing NULL node avoided\n");
		return;
	}

	free(node->inputs);
	free(node);
}


struct chromosome *createChromosome(struct parameters *params) {

	struct chromosome *chromo;
	int i;

	chromo = (struct chromosome*)malloc(sizeof(struct chromosome));

	chromo->nodes = (struct node**)malloc(params->numNodes * sizeof(struct node*));

	chromo->outputNodes = (int*)malloc(params->numOutputs * sizeof(int));

	chromo->activeNodes = (int*)malloc(params->numNodes * sizeof(int));

	chromo->outputValues = (double*)malloc(params->numOutputs * sizeof(double));

	for (i = 0; i < params->numNodes; i++) {
		chromo->nodes[i] = createNode(params->numInputs, params->numNodes, params->arity, params->numFunctions, i);
	}

	for (i = 0; i < params->numOutputs; i++) {
		chromo->outputNodes[i] = getRandomChromosomeOutput(params->numInputs, params->numNodes);
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


struct chromosome *createChromosomeFromArray(struct parameters *params, int *array) {

	struct chromosome *chromo;
	int i, j;

	chromo = (struct chromosome*)malloc(sizeof(struct chromosome));

	chromo->nodes = (struct node**)malloc(params->numNodes * sizeof(struct node*));

	chromo->outputNodes = (int*)malloc(params->numOutputs * sizeof(int));

	chromo->activeNodes = (int*)malloc(params->numNodes * sizeof(int));

	chromo->outputValues = (double*)malloc(params->numOutputs * sizeof(double));

	for (i = 0; i < params->numNodes; i++) {
		chromo->nodes[i] = (struct node*)malloc(sizeof(struct node));

		chromo->nodes[i]->inputs = (int*)malloc(params->arity * sizeof(int));

		chromo->nodes[i]->function = array[i * (params->arity + 1)];
		for (j = 0; j < params->arity; j++) {
			chromo->nodes[i]->inputs[j] = array[i * (params->arity + 1) + j + 1];
		}

		chromo->nodes[i]->active = 1;

		chromo->nodes[i]->output = 0;

		chromo->nodes[i]->maxArity = params->arity;
		chromo->nodes[i]->actArity = params->arity;
	}

	for (i = 0; i < params->numOutputs; i++) {
		chromo->outputNodes[i] = array[params->numNodes * (params->arity + 1) + i];
	}

	chromo->numInputs = params->numInputs;
	chromo->numNodes = params->numNodes;
	chromo->numOutputs = params->numOutputs;
	chromo->arity = params->arity;

	chromo->numActiveNodes = params->numNodes;

	chromo->fitness = -1;

	setChromosomeActiveNodes(chromo);

	chromo->nodeInputsHold = (double*)malloc(params->arity * sizeof(double));

	return chromo;
}


void copyChromosome(struct chromosome *dst, struct chromosome *src) {

	int i;

	dst->numInputs = src->numInputs;
	dst->numOutputs = src->numOutputs;
	dst->numNodes = src->numNodes;
	dst->numActiveNodes = src->numActiveNodes;
	dst->arity = src->arity;
	dst->fitness = src->fitness;

	for(i = 0; i < src->numNodes; i++) {
		copyNode(dst->nodes[i], src->nodes[i]);
		dst->activeNodes[i] = src->activeNodes[i];
	}

	for(i = 0; i < src->numOutputs; i++) {
		dst->outputNodes[i] = src->outputNodes[i];
		dst->outputValues[i] = src->outputValues[i];
	}

	//inputHold não precisa?
}


void freeChromosome(struct chromosome *chromo) {

	if (chromo == NULL) {
		printf("Warning: Avoiding free NULL chromosome.\n");
		return;
	}

	int i;

	for(i = 0; i < chromo->numNodes; i++) {
		freeNode(chromo->nodes[i]);
	}

	free(chromo->nodes);
	free(chromo->outputNodes);
	free(chromo->activeNodes);
	free(chromo->outputValues);
	free(chromo->nodeInputsHold);
	free(chromo);
}


int getRandomNodeInput(int numChromoInputs, int numNodes, int nodePosition) {
	/* pick any previous node including inputs */
	return randint(0, numChromoInputs + nodePosition);
}


int getRandomChromosomeOutput(int numInputs, int numNodes) {

	/* returns any previous node */
	return randint(0, numInputs + numNodes);
}


void setChromosomeActiveNodes(struct chromosome *chromo) {

	int i, j;

	chromo->numActiveNodes = 0;

	for(i = 0; i < chromo->numNodes; i++) {
		chromo->activeNodes[i] = 0;
	}

	for(i = 0; i < chromo->numOutputs; i++) {
		if (chromo->outputNodes[i] >= chromo->numInputs)
			chromo -> activeNodes[ chromo->outputNodes[i] - chromo->numInputs ] = 1;
	}

	for(i = chromo->numNodes-1; i >= 0; i--) {
		if (chromo->activeNodes[i]) {
			for (j=0; j<chromo->arity; j++) {
				if (chromo->nodes[i]->inputs[j] >= chromo->numInputs)
					chromo->activeNodes[ chromo->nodes[i]->inputs[j] - chromo->numInputs ] = 1;
			}
			chromo->numActiveNodes++;
		}
	}
}


void printChromosome(struct chromosome *chromo) {
	int i, j;

	if (chromo == NULL) {
		printf("Warning: Avoiding print NULL chromosome.\n");
		return;
	}

	for(i = 0; i < chromo->numNodes; i++) {

		printf("%d ", chromo->nodes[i]->function);

		for(j = 0; j < chromo->arity; j++) {

			printf("%d ", chromo->nodes[i]->inputs[j]);
		}
	}
	printf("| ");
	for(i = 0; i < chromo->numOutputs; i++) {

		printf("%d ", chromo->outputNodes[i]);
	}
	printf("| ");
	for(i = 0; i < chromo->numNodes; i++){

		printf("%d ", chromo->activeNodes[i]);
	}
	printf("| %.2f\n", chromo->fitness);
}


/* -------------------------------------------------- */


void singleMutation(struct chromosome *chromo, struct parameters *params) {

	int numFunctionGenes, numInputGenes, numOutputGenes;
	int numGenes;
	int geneToMutate;
	int nodeIndex;
	int nodeInputIndex;

	int mutatedActive = 0;
	int previousGeneValue;
	int newGeneValue;

	/* get the number of each type of gene */
	numFunctionGenes = params->numNodes;
	numInputGenes = params->numNodes * params->arity;
	numOutputGenes = params->numOutputs;

	/* set the total number of chromosome genes */
	numGenes = numFunctionGenes + numInputGenes + numOutputGenes;

	/* while active gene not mutated */
	while (mutatedActive == 0) {

		/* select a random gene */
		geneToMutate = randint(0, numGenes);

		/* mutate function gene */
		if (geneToMutate < numFunctionGenes) {

			nodeIndex = geneToMutate;

			previousGeneValue = chromo->nodes[nodeIndex]->function;

			chromo->nodes[nodeIndex]->function = randint(0, params->numFunctions);

			newGeneValue = chromo->nodes[nodeIndex]->function;

			if ((previousGeneValue != newGeneValue) && (chromo->nodes[nodeIndex]->active == 1)) {
				mutatedActive = 1;
			}

		}

		/* mutate node input gene */
		else if (geneToMutate < numFunctionGenes + numInputGenes) {

			nodeIndex = (int) ((geneToMutate - numFunctionGenes) / chromo->arity);
			nodeInputIndex = (geneToMutate - numFunctionGenes) % chromo->arity;

			previousGeneValue = chromo->nodes[nodeIndex]->inputs[nodeInputIndex];

			chromo->nodes[nodeIndex]->inputs[nodeInputIndex] = getRandomNodeInput(chromo->numInputs, chromo->numNodes, nodeIndex);

			newGeneValue = chromo->nodes[nodeIndex]->inputs[nodeInputIndex];

			if ((previousGeneValue != newGeneValue) && (chromo->nodes[nodeIndex]->active == 1)) {
				mutatedActive = 1;
			}
		}

		/* mutate output gene */
		else {
			nodeIndex = geneToMutate - numFunctionGenes - numInputGenes;

			previousGeneValue = chromo->outputNodes[nodeIndex];

			chromo->outputNodes[nodeIndex] = getRandomChromosomeOutput(chromo->numInputs, chromo->numNodes);

			newGeneValue = chromo->outputNodes[nodeIndex];

			if (previousGeneValue != newGeneValue) {
				mutatedActive = 1;
			}
		}
	}
}


/* -------------------------------------------------- */


void executeChromosome(struct chromosome *chromo, double *inputs) {

	int i, j;

	int nodeInputLocation;
	int currentActiveNode;
	int currentActiveNodeFuction;
	int nodeArity;

	int numInputs = chromo->numInputs;
	int numNodes = chromo->numNodes;
	int numOutputs = chromo->numOutputs;

	/* error checking */
	if (chromo == NULL) {
		printf("Error: cannot execute uninitialised chromosome.\n Terminating CGP-Library.\n");
		exit(0);
	}

	/* for all of the active nodes */
	for (i = 0; i < numNodes; i++) {

		if (chromo->activeNodes[i]) {

			currentActiveNode = i;
			nodeArity = chromo->nodes[currentActiveNode]->actArity;

			for(j = 0; j < nodeArity; j++) {
				
				nodeInputLocation = chromo->nodes[currentActiveNode]->inputs[j];

				/* verify if the input location is a node or real input */
				if(nodeInputLocation < numInputs) {
					chromo->nodeInputsHold[j] = inputs[nodeInputLocation];
				}
				else {
					chromo->nodeInputsHold[j] = chromo->nodes[nodeInputLocation - numInputs]->output;
				}
			}

			/* get the functionality of the active node under evaluation */
			currentActiveNodeFuction = chromo->nodes[currentActiveNode]->function;

			/* calculate the output of the active node under evaluation */
			//melhorar depois
			double output = 0.0;
			switch(currentActiveNodeFuction) {
				case ADD:
					output = chromo->nodeInputsHold[0] + chromo->nodeInputsHold[1];
					break;
				case SUB:
					output = chromo->nodeInputsHold[0] - chromo->nodeInputsHold[1];
					break;
				case MUL:
					output = chromo->nodeInputsHold[0] * chromo->nodeInputsHold[1];
					break;
				case DIV:
					output = chromo->nodeInputsHold[0] / chromo->nodeInputsHold[1];
					break;
			}

			if (isnan(output) != 0) output = 0;
			else if (isinf(output) != 0) output = (output > 0) ? DBL_MAX : DBL_MIN;

			chromo->nodes[currentActiveNode]->output = output;
		}
	}

	/* Set the chromosome outputs */
	for (i = 0; i < numOutputs; i++) {

		if (chromo->outputNodes[i] < numInputs) {
			chromo->outputValues[i] = inputs[chromo->outputNodes[i]];
		}
		else {
			chromo->outputValues[i] = chromo->nodes[chromo->outputNodes[i] - numInputs]->output;
		}
		//printf("output: %.2f\n", chromo->outputValues[i]);
	}
}


double calculateFitness(struct chromosome *chromo, struct dataset *data) {
	
	int i, j;
	double error = 0;

	for(i = 0; i < data->numSamples; i++) {

		executeChromosome(chromo, data->inputs[i]);

		for(j = 0; j < chromo->numOutputs; j++) {
			error += fabs(chromo->outputValues[j] - data->outputs[i][j]);
		}

	}

	chromo->fitness = error;
	return error;
}


/* -------------------------------------------------- */


struct chromosome *executeCGP(struct parameters *params, struct dataset *data, int numGens) {

	struct chromosome *chromo, *best;

	int i, j;
	int popSize = 5;

	/* creates popSize chromosomes and stores the best one */
	printf("Population\n");

	best = createChromosome(params);
	calculateFitness(best, data);

	for(i = 0; i < popSize-1; i++) {
		chromo = createChromosome(params);
		calculateFitness(chromo, data);
		if(chromo->fitness < best->fitness) {
			copyChromosome(best, chromo);
		} else if(chromo->fitness == best->fitness && chromo->numActiveNodes <= best->numActiveNodes) {
			copyChromosome(best, chromo);
		}
		/* if isn't the last iteration, free it */
		if (i < popSize-2) freeChromosome(chromo);
	}

	for(i = 0; i < numGens; i++) {
		for(j = 0; j < popSize-1; j++) {
			/* copies the best to mutate */
			copyChromosome(chromo, best);
			/* applies the mutation */
			singleMutation(chromo, params);
			calculateFitness(chromo, data);
			/* if a mutated chromosome is better than best, save it */
			if(chromo->fitness < best->fitness) {
				copyChromosome(best, chromo);
			} else if(chromo->fitness == best->fitness && chromo->numActiveNodes <= best->numActiveNodes) {
				copyChromosome(best, chromo);
			}
		}
	}
	freeChromosome(chromo);

	return best;
}


/* -------------------------------------------------- */


struct dataset *loadDataset(char *fileName) {

	FILE *file;
	int i, j;
	struct dataset *dset;
	char buffer[100];

	file = fopen(fileName, "r");
	if (!file) {
		printf("Error: file %s was not found.\nExiting.\n", fileName);
		exit(0);
	}

	fgets(buffer, 100, file);
	
	dset = (struct dataset *)malloc(sizeof(struct dataset));

	dset->numInputs  = atoi(strtok(buffer, ","));
	dset->numOutputs = atoi(strtok(NULL, ","));
	dset->numSamples = atoi(strtok(NULL, ","));

	dset->inputs = (double **) malloc(dset->numSamples * sizeof(double*));
	for(i = 0; i < dset->numSamples; i++) dset->inputs[i] = (double *) malloc(dset->numInputs * sizeof(double));

	dset->outputs = (double **) malloc(dset->numSamples * sizeof(double*));
	for(i = 0; i < dset->numSamples; i++) dset->outputs[i] = (double *) malloc(dset->numOutputs * sizeof(double));

	for(i = 0; i < dset->numSamples; i++) {

		fgets(buffer, 100, file);

		dset->inputs[i][0] = atof(strtok(buffer, ","));
		for(j = 1; j < dset->numInputs; j++) {
			dset->inputs[i][j] = atof(strtok(NULL, ","));
		}
		for(int j=0; j<dset->numOutputs; j++) {
			dset->outputs[i][j] = atof(strtok(NULL, ","));
		}
	}

	fclose(file);

	return dset;
}


void freeDataset(struct dataset *data) {
	int i, j;
	for(i = 0; i < data->numSamples; i++) {
		free(data->inputs[i]);
		free(data->outputs[i]);
	}
	free(data->inputs);
	free(data->outputs);
	free(data);
}


/* -------------------------------------------------- */


int randint(int min, int max) {

	return rand() % (max-min) + min;
}


float randfloat(float min, float max) {

	return (rand() / (float) RAND_MAX) * (max - min) + min;
}


struct parameters *initialiseParameters(int numNodes, int arity, int numFunctions, struct dataset *data) {
	
	struct parameters *params;

	if (data == NULL) {
		printf("Error: Dataset not initialised.\nExiting.\n");
		exit(0);
	}

	params = (struct parameters*)malloc(sizeof(struct parameters));

	params->numInputs = data->numInputs;
	params->numOutputs = data->numOutputs;
	params->numNodes = numNodes;
	params->arity = arity;
	params->numFunctions = numFunctions;

	return params;
}


void printParameters(struct parameters *params) {
	printf("Inputs: %d\n", params->numInputs);
	printf("Outputs: %d\n", params->numOutputs);
	printf("Nodes: %d\n", params->numNodes);
	printf("Max Arity: %d\n", params->arity);
	printf("Functions: %d\n", params->numFunctions);
}

#endif