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

int randint(int min, int max) {

	return rand() % (max-min) + min;
}


float randfloat(float min, float max) {

	return (rand() / (float) RAND_MAX) * (max - min) + min;
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
		chromo->outputNodes[i] = getRandomNodeOutput(params->numInputs, params->numNodes);
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


int getRandomNodeInput(int numChromoInputs, int numNodes, int nodePosition) {
	/* pick any previous node including inputs */
	return randint(0, numChromoInputs + nodePosition);
}


int getRandomNodeOutput(int numInputs, int numNodes) {

	/* returns any previous node */
	return randint(0, numInputs + numNodes);
}


void setChromosomeActiveNodes(struct chromosome *chromo) {

	int i, j;

	for(i = 0; i < chromo->numNodes; i++) {
		chromo->activeNodes[i] = 0;
	}

	for(i = 0; i < chromo->numOutputs; i++) {
		chromo -> activeNodes[ chromo->outputNodes[i] - chromo->numInputs ] = 1;
	}

	for(i = chromo->numNodes-1; i >= 0; i--) {

		if (chromo->activeNodes[i]) {

			for (j=0; j<chromo->arity; j++) {
				if (chromo->nodes[i]->inputs[j] >= chromo->numInputs)
					chromo->activeNodes[ chromo->nodes[i]->inputs[j] - chromo->numInputs ] = 1;
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
		printf("output: %d\n", chromo->outputValues[i]);
	}
}


double calculateFitness(struct chromosome *chromo, struct dataset *data) {
	
	int i, j;
	double error = 0;

	for(i = 0; i < data.numSamples; i++) {

		executeChromosome(chromo, data.inputs[i]);

		for(j = 0; j < chromo->numOutputs; j++) {
			error += fabs(chromo->outputValues[j] - data.outputs[i][j]);
		}

	}

	chromo->fitness = error;
	return error;
}


/* -------------------------------------------------- */


struct dataset *loadDataset(char *fileName) {
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


void printChromosome(struct chromosome *chromo) {
	int i, j;

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

#endif