#ifndef HEADER_CGP_
#define HEADER_CGP_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include <thrust/execution_policy.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>

#include "cgp.cuh"

#define NTHREADS 1024

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
		n->inputs[i] = getRandomNodeInput(numInputs, nodePosition);
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


int getRandomNodeInput(int numChromoInputs, int nodePosition) {
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

			chromo->nodes[nodeIndex]->inputs[nodeInputIndex] = getRandomNodeInput(chromo->numInputs, nodeIndex);

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


void executeChromosome(struct chromosome *chromo, double inputs) {

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
					chromo->nodeInputsHold[j] = inputs;
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

			// if (isnan(output) != 0) output = 0;
			// else if (isinf(output) != 0) output = (output > 0) ? DBL_MAX : DBL_MIN;

			chromo->nodes[currentActiveNode]->output = output;
		}
	}

	/* Set the chromosome outputs */
	for (i = 0; i < numOutputs; i++) {
		if (chromo->outputNodes[i] < numInputs) {
			chromo->outputValues[i] = inputs;
		}
		else {
			chromo->outputValues[i] = chromo->nodes[chromo->outputNodes[i] - numInputs]->output;
		}
	}
}


double calculateFitness(struct chromosome *chromo, struct dataset *data) {
	
	int i, j;
	double error = 0.0;

	for(i = 0; i < data->numSamples; i++) {
		executeChromosome(chromo, data->inputs[i]);
		for(j = 0; j < chromo->numOutputs; j++) {
			error += fabs(chromo->outputValues[j] - data->outputs[i]);
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
	int i;
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

	dset->inputs = (double*)malloc(dset->numSamples * sizeof(double));
	dset->outputs = (double*)malloc(dset->numSamples * sizeof(double));

	for(i = 0; i < dset->numSamples; i++) {
		fgets(buffer, 100, file);
		dset->inputs[i] = atof(strtok(buffer, ","));
		dset->outputs[i] = atof(strtok(NULL, ","));
	}

	fclose(file);

	return dset;
}


void freeDataset(struct dataset *data) {
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


/* ------------------------- */
/* CUDA PART */
/* ------------------------- */

__host__ int *createArrayChromosome(struct parameters *params) {

	int i, size; int *array;

	size = params->numNodes * (params->arity + 1) + params->numOutputs;
	array = (int*)malloc(size * sizeof(int));

	for(i = 0; i < params->numNodes; i++) {
		array[i * (params->arity + 1)] = randint(0, params->numFunctions);
		array[i * (params->arity + 1) + 1] = randint(0, params->numInputs + i);
		array[i * (params->arity + 1) + 2] = randint(0, params->numInputs + i);
	}
	for(i = 0; i < params->numOutputs; i++) {
		array[params->numNodes * (params->arity + 1) + i] = randint(0, params->numInputs + params->numNodes);
	}

	return array;
}


__host__ void CUDAcreateArrayFromChromosome(struct chromosome chromo, int *array) {
	int i;
	for(i = 0; i < chromo.numNodes; i++) {
		array[i * (chromo.arity + 1)] = chromo.nodes[i]->function;
		array[i * (chromo.arity + 1) + 1] = chromo.nodes[i]->inputs[0];
		array[i * (chromo.arity + 1) + 2] = chromo.nodes[i]->inputs[1];
	}
	for(i = 0; i < chromo.numOutputs; i++) {
		array[chromo.numNodes * (chromo.arity + 1) + i] = chromo.outputNodes[i];
	}
}


__global__ void CUDAcalculateChromosomeOutputs(int *solution, double *inputs, double *outputs, int numSamples, int numInputs, int numNodes) {

	int sample = blockIdx.x * blockDim.x + threadIdx.x;

	if (sample < numSamples) {
		int i;

		double *nodeOutputs = new double[numNodes];

		for(i = 0; i < 9; i++) {

			nodeOutputs[i] = 0.0;

			int function = solution[i*3];

			int inIdx1 = solution[i*3+1];
			int inIdx2 = solution[i*3+2];

			double inValue1 = (inIdx1 < numInputs) ? inputs[sample] : nodeOutputs[inIdx1 - numInputs];
			double inValue2 = (inIdx2 < numInputs) ? inputs[sample] : nodeOutputs[inIdx2 - numInputs];

			double output = 0.0;
			switch(function) {
				case ADD:
					output = inValue1 + inValue2;
					break;
				case SUB:
					output = inValue1 - inValue2;
					break;
				case MUL:
					output = inValue1 * inValue2;
					break;
				case DIV:
					output = inValue1 / inValue2;
					break;
			}

			nodeOutputs[i] = output;
		}
		outputs[sample] = nodeOutputs[solution[27]-numInputs];//output node
		delete nodeOutputs;
	}
}


__host__ void CUDAcalculateFitness(struct chromosome *chromo, struct dataset *data) {

	int chromoSize; int *chromoArray; int *dSolPtr;

	double *dInPtr, *dOutPtr;

	chromoSize = chromo->numNodes * (chromo->arity + 1) + chromo->numOutputs;
	chromoArray = (int*)malloc(chromoSize * sizeof(int));

	CUDAcreateArrayFromChromosome(*chromo, chromoArray);

	thrust::device_vector<int> d_solution(chromoArray, chromoArray + chromoSize);

	thrust::device_vector<double> d_inputs(data->inputs, data->inputs + data->numSamples);
	thrust::device_vector<double> d_outputs(data->outputs, data->outputs + data->numSamples);

	thrust::device_vector<double> chromoOutputs(data->numSamples);

	dSolPtr = thrust::raw_pointer_cast(d_solution.data());
	dInPtr = thrust::raw_pointer_cast(d_inputs.data());
	dOutPtr = thrust::raw_pointer_cast(chromoOutputs.data());

	int numThreads = 1024;
	int numBlocks = ceil((float)data->numSamples/numThreads);
	CUDAcalculateChromosomeOutputs<<<numBlocks, numThreads>>>(dSolPtr, dInPtr, dOutPtr, data->numSamples, 1, chromo->numNodes);

	thrust::transform(d_outputs.begin(), d_outputs.end(), 
		chromoOutputs.begin(), 
		chromoOutputs.begin(), 
		thrust::minus<double>());

	chromo->fitness = thrust::reduce(chromoOutputs.begin(), chromoOutputs.end());
}


__host__ double CUDAcalculateFitness2(int *chromoArray, 
	thrust::device_vector<double> &d_inputs, 
	thrust::device_vector<double> &d_outputs, 
	int numSamples,
	int numNodes) {

	int chromoSize; 

	int *dSolPtr;
	double *dInPtr;
	double *dOutPtr;

	chromoSize = 28;//chromo->numNodes * (chromo->arity + 1) + chromo->numOutputs;
	// chromoArray = (int*)malloc(chromoSize * sizeof(int));

	thrust::device_vector<int> d_solution(chromoArray, chromoArray + chromoSize);
	thrust::device_vector<double> chromoOutputs(numSamples);

	//thrust::device_vector<double> d_inputs(data->inputs, data->inputs + data->numSamples);
	//thrust::device_vector<double> d_outputs(data->outputs, data->outputs + data->numSamples);

	dSolPtr = thrust::raw_pointer_cast(d_solution.data());
	dInPtr = thrust::raw_pointer_cast(d_inputs.data());
	dOutPtr = thrust::raw_pointer_cast(chromoOutputs.data());

	int numThreads = NTHREADS;
	int numBlocks = ceil((float)numSamples/numThreads);
	CUDAcalculateChromosomeOutputs<<<numBlocks, numThreads>>>(dSolPtr, dInPtr, dOutPtr, numSamples, 1, numNodes);

	thrust::transform(d_outputs.begin(), d_outputs.end(), 
		chromoOutputs.begin(), 
		chromoOutputs.begin(), 
		thrust::minus<double>());

	return thrust::reduce(chromoOutputs.begin(), chromoOutputs.end());
}


__host__ int *CUDAexecuteCGP(struct parameters *params, struct dataset *data, int popSize, int numGens) {

	// struct chromosome *chromo, *best;
	int *temp, *best;

	int i, j, chromoSize;

	/* device array for the chromosome */
	// thrust::device_vector<int> d_solution(chromoArray, chromoArray + chromoSize);

	/* device arrays for the data inputs and outputs */
	thrust::device_vector<double> d_inputs(data->inputs, data->inputs + data->numSamples);
	thrust::device_vector<double> d_outputs(data->outputs, data->outputs + data->numSamples);

	/* device array that holds the solution outputs */
	// thrust::device_vector<double> chromoOutputs(data->numSamples);

	/* pointers created to pass thrust arrays to kernel */
	// dSolPtr = thrust::raw_pointer_cast(d_solution.data());
	// dInPtr = thrust::raw_pointer_cast(d_inputs.data());
	// dOutPtr = thrust::raw_pointer_cast(chromoOutputs.data());

	/* creates popSize chromosomes and stores the best one */
	printf("Population\n");

	// best = createArrayChromosome(params);
	int array[] = {1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 2, 4, 1, 4, 1, 2, 2, 5, 3, 2, 3, 6};

	// calculateFitness(best, data);
	double fit = CUDAcalculateFitness2(array, d_inputs, d_outputs, data->numSamples, params->numNodes);
	printf("Fit: %6.2f\n", fit);

	// for(i = 0; i < popSize-1; i++) {
	// 	chromo = createChromosome(params);
	// 	calculateFitness(chromo, data);
	// 	if(chromo->fitness < best->fitness) {
	// 		copyChromosome(best, chromo);
	// 	} else if(chromo->fitness == best->fitness && chromo->numActiveNodes <= best->numActiveNodes) {
	// 		copyChromosome(best, chromo);
	// 	}
	// 	/* if isn't the last iteration, free it */
	// 	if (i < popSize-2) freeChromosome(chromo);
	// }

	// for(i = 0; i < numGens; i++) {
	// 	for(j = 0; j < popSize-1; j++) {
	// 		/* copies the best to mutate */
	// 		copyChromosome(chromo, best);
	// 		/* applies the mutation */
	// 		singleMutation(chromo, params);
	// 		calculateFitness(chromo, data);
	// 		/* if a mutated chromosome is better than best, save it */
	// 		if(chromo->fitness < best->fitness) {
	// 			copyChromosome(best, chromo);
	// 		} else if(chromo->fitness == best->fitness && chromo->numActiveNodes <= best->numActiveNodes) {
	// 			copyChromosome(best, chromo);
	// 		}
	// 	}
	// }
	// freeChromosome(chromo);

	return best;
}

#endif