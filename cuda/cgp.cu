#ifndef HEADER_CGP_
#define HEADER_CGP_

#include "cgp.cuh"

#define NUMNODES 100
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

			if (isnan(output) != 0) output = 0;
			else if (isinf(output) != 0) output = (output > 0) ? DBL_MAX : DBL_MIN;

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
/*         CUDA PART         */
/* ------------------------- */


struct absminus {
	__host__ __device__
	double operator()(double &a, double &b) const {
		return fabs(a - b);
	}
};


__host__ void CUDAcreateArrayChromosome(thrust::host_vector<int> &array, struct parameters *params) {
	int i;
	for(i = 0; i < params->numNodes; i++) {
		array[i * (params->arity + 1)] = randint(0, params->numFunctions);
		array[i * (params->arity + 1) + 1] = randint(0, params->numInputs + i);
		array[i * (params->arity + 1) + 2] = randint(0, params->numInputs + i);
	}
	for(i = 0; i < params->numOutputs; i++) {
		array[params->numNodes * (params->arity + 1) + i] = randint(0, params->numInputs + params->numNodes);
	}
}


__global__ void CUDAcalculateChromosomeOutputs(int *solution, double *inputs, double *outputs, int numSamples, int numInputs, int numNodes) {

	int sample = blockIdx.x * blockDim.x + threadIdx.x;

	if (sample < numSamples) {
		int i;
		
		double nodeOutputs[NUMNODES];

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

			if (isnan(output) != 0) output = 0;
			else if (isinf(output) != 0) output = (output > 0) ? DBL_MAX : DBL_MIN;

			nodeOutputs[i] = output;
		}
		outputs[sample] = nodeOutputs[solution[27]-numInputs];
	}
}


__host__ double CUDAcalculateFitness(thrust::host_vector<int> &h_solution, thrust::device_vector<double> &d_inputs, thrust::device_vector<double> &d_outputs, int numSamples, int numInputs, int numNodes) {

	int *dSolPtr;
	double *dInPtr;
	double *dOutPtr;

	thrust::device_vector<int> d_solution(h_solution.begin(), h_solution.end());
	thrust::device_vector<double> chromoOutputs(numSamples);

	dSolPtr = thrust::raw_pointer_cast(d_solution.data());
	dInPtr  = thrust::raw_pointer_cast(d_inputs.data());
	dOutPtr = thrust::raw_pointer_cast(chromoOutputs.data());

	int numThreads = NTHREADS;
	int numBlocks = ceil((float)numSamples/numThreads);
	CUDAcalculateChromosomeOutputs<<<numBlocks, numThreads>>>(dSolPtr, dInPtr, dOutPtr, numSamples, numInputs, numNodes);

	thrust::transform(d_outputs.begin(), d_outputs.end(), 
		chromoOutputs.begin(), 
		chromoOutputs.begin(),
		absminus());

	return thrust::reduce(chromoOutputs.begin(), chromoOutputs.end());
}


__host__ void CUDAsingleMutation(thrust::host_vector<int> &solution, struct parameters *params) {
	
	int size = params->numNodes * (params->arity + 1) + params->numOutputs;
	int geneToMutate, nodeIndex, subIndex, oldValue, newValue;

	do {
		/* picks a random gene to mutate */
		geneToMutate = randint(0, size);

		nodeIndex = geneToMutate/(params->arity+1);
		subIndex = geneToMutate%(params->arity+1);

		/*store the old value */
		oldValue = solution[geneToMutate];

		/* mutate normal node */
		if (nodeIndex < params->numNodes) {
			/* mutate function gene */
			if (subIndex == 0) {
				newValue = randint(0, params->numFunctions);
			}
			/* mutate input gene */
			else {
				newValue = getRandomNodeInput(params->numInputs, nodeIndex);
			}
		}
		/* mutate output node */
		else {
			newValue = getRandomChromosomeOutput(params->numInputs, params->numNodes);
		}
	} while(oldValue == newValue);

	solution[geneToMutate] = newValue;
}


__host__ int *CUDAexecuteCGP(struct parameters *params, struct dataset *data, int popSize, int numGens) {

	int i, j, solSize; int *result;
	double bestFit, tempFit;

	solSize = params->numNodes * (params->arity + 1) + params->numOutputs;

	/* host arrays that handles solutions */
	thrust::host_vector<int> best(solSize);
	thrust::host_vector<int> temp(solSize);

	/* device arrays for the data inputs and outputs */
	thrust::device_vector<double> d_inputs(data->inputs, data->inputs + data->numSamples);
	thrust::device_vector<double> d_outputs(data->outputs, data->outputs + data->numSamples);

	/* creates popSize chromosomes and stores the best one */
	CUDAcreateArrayChromosome(best, params);

	bestFit = CUDAcalculateFitness(best, d_inputs, d_outputs, data->numSamples, params->numInputs, params->numNodes);

	for(i = 0; i < popSize-1; i++) {

		CUDAcreateArrayChromosome(temp, params);

		tempFit = CUDAcalculateFitness(temp, d_inputs, d_outputs, data->numSamples, params->numInputs, params->numNodes);

		if(tempFit <= bestFit) {
			/* copying temp --> best */
			thrust::copy(temp.begin(), temp.end(), best.begin());
			bestFit = tempFit;
		}
	}

	for(i = 0; i < numGens; i++) {

		for(j = 0; j < popSize-1; j++) {

			/* copies the best to mutate best --> temp */
			thrust::copy(best.begin(), best.end(), temp.begin());

			/* applies the mutation */
			CUDAsingleMutation(temp, params);
			tempFit = CUDAcalculateFitness(temp, d_inputs, d_outputs, data->numSamples, params->numInputs, params->numNodes);
			
			/* if a mutated chromosome is better than best, save it  */
			if(tempFit <= bestFit) {
				/* copying temp --> best */
				thrust::copy(temp.begin(), temp.end(), best.begin());
				bestFit = tempFit;
			}
		}
	}

	result = (int*)malloc(solSize * sizeof(int));
	thrust::copy(best.begin(), best.end(), result);
	return result;
}

#endif