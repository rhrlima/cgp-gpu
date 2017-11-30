#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/reduce.h>
#include <thrust/transform.h>

struct node {
	int function;
	int *inputs;
	int active;
	double output;
	int maxArity;
	int actArity;
};

struct chromosome {
	int numInputs;
	int numOutputs;
	int numNodes;
	int numActiveNodes;
	int arity;
	struct node **nodes;
	int *outputNodes;
	int *activeNodes;
	double fitness;
	double *outputValues;
	double *nodeInputsHold;
};

struct parameters {
	int numInputs;		//Inputs for the dataset
	int numOutputs;		//Outputs for the dataset
	int numNodes;		//Cols * Rows
	int arity;			//Max arity from functions (usually 2)
	int numFunctions;	//Available functions
};

struct dataset {
	int numInputs;
	int numOutputs;
	int numSamples;
	double *inputs;
	double *outputs;
};


/* chromosome creation */
struct node *createNode(int numInputs, int numNodes, int arity, int numFunctions, int nodePosition);

void copyNode(struct node *dst, struct node *src);

void freeNode(struct node *node);

struct chromosome *createChromosome(struct parameters *params);

struct chromosome *createChromosomeFromArray(struct parameters *params, int *array);

int getRandomNodeInput(int numChromoInputs, int nodePosition);

int getRandomChromosomeOutput(int numInputs, int numNodes);

void setChromosomeActiveNodes(struct chromosome *chromo);

void copyChromosome(struct chromosome *dst, struct chromosome *src);

void freeChromosome(struct chromosome *chromo);

void printChromosome(struct chromosome *chromo);


/* mutation */
void singleMutation(struct chromosome *chromo, struct parameters *params);


/* chromosome evaluation */
void executeChromosome(struct chromosome *chromo, double inputs);

double calculateFitness(struct chromosome *chromo, struct dataset *data);


/* evolutionary strategy */
struct chromosome *executeCGP(struct parameters *params, struct dataset *data, int numGens);


/* dataset */
struct dataset *loadDataset(char *fileName);

void freeDataset(struct dataset *data);


/* utils */
int randint(int min, int max);

float randfloat(float min, float max);

struct parameters *initialiseParameters(int numNodes, int arity, int numFunctions, struct dataset *data);

void printParameters(struct parameters *params);


/* ------------------------- */
/*         CUDA PART         */
/* ------------------------- */


__host__ void CUDAcreateArrayChromosome(thrust::host_vector<int> &array, struct parameters *params);

__global__ void CUDAcalculateChromosomeOutputs(int *solution, double *inputs, double *outputs, int numSamples, int numInputs, int numNodes);

__host__ double CUDAcalculateFitness(thrust::host_vector<int> &h_solution, thrust::device_vector<double> &d_inputs, thrust::device_vector<double> &d_outputs, int numSamples, int numInputs, int numNodes);

__host__ int *CUDAexecuteCGP(struct parameters *params, struct dataset *data, int popSize, int numGens);