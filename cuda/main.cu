#include <stdio.h>

#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/functional.h>

#include "cgp.cuh"

int main() {

	unsigned int seed = time(NULL);
	srand(seed);

	char dataset_file[100];

	struct dataset *data;
	struct parameters *params;
	// struct chromosome *chromo, *best;
	struct chromosome *best;

	strcpy(dataset_file, "datasets/symbolic2_1024.data");	//x*x + x+x | (9, 2, 4)
	//strcpy(dataset_file, "datasets/symbolic3.data");	//(x0+x1) + (x0*x1) + (-x0)*(x1*x1) | (8, 2, 4)
	//strcpy(dataset_file, "datasets/symbolic4.data");  //x0*x1*x1 + x2*x1 + x*3

	data = loadDataset(dataset_file);

	params = initialiseParameters(9, 2, 4, data);//numNodes, maxArity, numFunctions

	printf("Dataset: '%s'\n", dataset_file);

	printParameters(params);

	// printf("Running CGP\n");
	// chromo = executeCGP(params, data, 10000);

	// printf("Best solution found\n");
	// printChromosome(chromo);

	int array[] = {1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 2, 4, 1, 4, 1, 2, 2, 5, 3, 2, 3, 6}; //best for dataset 2
	// // int array[] = {0, 0, 1, 2, 0, 1, 1, 0, 0, 1, 4, 0, 2, 1, 1, 2, 5, 6, 0, 2, 3, 0, 8, 7, 9}; //best for dataset 3
	// // int array[] = {2, , 1, 2};
	best = createChromosomeFromArray(params, array);
	calculateFitness(best, data);

	// printf("Best hardcoded\n");
	printChromosome(best);

	/* CUDA */
	int size = best->numNodes * (best->arity + 1) + best->numOutputs;
	int *solution = (int*)malloc(size * sizeof(int));
	createArrayFromChromosome(*best, solution);
	for(int i=0; i<size; i++)
		printf("%d ", solution[i]);
	printf("\n");
	

	int *d_solution;
	double *d_data_inputs;//, *d_data_outputs;

	printf("Cuda Malloc\n");
	cudaMalloc(&d_solution, 28 * sizeof(int));
	cudaMalloc(&d_data_inputs, data->numSamples * sizeof(double));
	// cudaMalloc(&d_data_outputs, data->numSamples * sizeof(double));

	printf("Memory Copy\n");
	cudaMemcpy(d_solution, solution, 28 * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_data_inputs, data->inputs, data->numSamples * sizeof(double), cudaMemcpyHostToDevice);
	// // cudaMemcpy(d_data_outputs, data->outputs, data->numSamples * sizeof(double), cudaMemcpyHostToDevice);

	// //setUpChromosomeData<<<4, 256>>>(dInPtr, dOutPrt, d_data_inputs, d_data_outputs, data->numSamples);

	thrust::device_vector<double> outputs(data->numSamples);
	double *dOutPrt = thrust::raw_pointer_cast(outputs.data());
	printf("A\n");
	teste<<<4, 256>>>(d_solution, d_data_inputs, dOutPrt, data->numSamples, 1, best->numNodes);
	printf("B\n");
	double error = 0.0;
	for(int i=0; i<data->numSamples; i++) {
		error += fabs(outputs[i] - data->outputs[i]);
	}
	printf("Fitness: %f\n", error);

	// freeChromosome(chromo);
	freeChromosome(best);

	freeDataset(data);
	free(params);

	return 0;
}