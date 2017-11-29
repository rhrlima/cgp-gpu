#include <stdio.h>

#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

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

	// printf("Dataset: '%s'\n", dataset_file);

	// printParameters(params);

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

	/* test */
	best->fitness = -1;

	thrust::device_vector<double> outputs(data->numSamples);
	double *out = thrust::raw_pointer_cast(outputs.data());

	//4 x 256 = 1024 samples
	// thrust::fill(fitnesses.begin(), fitnesses.end(), 0);

	int numThreads = 256;
	int numBlocks = ceil((float)data->numSamples/numThreads);
	cudaCalculateFitnesses<<<4, 256>>>(*best, out, data->inputs, data->numSamples);
	
	for(int i=0; i<data->numSamples; i++) {
		printf("%4d %6.2f\n", i, (double) outputs[i]);
	}

	// double error = thrust::reduce(fitnesses.begin(), fitnesses.end());
	// printf("Fitness: %f\n", error);

	// freeChromosome(chromo);
	freeChromosome(best);

	freeDataset(data);
	free(params);

	//---------------------

	return 0;
}