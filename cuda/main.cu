#include <stdio.h>

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

	printf("Dataset: '%s'\n", dataset_file);

	printParameters(params);

	int array[] = {1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 2, 4, 1, 4, 1, 2, 2, 5, 3, 2, 3, 6};
	best = createChromosomeFromArray(params, array);
	calculateFitness(best, data);

	// printf("Best hardcoded\n");
	printChromosome(best);

	// /* CUDA */

	CUDAexecuteCGP(params, data, 5, 100);

	// freeChromosome(chromo);
	freeChromosome(best);

	freeDataset(data);
	free(params);

	return 0;
}