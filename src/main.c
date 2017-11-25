#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "cgp.h"

int main() {

	srand(time(NULL)); //No seed
	//srand(2);

	CONFIG params;
	params.numInputs = 2;
	params.numOutputs = 4;
	params.numFunctions = 4;
	params.functions;
	params.numCols = 3;
	params.numRows = 2;
	params.levelsBack = 5;
	params.numGenes = 6;

	//Y = X0+X1 + X0*X1 + -X0*(X1*X1) + 0
	DATASET dset = loadDataset("oi");
	
	execute(5, 100, params, dset);

	GENES G1 = hardCodedSolution();
	//printSolution(G1, 6, 2, 4);

	freeSolution(G1, 6);
	freeDataset(dset);

	return 0;
}