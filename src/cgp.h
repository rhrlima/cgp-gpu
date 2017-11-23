extern int teste;

typedef struct node {
	int func;
	int *links;
} NODE;

typedef struct genes {
	struct node *nodes;
	int *toEvaluate;
	double *nodesOutput;
	double fitness;
} GENES;

typedef struct config {
	int numInputs;
	int numOutputs;
	int numFunctions;
	char *functions;
	int numCols;
	int numRows;
	int levelsBack;
	int numGenes;
} CONFIG;

void hardCodedSolution(GENES *genes, int numGenes, int numOutputs);

double calculateFitness(GENES *genes, int numInputs, int numOutputs, int numGenes, double inputData[][2], double *outputData, int numCases);

void mutation(GENES parent, GENES *offspring, CONFIG params);

void freeSolution(GENES genes, int numGenes);

void printSolution(GENES genes, int numGenes, int numInputs, int numOutputs);

GENES * createPopulation(int popSize, CONFIG params);