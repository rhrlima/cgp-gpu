extern int teste;

typedef struct node {
	int func;
	int *links;
} NODE;

typedef struct genes {
	struct node *nodes;
	int *toEvaluate;
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

typedef struct dataset {
	int numInputs;
	int numOutputs;
	int numCases;
	double **inputs;
	double *outputs;
} DATASET;


/* DATASET */
DATASET loadDataset(char *fileName);


/* CGP */
GENES * createPopulation(int popSize, CONFIG params);

void calculateFitness(GENES *solution, int numInputs, int numOutputs, int numGenes, DATASET dataset);

GENES mutation(GENES parent, CONFIG params);


/* UTILS */
int randint(int min, int max);

float randfloat(float min, float max);

void printSolution(GENES genes, int numGenes, int numInputs, int numOutputs);

void freeSolution(GENES genes, int numGenes);

void hardCodedSolution(GENES *genes, int numGenes, int numOutputs);