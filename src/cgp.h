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

void execute(int popSize, int numGen, CONFIG params, DATASET dataset);


/* UTILS */
int randint(int min, int max);

float randfloat(float min, float max);

NODE copyNode(NODE source, int numInputs);

GENES copySolution(GENES source, int numGenes, int numInputs, int numOutputs);

void freeSolution(GENES solution, int numGenes);

void freeDataset(DATASET dataset);

void printSolution(GENES genes, int numGenes, int numInputs, int numOutputs);

GENES hardCodedSolution();