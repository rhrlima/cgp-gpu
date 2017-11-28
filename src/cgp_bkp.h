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
	int numNodes;		//Cols * Rows
	int numOutputs;		//Outputs for the dataset
	int arity;			//Max arity from functions (usually 2)
	int numFunctions;	//Available functions
	//int numCols;		//User defined number of columns
	//int numRows;		//User defined number of rows
	//int levelsBack;		//User defined back connections
};

typedef struct dataset {
	int numInputs;
	int numOutputs;
	int numCases;
	double **inputs;
	double **outputs;
} DATASET;


/* DATASET */
DATASET loadDataset(char *fileName);


/* CGP */
GENES createSolution(CONFIG p);

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

GENES createSolutionFromArray(int *vet, int numGenes, int arity, int numOutputs);