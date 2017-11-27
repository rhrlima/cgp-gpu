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

struct dataset {
	int numInputs;
	int numOutputs;
	int numSamples;
	double **inputs;
	double **outputs;
};


/* chromosome creation */
struct node *createNode(int numInputs, int numNodes, int arity, int numFunctions, int nodePosition);

int getRandomNodeInput(int numChromoInputs, int numNodes, int nodePosition);

int getRandomNodeOutput(int numInputs, int numNodes);

void setChromosomeActiveNodes(struct chromosome *chromo);

struct chromosome *createChromosome(struct parameters *params);


/* chromosome evaluation */
void executeChromosome(struct chromosome *chromo, double *inputs);

double calculateFitness(struct chromosome *chromo, struct dataset *data)

/* dataset */
struct dataset *loadDataset(char *fileName);


/* utils */
void printChromosome(struct chromosome *chromo);