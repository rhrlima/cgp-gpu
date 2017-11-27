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
	double **inputs;
	double **outputs;
};


/* chromosome creation */
struct chromosome *createChromosome(struct parameters *params);

struct chromosome *createChromosomeFromArray(struct parameters *params, int *array);

struct node *createNode(int numInputs, int numNodes, int arity, int numFunctions, int nodePosition);

int getRandomNodeInput(int numChromoInputs, int numNodes, int nodePosition);

int getRandomChromosomeOutput(int numInputs, int numNodes);

void setChromosomeActiveNodes(struct chromosome *chromo);


/* mutation */
void singleMutation(struct parameters *params, struct chromosome *chromo);


/* chromosome evaluation */
void executeChromosome(struct chromosome *chromo, double *inputs);

double calculateFitness(struct chromosome *chromo, struct dataset *data);


/* evolutionary strategy */

/* 4 + 1 Evolutionary Strategy */
struct chromosome *executeCGP(struct parameters *params, struct dataset *data, int numGens);


/* dataset */
struct dataset *loadDataset(char *fileName);


/* utils */
int randint(int min, int max);

float randfloat(float min, float max);

struct parameters *initialiseParameters(int numNodes, int arity, int numFunctions, struct dataset *data);

struct node *copyNode(struct node *node);

struct chromosome *copyChromosome(struct chromosome *chromo);

void freeNode(struct node *node);

void freeChromosome(struct chromosome *chromo);

void printChromosome(struct chromosome *chromo);

void printParameters(struct parameters *params);