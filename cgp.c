#include <stdio.h>

// typedef struct {
// 	int functionGene;
// 	int numInputs;
// 	int inputs[];
// } GENE;

// typedef struct {
// 	int numGenes;
// 	int 
// 	struct GENE *genes;

// }CHROMOSSOME;

// void nodesToProcess(G, NP) {
// 	for(int i=0; i<)
// }

int main() {

	int numGenes = 6;
	int numInputs = 2;
	int numOutputs = 4;

	int chromLenght = numGenes * (numInputs + 1) + numOutputs;// 22
	int chromossome[] = {0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 2, 5, 7, 3};

	int toEvaluate[numGenes+numInputs];

	printf("Cromossomo:\n");
	for(int i=0; i<chromLenght; i++) {
		if (i > 0 && i < chromLenght-numOutputs+1 && i % 3 == 0)
			printf("| ");
		printf("%d ", chromossome[i]);
	}
	printf("\nPara avaliar:\n");
	for(int i=0; i<numGenes+numInputs; i++) {
		toEvaluate[i] = 0;
		printf("%d ", toEvaluate[i]);
	}
	printf("\nAtivando outputs:\n");
	for(int i=chromLenght-numOutputs; i<chromLenght; i++) {
		toEvaluate[chromossome[i]] = 1;
	}

	for(int i=0; i<numGenes+numInputs; i++) {
		printf("%d ", toEvaluate[i]);
	}
	printf("\n");

	for(int i=numGenes+numInputs-1; i>=0; i--) {
		printf("Output: %d\n", i);
		if (toEvaluate[i]) {
			int index = (i) * (numInputs+1) - numInputs * (numInputs+1);
			if (i >= numInputs) {
				printf("Genes:\n");
				for(int j=0; j<numInputs+1; j++) {
					printf("%d:", chromossome[index+j]);
					if (j == 0) {
						printf(" function\n");
					} else {
						if (toEvaluate[chromossome[index+j]] == 0) {
							printf(" ativado\n");
							toEvaluate[chromossome[index+j]] = 1;
						} else {printf(" ok\n");}
					}
				}
			}
		} else {
			printf("Not activated\n");
		}
	}

	for(int i=0; i<numGenes+numInputs; i++) {
		printf("%d ", toEvaluate[i]);
	}
	printf("\n");

	return 0;
}