#include "head.h"


// Symmetric positive-definite sparse matrix
class matrix {

public:
	
	int getDimention() { return n; }
	int readAFromFile(std::ifstream& fin);
	void writeAToFile(std::ofstream& fout);
	int decomposeChol();

	void generateSparseMatrixA(int n_new, int max_width);
	void createHilbertMatrix(int size);
	void addConditionNumber(int k) { di[0] += pow(10, -k); }

protected:
	vector <real> di, al;
	vector <int> ia;
	int n;
};