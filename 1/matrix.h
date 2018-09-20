#include "head.h"


// Symmetric positive-definite matrix
class matrix {

public:
	int getDimention() { return n; }
	int readAFromFile(std::ifstream& fin);
	int decomposeChol();

	void generateSparseMatrixA(int n, int max_width);

protected:
	vector <real> di, al;
	vector <int> ia;
	int n;
};