#include "head.h"


// Symmetric positive-definite matrix
class matrix {

public:
	int getDimention() { return n; }
	
	int readFromFile(std::ifstream& fin);
	void writeToFile(std::ofstream& fout);

	int decomposeChol();
	vector <real> execDirectTraversal(vector<real> &F);
	vector <real> execReverseTraversal(vector<real> &y);

	void convAToDense();
	void convLToDense();

	vector <vector <real>> calcCompErrorForL();
	vector <real> calcCompErrorFory();
	vector <real> calcCompErrorForx();

private:
	vector <real> di, al;
	vector <int> ia;
	int n;

	vector <vector <real>> L, A;
};