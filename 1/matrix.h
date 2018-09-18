#include "head.h"


// Симметричная положительно определённая матрица
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

	vector <vector <real>> calcDiscarapancyForL();
	vector <real> calcDiscarapancyFory();
	vector <real> calcDiscarapancyForx();

private:
	vector <real> di, al;
	vector <int> ia;
	int n;

	vector <vector <real>> L, A;
};