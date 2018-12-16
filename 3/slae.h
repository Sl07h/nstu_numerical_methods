#include "head.h"



// System of linear algebraic equations
class SLAE {

public:
	int readSLAEfromFiles(const string &folderName);
	void writeDenseMatrixLToFile(std::ofstream& fout, const char *str);
	void writeDenseMatrixUToFile(std::ofstream& fout, const char *str);
	void convAToDense();
	void convLUToDense();
	void multAAndX();


	void getVectX(vector <real> &x) { x = F; };
	void generateVectX(int size);
	void writeXToFile(const char *fileName);
	void writeXToStream(std::ofstream& fout);
	void writeFToFile(const char *fileName);


	int getDimention() { return n; }
	void decomposionD();
	void decomposionChol();
	void decomposionLUsq();
	vector <real> execDirectTraversal(const vector <real> &_F);
	vector <real> execReverseTraversal(const vector <real> &_y);
	void createHilbertMatrix(int size);

	void LOS(std::ofstream& fout);
	void LOSfactD(std::ofstream& fout);
	void LOSfactLUsq(std::ofstream& fout);
	void clearAll();
	void setMaxiter(int new_maxiter) { maxiter = new_maxiter; }
	void setE(real new_E) { E = new_E; }

private:
	vector <real> multA(const vector <real>&x);
	vector <real> multD(const vector <real>&x);
	vector <real> multLUsq(const vector <real>&x);


	real calcRelativeDiscrepancy();
	real calcNormE(vector <real> &x);
	real calcScalarProduct(const vector <real> &a, const vector <real> &b);

	

	vector <vector <double>> L, U;
	vector <real> di, al, au, di_f, al_f, au_f;
	vector <int> ia, ja;
	int n, maxiter;
	real E;
	vector <real> x, r, z, p, F, Ftmp;

};