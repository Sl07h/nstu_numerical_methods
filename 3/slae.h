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


	void getVectX(vec &x) { x = F; };
	void generateVectX(int size);
	void writeXToFile(const char *fileName);
	void writeXToStream(std::ofstream& fout);
	void writeFToFile(const char *fileName);


	int getDimention() { return n; }
	void decomposionD();
	void decomposionChol();
	void decomposionLUsq();
	vec execDirectTraversal(const vec &_F);
	vec execReverseTraversal(const vec &_y);
	void createHilbertMatrix(int size);

	void LOS(std::ofstream& fout);
	void LOSfactD(std::ofstream& fout);
	void LOSfactLUsq(std::ofstream& fout);
	void clearAll();
	void setMaxiter(int new_maxiter) { maxiter = new_maxiter; }
	void setE(real new_E) { E = new_E; }

private:
	vec multA(const vec&x);
	vec multD(const vec&x);
	vec multLUsq(const vec&x);


	real calcRelativeDiscrepancy();
	real calcNormE(vec &x);

	

	vector <vector <double>> L, U;
	vec di, al, au, di_f, al_f, au_f;
	vector <int> ia, ja;
	int n, maxiter;
	real E;
	vec x, r, z, p, F, Ftmp;

};