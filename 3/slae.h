#include "head.h"



// —истема линейных алгебраических уравнений
class SLAE {

public:
	int readSLAEfromFiles(const string &folderName, bool firstNumberIsOne);
	void writeSLAEToFiles(const string &folderName);
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
	void createHilbertMatricies(int a, int b, int step, const string &folderNameTemplate);

	pair<int, real> LOS();
	pair<int, real> LOSfactD();
	pair<int, real> LOSfactLUsq();
	void clearAll();
	void setMaxiter(int new_maxiter) { maxiter = new_maxiter; }
	void setE(real new_E) { E = new_E; }

private:
	vec multA(const vec&x);
	vec multD(const vec&x);

	real calcRelativeDiscrepancy();
	real calcNormE(vec &x);

	vector <vector <double>> L, U;
	vec di, al, au, di_f, al_f, au_f;
	vector <int> ia, ja;
	int n, maxiter;
	real E;
	vec x, r, z, p, F, Ftmp;
};