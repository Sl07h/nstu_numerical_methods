#include "head.h"


// Векторы xk and F
class vect {

public:

	void getVectX(vector <real> &x) { x = F; };

	void generateVectX(int size);
	// Создаём начальное приближение (нулевой вектор)
	void generateInitualGuess(int size) { xk.clear(); xk.resize(size, real(0)); }
	void writeTableToFile(std::ofstream& fout, int presision, real w, int iterations, int condNumber);
	void writeFToFile(char *fileName);
	bool isXcorrect();

protected:
	vector <real> xk, F;
};