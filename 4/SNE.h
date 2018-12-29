#include "head.h"


class SNE {

public:
	real calcNormE(vec &x);
	void calcVectorF(vec &x, vec &F);
	int inputSNE(int newTestNumber);
	void calcMatrixJAnalitically(vec &x, vec2 &J);
	void calcMatrixJNumerally(vec &x, vec2 &J);
	void deleteRow(int columnToDeleteCount, vec &F_tmp, vec &F, vec2 &Tmp, vec2 &J);
	void svertka(int columnToDeleteCount, vec &F_tmp, vec &F, vec2 &Tmp, vec2 &J);
	void calcGauss(vec &dx, vec &F, vec2 &J);
	void solveSNE();


private:
	int m;				// число функций (число строк)  
	int n;				// число переменных (число столбцов)
	int maxiter;		// максимальное число итераций
	int maxiterBeta;	// максимальное число итераций бетты
	real E1;			// Ограничение точности результата
	real E2;			// Ограничение бетты
	real normF_0;		// Норма начального приближения
	int method;			// Номер задания {2,3,7}
	int testNumber;		// Номер теста {1,...7}
	vec x, F;
};