#include "head.h"
#include "slae.h"


void main() {

	SLAE slae;

	//slae.readMatrixFromFile("A.txt");
	slae.generateMatrixWith7Diagonals(10, 5);
	slae.convMatrixToDense();
	slae.generateVectX(slae.getDimention());
	slae.mult();
	
	slae.generateInitualGuess(slae.getDimention());
	//slae.generateMatrixWithDiagonal(8, 3);
	for (int i = 0; i < 10000; ++i)
	slae.calcIterative(1, 0.5);
	//slae.GaussSeildel(0.2);
	slae.setE(0.00001);
	slae.setMaxiter(1000);
	
	//cout << slae.findOptimalW(1) << endl;
	//slae.calcIterative(1, 0.5);
	slae.writeDenseMatrixToFile("x.txt");
	system("pause");
}
