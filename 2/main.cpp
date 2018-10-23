#include "head.h"
#include "slae.h"


void main() {

	std::ofstream fout;
	fout.open("research.txt");
	SLAE slae;

	//slae.readMatrixFromFile("A.txt");
	slae.generateMatrixWith7Diagonals(10, 5);
	slae.convMatrixToDense();
	slae.generateVectX(slae.getDimention());
	slae.mult();

	slae.generateInitualGuess(slae.getDimention());


	slae.setE(0.0001);
	slae.setMaxiter(100000);

	cout << slae.calcIterative(1, 0.1) << endl;
	

	cout << slae.findOptimalW(1, true, fout) << endl;
	
	slae.writeDenseMatrixToFile("x.txt");
	fout.close();
}
