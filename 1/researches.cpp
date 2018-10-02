#include "slae.h"


void main() {

	std::ifstream fin;
	std::ofstream fout;
	fout.open("x_condNum.txt");
	int K_max = 15;

	/*matrix m;
	m.generateSparseMatrixA(size, 4);
	fout.open("A.txt");
	m.writeAToFile(fout);
	fout.close();*/
	cout << "Condition number research:" << endl;

	for (int k = 0; k <= K_max; ++k) {

		fin.open("A.txt");
		cout << "k = " << k << endl;
		SLAE slae;

		slae.readAFromFile(fin);
		slae.addConditionNumber(k);
		slae.generateVectX(slae.getDimention());

		slae.mult();

		slae.decomposeChol();

		slae.execDirectTraversal();
		slae.execReverseTraversal();

		slae.writeTableToFile(fout);
		fin.close();
	}

	fout.close();



	cout << "Hilbert matrix research:" << endl;
	fout.open("x_Hilbert.txt");

	for (int k = 0; k <= K_max; ++k) {

		cout << "k = " << k << endl;
		SLAE slae;
		slae.createHilbertMatrix(k);
		slae.generateVectX(k);
		slae.mult();

		slae.decomposeChol();

		slae.execDirectTraversal();
		slae.execReverseTraversal();

		slae.writeTableToFile(fout);
	}
	fout.close();
}