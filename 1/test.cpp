#define CATCH_CONFIG_MAIN 
#include "catch.hpp"
#include "head.h"
#include "slae.h"



TEST_CASE("Condition number research") {

	std::ifstream fin;
	std::ofstream fout;
	fout.open("x_condNum.txt");
	int max_width = 10, K_max = 15;
	
	/*matrix m;
	m.generateSparseMatrixA(size, 4);
	fout.open("A.txt");
	m.writeAToFile(fout);
	fout.close();*/
	cout << "Condition number research:" << endl;

	for (int k = 0; k < K_max; ++k) {
	
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
		CHECK(slae.isXcorrect());
	}
	
	fout.close();
}



TEST_CASE("Hilbert matrix research") {

	std::ofstream fout;
	fout.open("x_Hilbert.txt");
	int K_max = 15;

	cout << "Hilbert matrix research:" << endl;
	for (int k = 0; k < K_max; ++k) {

		cout << "k = " << k << endl;
		SLAE slae;

		slae.createHilbertMatrix(k);
		slae.generateVectX(k);
		slae.mult();
		//slae.convAToDense();
		//slae.writeMatrixtoFile(fout, "H:");

		slae.decomposeChol();
		slae.execDirectTraversal();
		slae.execReverseTraversal();

		slae.writeTableToFile(fout);
		//slae.writeVectToFile(fout, "Vector x:");
		CHECK(slae.isXcorrect());
	}
	
	fout.close();
}




TEST_CASE("Gauss' method research") {

	std::ifstream fin;
	std::ofstream fout;
	fout.open("x_Gauss.txt");
	int K_max = 15;

	cout << "Gauss' method research:" << endl;
	for (int k = 0; k < K_max; ++k) {

		fin.open("A.txt");
		cout << "k = " << k << endl;
		SLAE slae;

		slae.readAFromFile(fin);
		slae.addConditionNumber(k);
		slae.generateVectX(slae.getDimention());

		slae.mult();

		slae.convAToDense();
		slae.calcGauss();

		slae.writeTableToFile(fout);
		//slae.writeVectToFile(fout, "Vector x:");
		fin.close();
		CHECK(slae.isXcorrect());
	}

	fout.close();
}


TEST_CASE("Advanced Gauss' method research") {

	std::ifstream fin;
	std::ofstream fout;
	fout.open("x_AdvGauss.txt");
	int K_max = 15;

	cout << "Advanced Gauss' method research:" << endl;
	for (int k = 0; k < K_max; ++k) {

		fin.open("A.txt");
		cout << "k = " << k << endl;
		SLAE slae;

		slae.readAFromFile(fin);
		slae.addConditionNumber(k);
		slae.generateVectX(slae.getDimention());

		slae.mult();

		slae.convAToDense();
		slae.calcGauss();

		slae.writeTableToFile(fout);
		//slae.writeVectToFile(fout, "Vector x:");
		fin.close();
		CHECK(slae.isXcorrect());
	}

	fout.close();
}