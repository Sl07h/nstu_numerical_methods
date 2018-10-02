#define CATCH_CONFIG_MAIN 
#include "catch.hpp"
#include "head.h"
#include "SLAE.h"



TEST_CASE("Random-generated sparse matrix") {

	std::ifstream fin;
	std::ofstream fout;
	int max_width = 10, K_max = 7, size = 10;;
	
	matrix m;
	m.generateSparseMatrixA(size, 4);
	fout.open("A.txt");
	m.writeAToFile(fout);
	fout.close();

	fout.open("x.txt");
	for (int k = 0; k <= K_max; ++k) {
	
		
		fin.open("A.txt");
	
		SLAE slae;
		cout << "k = " << k << endl;
		fout << "k = " << k << endl;
		cout << "Started generating sparse matrix A" << endl;
		//slae.generateSparseMatrixA(k, i, max_width);
		//slae.readAFromFile(fin);
		slae.readAFromFile(fin);
		slae.addConditionNumber(k);
		cout << "Started generating vector x" << endl;
		slae.generateVectX(size);
		//slae.writexToFile(fout);

		cout << "Started mult A*x = F" << endl;
		slae.mult();
		
		slae.convAToDense();
		
		slae.writeMatrixtoFile(fout, "Matrix A");
		//slae.writeVectToFile(fout, "Vector F");

		cout << "Started LL' decomposion" << endl;
		slae.decomposeChol();
		//slae.convLToDense();
		//slae.writeMatrixtoFile(fout, "Matrix L");

		cout << "Started direct traversal" << endl;
		slae.execDirectTraversal();

		cout << "Started reverse traversal" << endl;
		slae.execReverseTraversal();
		//slae.writeVectToFile(fout, "Vector x");
		slae.writexCompError(fout, "x - x*");
		cout << endl;
		fin.close();
		CHECK(slae.isXcorrect());

	}
	

	fout.close();
}