#define CATCH_CONFIG_MAIN 
#include "catch.hpp"
#include "head.h"
#include "SLAE.h"



TEST_CASE("n = 10, 100, 1000", "MAX WIDTH = 5, 10, 15") {
		
	int max_width = 10;
	for (int i = 5000; i <= 5000000; i *= 10, max_width += 10) {

		SLAE slae;
  		slae.generateSparseMatrixA(i, max_width);
		
		
		slae.generateVectX(i);
		vector <real> x_correct;
		slae.getVectX(x_correct);


		slae.mult();


		slae.decomposeChol();
		slae.execDirectTraversal();
		slae.execReverseTraversal();


		CHECK(slae.isXcorrect());
	}
}