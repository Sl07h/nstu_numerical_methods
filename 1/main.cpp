#include "slae.h"


int main() {

	std::ifstream fin;
	std::ofstream fout;

	SLAE slae;

	// Input of matrix A and vector F
	fin.open("A.txt");
	if (slae.readAFromFile(fin) != 0) {
		cout << "Can't read matrix A";
		return -1;
	}
	fin.close();

	fin.open("F.txt");
	if (slae.readFFromFile(fin, slae.getDimention()) != 0) {
		cout << "Can't read vector F";
		return -1;
	}
	fin.close();


	// LL' decomposion
	if (slae.decomposeChol() != 0) {
		cout << "Can't use Cholesky decomposition";
		return -1;
	}


	// Direct and reverse traversal
	slae.execDirectTraversal();
	slae.execReverseTraversal();


	// Output of vector x
	fout.open("x.txt");
	slae.writexToFile(fout);
	fout.close();
}