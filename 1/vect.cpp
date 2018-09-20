#include "head.h"
#include "vect.h"


// Input of vector F
int vect::readFFromFile(std::ifstream& fin, int size) {

	F.resize(size);
	for (int i = 0; i < size; ++i) {
		fin >> F[i];
	}

	return 0;
}


// Creating vector X = (1,2,...n)'
void vect::generateVectX(int size) {

	F.resize(size);
	for (int i = 1; i <= size; ++i) {
		F[i] = i;
	}
}


// Output of results
void vect::writexToFile(std::ofstream& fout) {

	fout << std::fixed << std::setprecision(std::numeric_limits<real>::digits10 + 1) << "x:" << endl;
	for (int i = 0; i < F.size();++i) {
		fout << F[i] << endl;
	}
}