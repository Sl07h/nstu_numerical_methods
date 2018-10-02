#include "head.h"
#include "vect.h"


// Input of vector
int vect::readVectFromFile(std::ifstream& fin, int size) {

	F.resize(size);
	for (int i = 0; i < size; ++i) {
		fin >> F[i];
	}

	return 0;
}


// Creating vector X = (1,2,...n)'
void vect::generateVectX(int size) {

	F.resize(size);
	for (int i = 0; i < size; ++i) {
		F[i] = i + 1;
	}
}


// Output of results
void vect::writeVectToFile(std::ofstream& fout, char *str) {

	fout << std::fixed << std::setprecision(std::numeric_limits<real>::digits10 + 1) << str << endl;
	for (int i = 0; i < F.size();++i) {
		fout << F[i] << endl;
	}

	fout << endl;
}


// Output of computational error
void vect::writexCompError(std::ofstream& fout, char *str) {

	fout << std::fixed << std::setprecision(std::numeric_limits<real>::digits10 + 1) << str << endl;
	for (int i = 0; i < F.size();++i) {
		fout << F[i] - real(i + 1) << endl;
	}
	
	fout << endl;
}


// Check if X is almost equal to (1,2,...,n)'
bool vect::isXcorrect() {

	for (int i = 0; i < F.size();++i) {

		if (abs(F[i] - (real)(i + 1)) > std::numeric_limits<real>::digits10 + 2)
			return false;
	}

	return true;
}
