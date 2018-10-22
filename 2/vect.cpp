#include "head.h"
#include "vect.h"



// Creating vector xk* = (1,2,...n)'
void vect::generateVectX(int size) {

	xk.resize(size);
	for (int i = 0; i < size; ++i) {
		xk[i] = i + 1;
	}
}




// Output of results
void vect::writeVectToFile(std::ofstream& fout, char *str) {

	fout << str << endl;
	fout << std::fixed << std::setprecision(std::numeric_limits<real>::digits10 + 1);
	for (int i = 0; i < xk.size();++i) {
		fout << xk[i] << endl;
	}
}


// generates 1/3 of table in research
void vect::writeTableToFile(std::ofstream& fout) {

	fout << std::fixed << std::setprecision(std::numeric_limits<real>::digits10 + 1);
	for (int i = 0; i < xk.size();++i) {
		fout << xk[i] << " ";
	}
	fout << " \t";

	fout << std::scientific;
	for (int i = 0; i < xk.size();++i) {
		fout << xk[i] - real(i + 1) << " ";
	}
	fout << " \t" << endl;
}


// Output of computational error
void vect::writexCompError(std::ofstream& fout, char *str) {

	fout << str << endl;
	fout << std::scientific;
	for (int i = 0; i < xk.size(); ++i) {
		fout << xk[i] - real(i + 1) << endl;
	}
}


// Check if xk is almost equal to (1,2,...,n)'
bool vect::isXcorrect() {

	for (int i = 0; i < xk.size();++i) {

		if (abs(xk[i] - (real)(i + 1)) > std::numeric_limits<real>::digits10 + 2)
			return false;
	}

	return true;
}
