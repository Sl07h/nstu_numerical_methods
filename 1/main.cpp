#include "matrix.h"


int inputVectorF(
	std::ifstream& fin,
	vector<real> &F, int size
) {
	F.resize(size);
	for (int i = 0; i < F.size(); ++i) {
		fin >> F[i];
	}

	return 0;
}


void createHilbertMatrix(vector <vector <real>> & H) {

	for (int i = 0; H.size(); ++i) {
		for (int j = 0; H.size(); ++j) {
			H[i][j] = 1 / (i + j - 1);
		}
	}
}


int main() {

	std::ifstream fin;
	std::ofstream fout;
	matrix m;
	vector<real> F;


	// ¬вод матрицы A и вектора F
	fin.open("A.txt");
	if (m.readFromFile(fin) != 0) {
		cout << "Can't read matrix A";
		return -1;
	}
	fin.close();
	//m.convAToDense();

	fin.open("F.txt");
	if (inputVectorF(fin, F, m.getDimention()) != 0) {
		cout << "Can't read vector F";
		return -1;
	}
	fin.close();


	if (m.decomposeChol() != 0) {
		cout << "Can't use Cholesky decomposition";
		return -1;
	}
	m.convLToDense();


	vector<real> x, y;
	y = m.execDirectTraversal(F);
	x = m.execReverseTraversal(y);


	fout.open("out.txt");
	fout << std::fixed << std::setprecision(15);
	m.writeToFile(fout);
	fout << endl << "y:" << endl;
	for (int i = 0; i < m.getDimention(); ++i)
		fout << y[i] << endl;

	fout << endl << "X:" << endl;
	for (int i = 0; i < m.getDimention(); ++i)
		fout << x[i] << endl;



	fout.close();
}