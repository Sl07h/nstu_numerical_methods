#include "head.h"
#include "slae.h"


// Output dense matrix A
void SLAE::writeMatrixtoFile(std::ofstream& fout, char *str) {

	fout << std::fixed << std::setprecision(std::numeric_limits<real>::digits10 + 1) << str << endl;
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < A.size(); ++j) {
			fout << A[i][j] << "\t";
		}
		fout << ";" << endl;
	}
	fout << endl;
}


// Converting sparse matrix to dense format
void SLAE::convAToDense() {

	A.clear();
	A.resize(n);
	for (int i = 0; i < n; ++i) {
		A[i].resize(n, 0);
	}

	for (int i = 0; i < n; ++i) {

		A[i][i] = di[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0; k < i1; ++k, ++j) {

			A[i][j] = al[k];
			A[j][i] = al[k];
		}
	}
}


// Converting sparse matrix to dense format
void SLAE::convLToDense() {

	A.clear();
	A.resize(n);
	for (int i = 0; i < n; ++i) {
		A[i].resize(n, 0);
	}

	for (int i = 0; i < n; ++i) {

		A[i][i] = di[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0; k < i1; ++k, ++j) {

			A[i][j] = al[k];
		}
	}
}


// A*x = F
void SLAE::mult() {

	vector <real_sum> F_tmp(n, 0);

	for (int i = 0; i < n; ++i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0; k < i1; ++k, ++j) {

			F_tmp[i] += al[k] * F[j]; // В F лежит x
			F_tmp[j] += al[k] * F[i];
		}
	}


	for (int i = 0; i < n; ++i) {
		F_tmp[i] += di[i] * F[i];
	}

	for (int i = 0; i < n; ++i)
		F[i] = real(F_tmp[i]);
}


// Direct traversal (calculating y):	Ly = F		y = L\F
void SLAE::execDirectTraversal() {

	real_sum tmp;

	for (int i = 0; i < n; ++i) {

		tmp = 0.0;
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0; k < i1; ++k, ++j) {
			tmp += al[k] * F[j];
		}

		F[i] = (F[i] - tmp) / di[i];
	}
}


// Reverse traversal (calculating x):	L'x = y		x = L'\y
void SLAE::execReverseTraversal() {

	vector <real> x;
	x.resize(n, 0);

	for (int i = n - 1; i >= 0; --i) {

		x[i] = F[i] / di[i];

		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0;k < i1; ++k, ++j) {
			F[j] -= al[k] * x[i];
		}
	}

	F = x;
}



// Gaussian elemination
void SLAE::calcGauss() {

	// Приведение к верхне-треугольному виду
	for (int j = 0; j < A.size(); ++j) {
		for (int i = j + 1; i < A.size(); ++i) {

			real toMult = A[i][j] / A[j][j]; // Коэффициент, на который надо умножить строку

			for (int k = 0; k < A.size(); ++k)	// Отняли стоку
				A[i][k] -= toMult * A[j][k];

			F[i] -= toMult * F[j];
		}
	}


	// Обратный обход
	vector <real> x;
	x.resize(A.size(), 0);

	for (int i = n - 1; i >= 0; --i) {

		real_sum tmp = 0.0;
		for (int j = i + 1;j < A.size(); ++j) {
			tmp += A[i][j] * x[j];
		}
		x[i] = (F[i] - tmp) / A[i][i];
	}
	F = x;
}


// Gaussian elemination with leading element selection
void SLAE::calcGaussAdvanced() {

	// Приведение к верхне-треугольному виду
	for (int j = 0; j < A.size(); ++j) {
		for (int i = j + 1; i < A.size(); ++i) {

			int max = -DBL_MAX, row = i;
			for (int k = i; k < A.size(); ++k) // Ищем строку с ведущим элементом
				if (A[k][j] > max) {
					max = A[k][j];
					row = k;
				}
			
			std::swap(A[row], A[i]); // Меняем строку с  ведущим элементом на i-ю
			std::swap(F[row], F[i]);

			real toMult = A[i][j] / A[j][j]; // Коэффициент, на который надо умножить строку

			for (int k = 0; k < A.size(); ++k)	// Отняли стоку
				A[i][k] -= toMult * A[j][k];

			F[i] -= toMult * F[j];
		}
	}


	// Обратный обход
	vector <real> x;
	x.resize(A.size(), 0);

	for (int i = n - 1; i >= 0; --i) {

		real_sum tmp = 0.0;
		for (int j = i + 1;j < A.size(); ++j) {
			tmp += A[i][j] * x[j];
		}
		x[i] = (F[i] - tmp) / A[i][i];
	}
	F = x;
}