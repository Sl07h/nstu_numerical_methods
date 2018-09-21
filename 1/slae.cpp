#include "head.h"
#include "slae.h"


// Converting sparse matrix to dense format
void SLAE::convToDense() {

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


// A*x = F
void SLAE::mult_dense() {

	vector <real> F_tmp(n, 0);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j)
			F_tmp[i] = A[i][j] * F[j];// � F ����� x
	}

	F = F_tmp;
}


// A*x = F
void SLAE::mult() {

	vector <real> F_tmp(n, 0);

	for (int i = 0; i < n; ++i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0; k < i1; ++k, ++j) {

			F_tmp[i] += al[k] * F[j]; // � F ����� x
			F_tmp[j] += al[k] * F[i];
		}
	}


	for (int i = 0; i < n; ++i) {
		F_tmp[i] += di[i] * F[i];
	}

	F = F_tmp;
}


// Create D. Hilbert's matrix
void SLAE::createHilbertMatrix() {

	for (int i = 0; A.size(); ++i) {
		for (int j = 0; A.size(); ++j) {
			A[i][j] = 1 / (i + j - 1);
		}
	}
}


// Direct traversal (calculating y):	Ly = F		y = L\F
void SLAE::execDirectTraversal() {

	real tmp;

	for (int i = 0; i < n; ++i) {

		tmp = 0.0;
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int k = i - (i1 - i0);

		for (int j = i0;j < i1;++j, ++k) {
			tmp += al[j] * F[k];
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
