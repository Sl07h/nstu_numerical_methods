#include "profMatrix.h"


// input sparse matrix A (n, di, ia, al)
int matrix::readAFromFile(std::ifstream& fin) {

	fin >> n;

	di.resize(n);
	for (int i = 0; i < di.size(); ++i) {
		fin >> di[i];
	}

	ia.resize(n + 1);
	for (int i = 0; i < ia.size(); ++i) {
		fin >> ia[i];
	}

	al.resize(ia.back());
	for (int i = 0; i < al.size(); ++i) {
		fin >> al[i];
	}

	return 0;
}


// Output sparse matrix A (n, di, ia, al)
void matrix::writeAToFile(std::ofstream& fout) {

	fout << n << endl;

	for (int i = 0; i < di.size(); ++i) {
		fout << di[i] << " ";
	}
	fout << endl;

	for (int i = 0; i < ia.size(); ++i) {
		fout << ia[i] << " ";
	}
	fout << endl;

	for (int i = 0; i < al.size(); ++i) {
		fout << al[i] << " ";
	}

	fout << endl;
}


// LL' decomposion of the A matrix
int matrix::decomposeChol() {

	real_sum tmp;

	// Идём построчно в верхнем треугольнике, что экививалентно
	// Обходу нижнего треугольника по столбцам вниз, начиная с первого
	for (int i = 0; i < n; ++i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);

		// Рассчёт элементов нижнего треугольника
		for (int k = i0; k < i1; ++k, ++j) {

			tmp = 0.0;
			int elem_i = ia[i]; // номер элемента i-й строки
			int elem_j = ia[j]; // номер элемента j-й строки
			int beg_i = i - (ia[i + 1] - ia[i]); // индекс первого элемента i-й строки
			int beg_j = j - (ia[j + 1] - ia[j]); // индекс первого элемента j-й строки

			int length_dif = beg_j - beg_i;

			if (length_dif >= 0)
				elem_i += length_dif;
			else
				elem_j += abs(length_dif);


			for (elem_i; elem_i < k; ++elem_i, ++elem_j)
				tmp += al[elem_i] * al[elem_j];


			al[k] = (al[k] - tmp) / di[j];
		}


		// Рассчёт диагонального элемента
		tmp = 0.0;
		for (int k = i0; k < i1; ++k)
			tmp += al[k] * al[k];
		di[i] = sqrt(di[i] - tmp);
	}

	return 0;
}


// Genereate sparse matrix A
void matrix::generateSparseMatrixA(int n_new, int max_width) {

	n = n_new;
	di.resize(n, 0);
	ia.resize(n + 1, 0);

	int prev = 0;
	int tmp_width; // Длина профиля строки

	for (int i = 0; i < n; ++i) {


		if (max_width >= i)
			tmp_width = i;
		else
			tmp_width = max_width;
		
		prev += tmp_width;

		ia[i + 1] = prev;

		int i0 = ia[i];
		int i1 = ia[i + 1];

		for (int k = i0; k < i1; ++k) {

			al.push_back(-(rand() % 5));
		}
	}


	vector <real_sum> tmp(n, 0);

	for (int i = 0; i < n; ++i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);

		for (int k = i0; k < i1; ++k, ++j) {
			tmp[i] += al[k];
			tmp[j] += al[k];
		}
	}

	for (int i = 0; i < n; ++i)
		di[i] = -tmp[i];
}


// Create D. Hilbert's matrix
void matrix::createHilbertMatrix(int size) {

	n = size;
	di.resize(n);
	ia.resize(n + 1);
	al.resize(n*(n - 1) / 2); // число элементов нижнего треугольника

	ia[0] = 0;
	ia[1] = 0;
	for (int i = 0; i < ia.size() - 1; ++i)
		ia[i + 1] = ia[i] + i;
	

	for (int i = 0; i < n; ++i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);
		for (int k = i0; k < i1; ++k, ++j) {

			al[k] = 1 / real(i + j + 1);
		}
		di[i] = 1 / real(i + j + 1);
	}
}
