#include "matrix.h"


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


// LL' decomposion of the A matrix
int matrix::decomposeChol() {

	real tmp;

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
void matrix::generateSparseMatrixA(int n, int max_width) {

	di.resize(n);
	ia.resize(n + 1);
	ia[0] = 0;
	ia[1] = 0;
	int prev = 0;

	for (int i = 0; i < n; ++i) {

		di[i] = (i % 3) * 100000;
		if (max_width >= i)
			max_width = i - 1;

		prev += rand() % max_width; // Чтобы не обращаться к массиву 2 раза
		ia[i + 1] = prev;

		int i0 = ia[i];
		int i1 = ia[i + 1];

		for (int k = i0; k < i1; ++k) {

			al[k] = rand() % 10 - 5;
		}
	}
}
