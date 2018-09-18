#include "matrix.h"


// Ввод входных данных
int matrix::readFromFile(std::ifstream& fin) {

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

	L.resize(n);
	for (int i = 0; i < n; ++i) {
		L[i].resize(n, 0);
	}

	return 0;
}


// Вывод результатов в файл
void matrix::writeToFile(std::ofstream& fout) {

	/*fout << n << endl;

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
	fout << endl;*/

	fout << "Matrix L:" << endl;
	for (int i = 0; i < n;++i) {
		for (int j = 0; j < n;++j) {
			fout << L[i][j] << " ";
		}
		fout << ";" << endl;
	}
}


// Преобразование разряженной матрицы A в плотную
void matrix::convAToDense() {

	for (int i = 0; i < n; ++i) {
		A[i][i] = di[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);
		for (int k = i0; k < i1; ++k, ++j)
			A[i][j] = al[k];
	}
}


// Преобразование разряженной матрицы L в плотную
void matrix::convLToDense() {

	for (int i = 0; i < n; ++i) {
		L[i][i] = di[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);
		for (int k = i0; k < i1; ++k, ++j)
			L[i][j] = al[k];
	}
}


// LL' разложнение матрицы A
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


// Вычисление невязки				A - L*L'
vector <vector <real>> matrix::calcDiscarapancyForL() {

	vector <vector <real>> dis;


	return dis;
}


// Вычисление невязки				F - L*y
vector <real> matrix::calcDiscarapancyFory() {

	vector <real> dis;


	return dis;
}


// Вычисление невязки				y - L'*x
vector <real> matrix::calcDiscarapancyForx() {

	vector <real> dis;


	return dis;
}


// Выполнеие прямого обхода:		Ly = F		y = L\F
vector <real> matrix::execDirectTraversal(vector<real> &F) {

	real tmp;
	vector <real> y;
	y.resize(n, 0);

	for (int i = 0; i < n; ++i) {

		tmp = 0.0;
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int k = i - (i1 - i0);

		for (int j = i0;j < i1;++j, ++k) {
			tmp += al[j] * y[k];
		}

		y[i] = (F[i] - tmp) / di[i];
	}

	return y;
}


// Выполнение обратного обхода:		L'x = y		x = L'\y
vector <real> matrix::execReverseTraversal(vector<real> &y) {
	// TO DO:
	/*
	Попробовать использовать дополнительную память tmp[y]
	*/


	vector <real> x, tmp;
	x.resize(n, 0);
	tmp.resize(n, 0);

	for (int i = n - 1; i >= 0; --i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];
		int k = i - (i1 - i0);

		for (int j = i0;j < i1; ++j, ++k) {
			//y[k] -= al[j] * x[i];
			tmp[k] += al[j] * x[i];
		}


	}

	for (int i = n - 1; i >= 0; --i) {
		x[i] = (y[i] - tmp[i]) / di[i];
	}
	return x;
}
