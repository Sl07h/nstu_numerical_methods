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
}


// LL' decomposion of the A matrix
int matrix::decomposeChol() {

	real_sum tmp;

	// ��� ��������� � ������� ������������, ��� �������������
	// ������ ������� ������������ �� �������� ����, ������� � �������
	for (int i = 0; i < n; ++i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);

		// ������� ��������� ������� ������������
		for (int k = i0; k < i1; ++k, ++j) {

			tmp = 0.0;
			int elem_i = ia[i]; // ����� �������� i-� ������
			int elem_j = ia[j]; // ����� �������� j-� ������
			int beg_i = i - (ia[i + 1] - ia[i]); // ������ ������� �������� i-� ������
			int beg_j = j - (ia[j + 1] - ia[j]); // ������ ������� �������� j-� ������

			int length_dif = beg_j - beg_i;

			if (length_dif >= 0)
				elem_i += length_dif;
			else
				elem_j += abs(length_dif);


			for (elem_i; elem_i < k; ++elem_i, ++elem_j)
				tmp += al[elem_i] * al[elem_j];


			al[k] = (al[k] - tmp) / di[j];
		}


		// ������� ������������� ��������
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
	int tmp_width; // ����� ������� ������

	for (int i = 0; i < n; ++i) {


		if (max_width >= i)
			tmp_width = i;
		else
			tmp_width = max_width;

		//if (tmp_width > 0)
			//prev += rand() % tmp_width; // ����� �� ���������� � ������� 2 ����

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
