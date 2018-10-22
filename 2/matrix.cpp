#include "matrix.h"


// ���� 7-������������ ������� �� �����
int matrix::readMatrixFromFile(char * fileName) {

	std::ifstream fin;
	fin.open(fileName);

	fin >> n >> m;
	di.resize(n);
	for (int i = 0; i < n; ++i)
		fin >> di[i];

	al1.resize(n - 1);
	for (int i = 0; i < al1.size(); ++i)
		fin >> al1[i];

	al2.resize(n - m - 2);
	for (int i = 0; i < al2.size(); ++i)
		fin >> al2[i];

	al3.resize(n - m - 3);
	for (int i = 0; i < al3.size(); ++i)
		fin >> al3[i];


	au1.resize(n - 1);
	for (int i = 0; i < au1.size(); ++i)
		fin >> au1[i];

	au2.resize(n - m - 2);
	for (int i = 0; i < au2.size(); ++i)
		fin >> au2[i];

	au3.resize(n - m - 3);
	for (int i = 0; i < au3.size(); ++i)
		fin >> au3[i];




	fin.close();
	return 0;
}


// ����� 7-�� ������������ ������� � ����
void matrix::writeMatrixToFile(char * fileName) {

	std::ofstream fout;
	fout.open(fileName);


	fout << n << m;
	for (int i = 0; i < n; ++i)
		fout << di[i] << " ";
	fout << endl;


	for (int i = 0; i < n - 1; ++i)
		fout << al1[i] << " ";
	fout << endl;

	for (int i = 0; i < n - m - 2; ++i)
		fout << al2[i] << " ";
	fout << endl;

	for (int i = 0; i < n - m - 3; ++i)
		fout << al3[i] << " ";
	fout << endl;


	for (int i = 0; i < n - 1; ++i)
		fout << au1[i] << " ";
	fout << endl;

	for (int i = 0; i < n - m - 2; ++i)
		fout << au2[i] << " ";
	fout << endl;

	for (int i = 0; i < n - m - 3; ++i)
		fout << au3[i] << " ";
	fout << endl;


	fout.close();
}


// ������ ������� A(k)
void matrix::generateMatrixWith7Diagonals(int new_n, int new_m) {

	n = new_n;
	m = new_m;
	di.resize(n, 0);
	au1.resize(n - 1);
	al1.resize(n - 1);

	au2.resize(n - m - 2);
	al2.resize(n - m - 2);

	au3.resize(n - m - 3);
	al3.resize(n - m - 3);

	for (int i = 0; i < al1.size(); ++i) {
		au1[i] = -rand() % 5;
		al1[i] = -rand() % 5;
	}

	for (int i = 0; i < al2.size(); ++i) {
		au2[i] = -rand() % 5;
		al2[i] = -rand() % 5;
	}

	for (int i = 0; i < al3.size(); ++i) {
		au3[i] = -rand() % 5;
		al3[i] = -rand() % 5;
	}


	for (int i = 0; i < di.size(); ++i)
		di[i] -= calcAii(i);

	di[0]++;
}


// ������� ����� ��������� ������ ��� ��������� ������� A(k)
real matrix::calcAii(int i) {

	real_sum sum = 0;

	if (i >= 1) { // ������ �����������
		sum += al1[i - 1];
		if (i >= m + 2) {
			sum += al2[i - m - 2];
			if (i >= m + 3)
				sum += al3[i - m - 3];
		}
	}

	sum += di[i];

	if (i < n - 1) { // ������� �����������
		sum += au1[i];
		if (i < n - m - 2) {
			sum += au2[i];
			if (i < n - m - 3)
				sum += au3[i];
		}
	}

	return sum;
}
