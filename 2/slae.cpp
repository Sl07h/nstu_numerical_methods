#include "slae.h"

// Преобразование 7-ми диагональной матрицы в плотный формат
void SLAE::convMatrixToDense() {

	A.clear();
	A.resize(n);

	for (int i = 0; i < n; ++i) {
		A[i].resize(n, 0);
		A[i][i] = di[i];
	}


	int j = 1;
	for (int i = 0; i < al1.size(); ++i, ++j) {

		A[i][j] = au1[i];
		A[j][i] = al1[i];
	}

	j = m + 2;
	for (int i = 0; i < al2.size(); ++i, ++j) {

		A[i][j] = au2[i];
		A[j][i] = al2[i];
	}

	j = m + 3;
	for (int i = 0; i < al3.size(); ++i, ++j) {

		A[i][j] = au3[i];
		A[j][i] = al3[i];
	}
}


// Вывод плотной матрицы в файл
void SLAE::writeDenseMatrixToFile(char * fileName) {

	std::ofstream fout;
	fout.open(fileName);


	for (int i = 0; i < n;++i) {
		for (int j = 0; j < n; ++j)
			fout << A[i][j] << "\t";
		fout << endl;
	}

	fout.close();
}


// Умножение i-й строки матрицы на вектор
real SLAE::multLine(vector <real> &line, int i, int mode) {

	real_sum sum = 0;
	if (mode == 1 || mode == 3) {	// Нижний треугольник

		if (i > 0) {

			sum += al1[i - 1] * line[i - 1];
			if (i > m + 1) {
				sum += al2[i - m - 2] * line[i - m - 2];
				if (i > m + 2)
					sum += al3[i - m - 3] * line[i - m - 3];
			}
		}
	}


	if (mode == 2 || mode == 3) {	// Главная диагональ
									// и верхний треугольник
		sum += di[i] * line[i];
		if (i < n - 1) {

			sum += au1[i] * line[i + 1];

			if (i < n - m - 2) {
				sum += au2[i] * line[i + m + 2];
				if (i < n - m - 3)
					sum += au3[i] * line[i + m + 3];
			}
		}
	}

	return sum;
}


// Умножение матрицы на вектор
void SLAE::mult() {

	F.resize(n, 0);
	// Нижний треугольник
	for (int i = 0; i < al1.size(); ++i)
		F[i + 1] += al1[i] * xk[i];
	for (int i = 0; i < al2.size(); ++i)
		F[i + m + 2] += al2[i] * xk[i];
	for (int i = 0; i < al3.size(); ++i)
		F[i + m + 3] += al3[i] * xk[i];

	// Главная диагональ
	for (int i = 0; i < di.size(); ++i)
		F[i] += di[i] * xk[i];

	// Верхний треугольник
	for (int i = 0; i < au1.size(); ++i)
		F[i] += au1[i] * xk[i + 1];
	for (int i = 0; i < au2.size(); ++i)
		F[i] += au2[i] * xk[i + m + 2];
	for (int i = 0; i < au3.size(); ++i)
		F[i] += au3[i] * xk[i + m + 3];
}


// Метод Якоби. 0 < w < 1
// Используется общая память для xk и xk1
void SLAE::Jacobi(real w) {

	real_sum sum;

	for (int i = 0; i < n; ++i) {
		sum = multLine(xk, i, 3);
		xk[i] += w * (F[i] - sum) / di[i];
		//xk1[i] = xk[i] + w * (F[i] - sum) / di[i];
	}
}


// Метод Гаусса-Зейделя. 0 < w < 2
void SLAE::GaussSeildel(real w) {

	real_sum sum;
	vector <real> xk1 = xk;

	for (int i = 0; i < n; ++i) {

		sum = multLine(xk1, i, 1);
		sum += multLine(xk, i, 2);

		xk1[i] = xk[i] + w * (F[i] - sum) / di[i];
	}
	xk = xk1;
}


// Решение СЛАУ итерационным методом
// 1  Метод Якоби
// 2  Метод Гаусса-Зейделя
int SLAE::calcIterative(int mode, real w) {
	int i = 0;
	while (i < maxiter && calcRelativeDiscrepancy() >= E) {

		if (mode == 1)
			Jacobi(w);
		else
			GaussSeildel(w);

		++i;
	}
	return i;
}

real SLAE::findOptimalW(int mode) {

	real maxIter;
	if (mode == 1) maxIter = 100;
	else maxIter = 201;

	for (int i = 0; i < maxIter; ++i)
		calcIterative(mode, real(i) / 100);

	return real();
}



vector <real_sum>& SLAE::sumAtlowerTrianByDiag(vector <real>& line) {

	vector <real_sum> sum;
	sum.resize(n, 0);

	// Нижний треугольник
	for (int i = 0; i < al1.size(); ++i)
		sum[i + 1] += al1[i] * line[i];
	for (int i = 0; i < al2.size(); ++i)
		sum[i + m + 2] += al2[i] * line[i];
	for (int i = 0; i < al3.size(); ++i)
		sum[i + m + 3] += al3[i] * line[i];

	return sum;
}



vector <real_sum>& SLAE::sumAtUpperTrianByDiag(vector <real>& line) {

	vector <real_sum> sum;
	sum.resize(n, 0);

	// Главная диагональ
	for (int i = 0; i < di.size(); ++i)
		F[i] += di[i] * line[i];

	// Верхний треугольник
	for (int i = 0; i < au1.size(); ++i)
		F[i] += au1[i] * line[i + 1];
	for (int i = 0; i < au2.size(); ++i)
		F[i] += au2[i] * line[i + m + 2];
	for (int i = 0; i < au3.size(); ++i)
		F[i] += au3[i] * line[i + m + 3];

	return sum;
}


// Метод Якоби с кешированием. 0 < w < 1
void SLAE::BoostedJacobi(real w) {

	vector <real> xk1, sum1, sum2;
	xk1.resize(n, 0);
	sum1 = sumAtlowerTrianByDiag(xk);
	sum2 = sumAtUpperTrianByDiag(xk);

	generateInitualGuess(n); // Создаём начальное приближение

	for (int i = 0; i < n; ++i) {
		xk1[i] = xk[i] + w * (F[i] - sum1[i] - sum2[i]) / di[i];

	}

	xk = xk1;
}



// Метод Гаусса-Зейделя с кешированием. 0 < w < 2
void SLAE::BoostedGaussSeildel(real w) {

	real_sum sum = 0;
	vector <real> xk1, vecOfSum;
	xk1.resize(n, 0);
	vecOfSum = sumAtUpperTrianByDiag(xk);


	for (int i = 0; i < n; ++i) {

		sum = multLine(xk1, i, 1);

		xk1[i] = xk[i] + w * (F[i] - vecOfSum[i] - sum) / di[i];
	}
	xk = xk1;
}



// Вычисление нормы в Евклидовом пространстве
real SLAE::calcNormE(vector <real> &x) {

	real_sum norma = 0;
	for (int i = 0; i < n; i++)
		norma += x[i] * x[i];

	return sqrt(norma);
}


// Рассчёт относительной невязки
real SLAE::calcRelativeDiscrepancy() {

	vector <real> numerator, denominator = F;
	numerator.resize(n);
	mult(); // F = A*xk
	for (int i = 0; i < n; ++i)
		numerator[i] = denominator[i] - F[i];

	return calcNormE(numerator) / calcNormE(denominator);
}


// Рассчёт числа обусловленности
int SLAE::calcCondNumber() {

	vector <real> numerator, denominator;
	numerator.resize(n);
	denominator.resize(n);
	for (int i = 0; i < n; ++i) {
		numerator[i] = xk[i] - (i + 1);
		denominator[i] = i + 1;
	}


	return calcNormE(numerator) / (calcNormE(denominator) * calcRelativeDiscrepancy());
}

