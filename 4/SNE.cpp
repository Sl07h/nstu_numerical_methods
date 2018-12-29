#include "SNE.h"




// Ввод данных
int SNE::inputSNE(int newTestNumber) {

	testNumber = newTestNumber;
	if (testNumber < 1 || testNumber > 7) {
		cout << "This test doesn't exist" << endl;
		return 1;
	}
	ifstream info("settings.txt"), x_0("x_0.txt");
	info >> m >> n >> maxiter >> maxiterBeta >> E1 >> E2 >> method;
	if (n > m) {
		cout << "Input error" << endl;
		return 2;
	}
	x.resize(n);
	for (int i = 0; i < n; i++)
		x_0 >> x[i];

	calcVectorF(x, F);
	normF_0 = calcNormE(F);
	return 0;
}


// Задание функций 
void SNE::calcVectorF(vec &x_tmp, vec &F_tmp) {

	F_tmp.resize(m);
	switch (testNumber)
	{
	case 1: {
		// Тест 1 (окружности не пересекаются) 
		F_tmp[0] = pow(x_tmp[0] + 3, 2) + pow(x_tmp[1] - 2, 2) - 4;
		F_tmp[1] = pow(x_tmp[0] - 3, 2) + pow(x_tmp[1] - 2, 2) - 4;
		break;
	}

	case 2: {
		// Тест 2 (одна точка пересечения) 
		F_tmp[0] = pow(x_tmp[0] + 2, 2) + pow(x_tmp[1] - 2, 2) - 4;
		F_tmp[1] = pow(x_tmp[0] - 2, 2) + pow(x_tmp[1] - 2, 2) - 4;;
		break;
	}

	case 3: {
		// Тест 3 (2 точки пересения) 
		F_tmp[0] = pow(x_tmp[0] + 1, 2) + pow(x_tmp[1] - 2, 2) - 4;
		F_tmp[1] = pow(x_tmp[0] - 1, 2) + pow(x_tmp[1] - 2, 2) - 4;
		break;
	}

	case 4: {
		// Тест 4 (2 окружности и прямая) 
		F_tmp[0] = pow(x_tmp[0] + 2, 2) + pow(x_tmp[1] - 2, 2) - 4;
		F_tmp[1] = pow(x_tmp[0] - 2, 2) + pow(x_tmp[1] - 2, 2) - 4;
		F_tmp[2] = x_tmp[0] + x_tmp[1] - 2;
		break;
	}

	case 5: {
		// Тест 5 (три попарно персекающиеся прямые) 
		F_tmp[0] = x_tmp[0];
		F_tmp[1] = x_tmp[1];
		F_tmp[2] = x_tmp[0] + x_tmp[1] - 4;
		break;
	}

	case 6: {
		// Тест 6 (три попарно персекающиеся прямые, одна взвешена)  
		F_tmp[0] = x_tmp[0];
		F_tmp[1] = x_tmp[1];
		F_tmp[2] = 100 * (x_tmp[0] + x_tmp[1] - 4);
		break;
	}

	case 7: {
		// Тест 7 (синусоида и прямая)  
		F_tmp[0] = sin(x_tmp[0]) - x_tmp[1];
		F_tmp[1] = x_tmp[0] - x_tmp[1] + 1;
		break;
	}
	}
}


// Задание матрицы Якоби аналитически
void SNE::calcMatrixJAnalitically(vec &x, vec2 &J) {

	switch (testNumber)
	{
	case 1: {
		// Тест 1 (окружности не пересекаются) 
		J[0][0] = 2 * (x[0] + 3); J[0][1] = 2 * (x[1] - 2);
		J[1][0] = 2 * (x[0] - 3); J[1][1] = 2 * (x[1] - 2);
		break;
	}

	case 2: {
		// Тест 2 (окружности пересекаются в одной точке)
		J[0][0] = 2 * (x[0] + 2); J[0][1] = 2 * (x[1] - 2);
		J[1][0] = 2 * (x[0] - 2); J[1][1] = 2 * (x[1] - 2);
		break;
	}

	case 3: {
		// Тест 3 (окружности пересекаются в двух точках) 
		J[0][0] = 2 * (x[0] + 1); J[0][1] = 2 * (x[1] - 2);
		J[1][0] = 2 * (x[0] - 1); J[1][1] = 2 * (x[1] - 2);
		break;
	}

	case 4: {
		// Тест 4 (2 окружности и прямая) 
		J[0][0] = 2 * (x[0] + 2); J[0][1] = 2 * (x[1] - 2);
		J[1][0] = 2 * (x[0] - 2); J[1][1] = 2 * (x[1] - 2);
		J[2][0] = 1; J[2][0] = 1;
		break;
	}

	case 5: {
		// Тест 5 (три попарно персекающиеся прямые)  
		J[0][0] = 1; J[0][1] = 0;
		J[1][0] = 0; J[1][1] = 1;
		J[2][0] = 1; J[2][1] = 1;
		break;
	}

	case 6: {
		// Тест 6 (три попарно персекающиеся прямые, одна взвешена)  
		J[0][0] = 1; J[0][1] = 0;
		J[1][0] = 0; J[1][1] = 1;
		J[2][0] = 100; J[2][1] = 100;
		break;
	}

	case 7: {
		// Тест 7 (синусоида и прямая)  
		J[0][0] = cos(x[0]);  J[0][1] = -1;
		J[1][0] = 1; J[1][1] = -1;
		break;
	}
	}
}


// Вычисление матрицы Якоби численно
void SNE::calcMatrixJNumerally(vec &x, vec2 &J) {

	real h = 1e-6;
	vec tmp, l1, l2;
	tmp.resize(n);
	l1.resize(m);
	l2.resize(m);
	for (int i = 0; i < n; i++)
		tmp[i] = x[i];

	calcVectorF(x, l1);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			tmp[j] += h;
			calcVectorF(tmp, l2);
			for (int k = 0; k < m; k++)
				l2[k] -= l1[k];
			J[i][j] = l2[i] / h;
			tmp[j] = x[j];
		}
	}
}


// Удаление ряда
void SNE::deleteRow(int columnToDeleteCount, vec &F_tmp, vec &F, vec2 &Tmp, vec2 &J) {
	int str_num; // номер убираемой строки  
	real min;

	for (int j = 0; j < columnToDeleteCount; j++) {
		min = fabs(F_tmp[0]);
		str_num = 0;
		for (int i = 1; i < m - j; i++) {
			if (fabs(F_tmp[i]) < min) {
				min = fabs(F_tmp[i]);
				str_num = i;
			}
		}
		swap(Tmp[str_num], Tmp[m - j - 1]);  //меняем строки матрицы 
		swap(F_tmp[str_num], F_tmp[m - j - 1]); //меняем функции 
	}

	for (int i = 0; i < n; i++) {
		F[i] = F_tmp[i];
		J[i] = Tmp[i];
	}
}


// Свёртка
void SNE::svertka(int columnToDeleteCount, vec &F_tmp, vec &F, vec2 &Tmp, vec2 &J) {

	int str_num; // номер убираемой строки
	real min, sum = 0;
	vec sum1;
	sum1.resize(n, 0);

	for (int j = 0; j < columnToDeleteCount; j++) {
		min = fabs(F_tmp[0]);
		str_num = 0;
		for (int i = 1; i < m - j; i++) {
			if (fabs(F_tmp[i]) < min) {
				min = fabs(F_tmp[i]);
				str_num = i;
			}
		}

		sum += pow(F_tmp[str_num], 2);
		for (int k = 0; k < n; k++)
			sum1[k] += F_tmp[str_num] * Tmp[str_num][k];
		swap(Tmp[str_num], Tmp[m - j - 1]);  //меняем строки матрицы  
		swap(F_tmp[str_num], F_tmp[m - j - 1]); //меняем функции 
	}
	F[m - columnToDeleteCount] = sum;
	for (int k = 0; k < n; k++)
		J[m - columnToDeleteCount][k] = 2 * sum1[k];

	for (int i = 0; i < n - 1; i++) {
		F[i] = F_tmp[i];
		J[i] = Tmp[i];
	}
}


// Метод Гаусса
void SNE::calcGauss(vec &dx, vec &F, vec2 &J) {

	int maxLine;
	real max, m;
	for (int i = 0; i < n; i++)
		dx[i] = -F[i];

	for (int k = 1; k < n; k++)
	{
		max = J[k - 1][k - 1];
		maxLine = k - 1;
		for (int l = k; l < n; l++)
			if (J[l][k - 1] > max) {
				max = J[l][k - 1];
				maxLine = l;
			}
		swap(J[maxLine], J[k - 1]);  //поменяли строки в матрице   
		swap(dx[maxLine], dx[k - 1]); //поменяли строки в приращении 
		for (int j = k; j < n; j++) {
			m = J[j][k - 1] / J[k - 1][k - 1];
			for (int i = 0; i < n; i++)
				J[j][i] = J[j][i] - m * J[k - 1][i];
			dx[j] = dx[j] - m * dx[k - 1];
		}
	}

	for (int i = n - 1; i >= 0; i--) {
		for (int j = i + 1; j < n; j++)
			dx[i] -= J[i][j] * dx[j];
		dx[i] = dx[i] / J[i][i];
	}
}


// Решение системы нелинейных уравнений
void SNE::solveSNE() {

	ofstream fout("result.txt");
	cout.precision(std::numeric_limits<real>::digits10 + 1);
	fout.precision(std::numeric_limits<real>::digits10 + 1);
	fout << "k\tx\ty\tbeta\tDiscrapancy\n";
	cout << "k\tx\ty\tbeta\tDiscrapancy\n";

	int columnToDeleteCount;
	vec F_tmp, x_k, x_k1, dx;
	vec2 J, J_tmp;
	J.resize(m);
	for (int i = 0; i < m; i++)
		J[i].resize(n);
	J_tmp.resize(m);
	for (int i = 0; i < m; i++)
		J_tmp[i].resize(n);
	F_tmp.resize(n);
	x_k = x;
	x_k1.resize(n);
	dx.resize(n);

	vec F_tmp1 = F;

	real normF_k = normF_0, normF = normF_0, normF1 = normF_0, beta;
	real nev = normF_k / normF_0;
	if (m - n != 0) {
		switch (method) {
		case 2: columnToDeleteCount = m - n; calcMatrixJAnalitically(x, J_tmp); break;
		case 3: columnToDeleteCount = m - n + 1; calcMatrixJAnalitically(x, J_tmp); break;
		case 7: columnToDeleteCount = m - n + 1; calcMatrixJNumerally(x, J_tmp); break;
		}
	}
	else //для первой итерации m=n 
	{
		columnToDeleteCount = 0;
		calcMatrixJAnalitically(x, J);
		/*calcMatrixJNumerally(x, J);*/
		for (int i = 0; i < m; i++)
			F_tmp[i] = F_tmp1[i];
	}

	
	for (int k = 0; k < maxiter && nev > E2; k++) {
		if (normF_k >= normF1 && k)
			break;
		if (m - n != 0 && k == 0) //для первой итерации 
		{
			switch (method) {
			case 2: deleteRow(columnToDeleteCount, F_tmp1, F_tmp, J_tmp, J); break;
			case 3: svertka(columnToDeleteCount, F_tmp1, F_tmp, J_tmp, J); break;
			case 7: svertka(columnToDeleteCount, F_tmp1, F_tmp, J_tmp, J); break;
			}
		}
		else {
			if (m - n != 0) {
				switch (method) {
				case 2:  calcMatrixJAnalitically(x_k, J_tmp); deleteRow(columnToDeleteCount, F_tmp1, F_tmp, J_tmp, J); break;
				case 3:  calcMatrixJAnalitically(x_k, J_tmp); svertka(columnToDeleteCount, F_tmp1, F_tmp, J_tmp, J); break;
				case 7:  calcMatrixJNumerally(x_k, J_tmp); svertka(columnToDeleteCount, F_tmp1, F_tmp, J_tmp, J); break;
				}
			}
			else if (k) {
				calcMatrixJAnalitically(x_k, J);
				/*calcMatrixJNumerally(x_k, J);*/
				for (int i = 0; i < m; i++)
					F_tmp[i] = F_tmp1[i];
			}
		}
		calcGauss(dx, F_tmp, J);
		normF = normF_k;
		normF1 = normF_k;
		beta = 1;
		for (int k_b = 0; k_b < maxiterBeta; k_b++) {
			for (int i = 0; i < n; i++)
				x_k1[i] = x_k[i] + beta * dx[i];
			calcVectorF(x_k1, F_tmp1);
			normF_k = calcNormE(F_tmp1);
			if (normF_k < normF) break;
			beta /= 2;
			if (E1 > beta) break;
		}
		for (int i = 0; i < n; i++)
			x_k[i] = x_k1[i];


		cout << k << "\t" << x_k << beta << "\t" << normF_k << endl;
		fout << k << "\t" << x_k << beta << "\t" << normF_k << endl;


		nev = normF_k / normF_0;
	}  fout.close();
}


// Вычисление Евклидовой нормы
real SNE::calcNormE(vec &x) {
	return sqrt(x*x);
}