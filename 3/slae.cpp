#include "head.h"
#include "slae.h"


// Ввод СЛАУ: матрицы A и вектора y
int SLAE::readSLAEfromFiles(const string &folderName) {

	std::ifstream fin;
	fin.open(folderName + "/" + "kuslau.txt");
	fin >> n >> maxiter >> E;
	cout << E << endl;
	fin.close();
	x.resize(n, 0);
	r.resize(n, 0);
	z.resize(n, 0);
	p.resize(n, 0);
	F.resize(n, 0);
	Ftmp.resize(n, 0);

	fin.open(folderName + "/" + "di.txt");
	di.resize(n);
	for (int i = 0; i < di.size(); ++i) {
		fin >> di[i];
	}
	fin.close();

	fin.open(folderName + "/" + "ig.txt");
	ia.resize(n + 1);
	for (int i = 0; i < ia.size(); ++i) {
		fin >> ia[i];
	}
	fin.close();

	fin.open(folderName + "/" + "jg.txt");
	ja.resize(ia.back());
	for (int i = 0; i < ja.size(); ++i) {
		fin >> ja[i];
	}
	fin.close();


	fin.open(folderName + "/" + "ggl.txt");
	al.resize(ia.back());
	for (int i = 0; i < al.size(); ++i) {
		fin >> al[i];
	}
	fin.close();

	fin.open(folderName + "/" + "ggu.txt");
	au.resize(ia.back());
	for (int i = 0; i < au.size(); ++i) {
		fin >> au[i];
	}
	fin.close();



	fin.open(folderName + "/" + "pr.txt");
	F.resize(n);
	for (int i = 0; i < F.size(); ++i) {
		fin >> F[i];
	}
	fin.close();

	return 0;
}


// Вывод плотной матрицы L
void SLAE::writeDenseMatrixLToFile(std::ofstream& fout, const char *str) {

	//fout << std::fixed << std::setprecision(std::numeric_limits<real>::digits10 + 1);
	fout << str;
	for (int i = 0; i < L.size(); ++i) {
		for (int j = 0; j < L.size(); ++j) {
			fout << L[i][j] << "  ";
		}
		fout << ";" << endl;
	}
	fout << "]" << endl << endl;
}


// Вывод плотной матрицы U
void SLAE::writeDenseMatrixUToFile(std::ofstream& fout, const char *str) {

	//fout << std::fixed << std::setprecision(std::numeric_limits<real>::digits10 + 1);
	fout << str;
	for (int i = 0; i < U.size(); ++i) {
		for (int j = 0; j < U.size(); ++j) {
			fout << U[i][j] << "  ";
		}
		fout << ";" << endl;
	}
	fout << "]" << endl << endl;
}


// Преобразуем разряженную матрицу в плотный формат
void SLAE::convAToDense() {

	L.clear();
	L.resize(n);
	for (int i = 0; i < n; ++i) {
		L[i].resize(n, 0);
	}

	for (int i = 0; i < n; ++i) {

		L[i][i] = di[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];

		for (int k = i0; k < i1; ++k) { // Идём по всему профилю

			L[i][ja[k]] = al[k];
			L[ja[k]][i] = au[k];
		}
	}
}


// Преоразуем матрицы L и U в плотный формат
void SLAE::convLUToDense() {

	L.clear();
	L.resize(n);
	U.clear();
	U.resize(n);
	for (int i = 0; i < n; ++i) {
		L[i].resize(n, 0);
		U[i].resize(n, 0);
	}

	for (int i = 0; i < n; ++i) {

		L[i][i] = di_f[i];
		U[i][i] = di_f[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];

		for (int k = i0; k < i1; ++k) {

			L[i][ja[k]] = al_f[k];
			U[ja[k]][i] = au_f[k];
		}
	}
}


// A*x = F
void SLAE::multAAndX() {
	F = multA(x);
}


// A*x = b, где x - произвольный вектор
vec SLAE::multA(const vec&x) {

	vec result(n);

	for (int i = 0; i < n; ++i) {

		result[i] = di[i] * x[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];

		for (int k = i0; k < i1; ++k) {

			result[i] += al[k] * x[ja[k]];
			result[ja[k]] += au[k] * x[i];
		}
	}

	return result;
}

// A*x = b, где x - произвольный вектор
vec SLAE::multD(const vec&x) {

	vec result(n);

	for (int i = 0; i < n; ++i)
		result[i] = di_f[i] * x[i];

	return result;
}


// A*x = b, где x - произвольный вектор
vec SLAE::multLUsq(const vec&x) {

	vec result(n, 0);

	for (int i = 0; i < n; ++i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];

		for (int k = i0; k < i1; ++k) {

			result[i] += al_f[k] * x[ja[k]];
			result[ja[k]] += au_f[k] * x[i];
		}
	}

	for (int i = 0; i < n; ++i) {
		result[i] += di_f[i] * x[i];
	}

	return result;
}


// Создаём вектор x* = (1,2,...n)'
void SLAE::generateVectX(int size) {

	x.resize(size);
	for (int i = 0; i < size; ++i) {
		x[i] = i + 1;
	}
}


// Вывод вектора F в файл 
void SLAE::writeFToFile(const char *fileName) {

	std::ofstream fout;
	fout.open(fileName);

	for (int i = 0; i < F.size(); ++i)
		fout << F[i] << endl;
	fout.close();
}








void SLAE::writeXToFile(const char * fileName) {

	std::ofstream fout;
	fout.open(fileName);
	for (int i = 0; i < x.size(); ++i)
		fout << x[i] << " ";
	fout << " \t";
	fout.close();
}


void SLAE::writeXToStream(std::ofstream& fout) {

	for (int i = 0; i < x.size(); ++i)
		fout << x[i] << "\n";
	fout << "\n";
}





// Диагональное предобуславливание M = D
void SLAE::decomposionD() {

	di_f.clear();
	di_f.resize(n);
	for (int i = 0; i < n; ++i)
		di_f[i] = 1.0 / sqrt(di[i]);
}



// LU_sq разложение матрицы А
void SLAE::decomposionLUsq() {

	real sum_u, sum_l, sum_d;
	di_f = di;
	al_f = al;
	au_f = au;

	// Идём построчно в верхнем треугольнике, что экививалентно
	// Обходу нижнего треугольника по столбцам вниз, начиная с первого
	for (int i = 0; i < n; ++i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];

		// Рассчёт элементов нижнего треугольника
		for (int k = i0; k < i1; ++k) {

			int j = ja[k]; // текущий j
			int j0 = ia[j]; // i0 строки j
			int j1 = ia[j + 1]; // i1 строки j
			sum_l = 0;
			sum_u = 0;
			int ki = i0; // Индекс l_ik
			int kj = j0; // Индекс u_kj

			while (ki < k && kj < j1) {

				if (ja[ki] == ja[kj]) { // l_ik * u_kj
					sum_l += al_f[ki] * au_f[kj];
					sum_u += au_f[ki] * al_f[kj];
					ki++;
					kj++;
				}
				else { // Ищем следующие элементы i и j строки, которые можем перемножить
					if (ja[ki] > ja[kj]) kj++;
					else ki++;
				}
			}

			al_f[k] = (al_f[k] - sum_l) / di_f[j];
			au_f[k] = (au_f[k] - sum_u) / di_f[j];
		}


		// Рассчёт диагонального элемента
		sum_d = 0.0;
		for (int k = i0; k < i1; ++k)
			sum_d += al_f[k] * au_f[k];
		di_f[i] = sqrt(di_f[i] - sum_d);
	}
}



// LL' разложение матрицы А
void SLAE::decomposionChol() {

	real sum;
	di_f = di;
	al_f = al;
	au_f = au;

	// Идём построчно в верхнем треугольнике, что экививалентно
	// Обходу нижнего треугольника по столбцам вниз, начиная с первого
	for (int i = 0; i < n; ++i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];

		// Рассчёт элементов нижнего треугольника
		for (int k = i0; k < i1; ++k) {

			int j = ja[k]; // текущий j
			int j0 = ia[j]; // i0 строки j
			int j1 = ia[j + 1]; // i1 строки j
			sum = 0;
			int ki = i0; // Индекс l_ik
			int kj = j0; // Индекс u_kj

			while (ki < k && kj < j1) {

				if (ja[ki] == ja[kj]) { // l_ik * u_kj
					sum += al_f[ki] * al_f[kj];

					ki++;
					kj++;
				}
				else { // Ищем следующие элементы i и j строки, которые можем перемножить
					if (ja[ki] > ja[kj]) kj++;
					else ki++;
				}
			}

			al_f[k] = (al_f[k] - sum) / di_f[j];

		}


		// Рассчёт диагонального элемента
		sum = 0.0;
		for (int k = i0; k < i1; ++k)
			sum += al_f[k] * al_f[k];
		di_f[i] = sqrt(di_f[i] - sum);
	}
}


// Прямой ход    L y = F    ==>    y = L^-1 F
vec SLAE::execDirectTraversal(const vec &_F) {

	vec y;
	y.resize(n, 0);

	for (int i = 0; i < n; ++i) {
		real sum = 0;
		int i0 = ia[i];
		int i1 = ia[i + 1];
		for (int k = i0; k < i1; ++k)
			sum += al_f[k] * y[ja[k]];

		y[i] = (_F[i] - sum) / di_f[i];
	}
	return y;
}


// Обратный ход    U(sq) x = y    ==>    x = U(sq)^-1 y
vec SLAE::execReverseTraversal(const vec &_y) {

	vec x, y = _y;
	x.resize(n);
	for (int i = n - 1; i >= 0; --i) {

		x[i] = y[i] / di_f[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];
		for (int k = i0; k < i1; ++k)
			y[ja[k]] -= au_f[k] * x[i];
	}

	return x;
}


// Генерация матрицы Гильберта
void SLAE::createHilbertMatrix(int size) {

	clearAll();
	n = size;
	di.resize(n);
	ia.resize(n + 1);
	ja.resize(n*(n - 1) / 2);
	al.resize(n*(n - 1) / 2);
	au.resize(n*(n - 1) / 2);

	ia[0] = 0;
	ia[1] = 0;
	for (int i = 0; i < ia.size() - 1; ++i)
		ia[i + 1] = ia[i] + i;


	for (int i = 0; i < n; ++i) {
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = 0;
		for (int k = i0; k < i1; ++k, ++j) {
			ja[k] = j;
			al[k] = 1.0 / real(i + j + 1);
			au[k] = 1.0 / real(i + j + 1);
		}
		di[i] = 1.0 / real(i + j + 1);
	}
}


void SLAE::clearAll() {

	n = 0;

	di.clear();
	ia.clear();
	ja.clear();
	al.clear();
	au.clear();

	di_f.clear();
	al_f.clear();
	au_f.clear();

	x.clear();
	r.clear();
	z.clear();
	p.clear();
	F.clear();
	Ftmp.clear();
}


// Вычисление нормы в Евклидовом пространстве
real SLAE::calcNormE(vec &x) {
	return sqrt(x * x);
}



// Рассчёт относительной невязки
real SLAE::calcRelativeDiscrepancy() {
	return calcNormE(r) / calcNormE(F);
	//return (r*r);
}


// Локально - оптимальная схема
void SLAE::LOS(std::ofstream& fout) {

	int i;
	x.clear();			// Задаём начальное приближение
	x.resize(n, 0);		// x_0 = (0, 0, ...)
	r.resize(n);

	vec xprev = x;
	r = F - multA(x);	// r_0 = f - A*x_0
	z = r;				// z_0 = r_0
	p = multA(z);		// p_0 = A*z_0


	auto begin = std::chrono::steady_clock::now();
	for (i = 0; i < maxiter && calcRelativeDiscrepancy() >= E; ++i) {

		real pp = (p * p);
		real alpha = (p * r) / pp;

		x = x + alpha * z;
		r = r - alpha * p;

		Ftmp = multA(r);
		real beta = -(p * Ftmp) / pp;

		z = r + beta * z;
		p = Ftmp + beta * p;

		if (x == xprev)
			break;
		xprev = x;
		cout << alpha << "\t" << beta << endl;
	}


	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	fout << "LOS" << "\t" << i << "\t" << elapsed_ms.count() << "\t" << calcRelativeDiscrepancy() << endl;
}


// Локально - оптимальная схема c неполной диагональной факторизацией
void SLAE::LOSfactD(std::ofstream& fout) {

	int i;
	// Задаём начальное приближение
	x.clear();
	x.resize(n, 0);		// x_0 = (0, 0, ...)
	vec xprev = x;
	decomposionD();

	//vec xprev = x;
	r = F - multA(x);	// r_0 = f - A*x_0
	r = multD(r);

	z = multD(r);		// z = U^-1 r

	p = multA(z);		// p = A*z
	p = multD(p);		// p = L^-1 A*z


	auto begin = std::chrono::steady_clock::now();
	for (i = 0; i < maxiter && calcRelativeDiscrepancy() >= E; ++i) {

		real pp = p * p;
		real alpha = (p*r) / pp;
		x = x + alpha * z;
		r = r - alpha * p;

		vec tmp = multD(r);
		tmp = multA(tmp);
		tmp = multD(tmp);
		real beta = -(p * tmp) / pp;
		p = tmp + beta * p;

		tmp = multD(r);
		z = tmp + beta * z;

		/*if (x == xprev)
			break;*/
		xprev = x;
	}


	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	fout << "LOS + diag" << "\t" << i << "\t" << elapsed_ms.count() << "\t" << calcRelativeDiscrepancy() << endl;
}





// Локально - оптимальная схема с неполной факторизацией LU(sq)
void SLAE::LOSfactLUsq(std::ofstream& fout) {

	int i;
	x.clear();						// Задаём начальное приближение
	x.resize(n, 0);					// x_0 = (0, 0, ...)
	//vec xprev = x;
	decomposionLUsq();

	r = F - multA(x);				// r_0 = f - A*x_0
	r = execDirectTraversal(r);		// r = L^-1 (f - A*x_0)

	z = execReverseTraversal(r);	// z = U^-1 r

	p = multA(z);					// p = A*z
	p = execDirectTraversal(p);		// p = L^-1 A*z


	auto begin = std::chrono::steady_clock::now();
	for (i = 0; i < maxiter && calcRelativeDiscrepancy() >= E; ++i) {

		real pp = p * p;
		real alpha = (p*r) / pp;
		x = x + alpha * z;
		r = r - alpha * p;

		vec tmp = execReverseTraversal(r);
		tmp = multA(tmp);
		tmp = execDirectTraversal(tmp);
		real beta = -(p * tmp) / pp;
		p = tmp + beta * p;

		tmp = execReverseTraversal(r);
		z = tmp + beta * z;

		/*if (x == xprev)
			break;*/
			//xprev = x;
	}


	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	fout << "LOS + LU(sq)" << "\t" << i << "\t" << elapsed_ms.count() << "\t" << calcRelativeDiscrepancy() << endl;
}
