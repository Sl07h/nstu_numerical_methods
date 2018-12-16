#include "head.h"
#include "slae.h"


// ���� ����: ������� A � ������� y
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


// ����� ������� ������� L
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


// ����� ������� ������� U
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


// ����������� ����������� ������� � ������� ������
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

		for (int k = i0; k < i1; ++k) { // ��� �� ����� �������

			L[i][ja[k]] = al[k];
			L[ja[k]][i] = au[k];
		}
	}
}


// ���������� ������� L � U � ������� ������
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


// A*x = b, ��� x - ������������ ������
vector <real> SLAE::multA(const vector <real>&x) {

	vector <real> result(n);

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

// A*x = b, ��� x - ������������ ������
vector <real> SLAE::multD(const vector <real>&x) {

	vector <real> result(n);

	for (int i = 0; i < n; ++i)
		result[i] = di_f[i] * x[i];

	return result;
}


// A*x = b, ��� x - ������������ ������
vector <real> SLAE::multLUsq(const vector <real>&x) {

	vector <real> result(n, 0);

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


// ������ ������ x* = (1,2,...n)'
void SLAE::generateVectX(int size) {

	x.resize(size);
	for (int i = 0; i < size; ++i) {
		x[i] = i + 1;
	}
}


// ����� ������� F � ���� 
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





// ������������ ������������������ M = D
void SLAE::decomposionD() {

	di_f.clear();
	di_f.resize(n);
	for (int i = 0; i < n; ++i)
		di_f[i] = 1.0 / sqrt(di[i]);
}



// LU_sq ���������� ������� �
void SLAE::decomposionLUsq() {

	real sum_u, sum_l, sum_d;
	di_f = di;
	al_f = al;
	au_f = au;

	// ��� ��������� � ������� ������������, ��� �������������
	// ������ ������� ������������ �� �������� ����, ������� � �������
	for (int i = 0; i < n; ++i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];

		// ������� ��������� ������� ������������
		for (int k = i0; k < i1; ++k) {

			int j = ja[k]; // ������� j
			int j0 = ia[j]; // i0 ������ j
			int j1 = ia[j + 1]; // i1 ������ j
			sum_l = 0;
			sum_u = 0;
			int ki = i0; // ������ l_ik
			int kj = j0; // ������ u_kj

			while (ki < k && kj < j1) {

				if (ja[ki] == ja[kj]) { // l_ik * u_kj
					sum_l += al_f[ki] * au_f[kj];
					sum_u += au_f[ki] * al_f[kj];
					ki++;
					kj++;
				}
				else { // ���� ��������� �������� i � j ������, ������� ����� �����������
					if (ja[ki] > ja[kj]) kj++;
					else ki++;
				}
			}

			al_f[k] = (al_f[k] - sum_l) / di_f[j];
			au_f[k] = (au_f[k] - sum_u) / di_f[j];
		}


		// ������� ������������� ��������
		sum_d = 0.0;
		for (int k = i0; k < i1; ++k)
			sum_d += al_f[k] * au_f[k];
		di_f[i] = sqrt(di_f[i] - sum_d);
	}
}



// LL' ���������� ������� �
void SLAE::decomposionChol() {

	real sum;
	di_f = di;
	al_f = al;
	au_f = au;

	// ��� ��������� � ������� ������������, ��� �������������
	// ������ ������� ������������ �� �������� ����, ������� � �������
	for (int i = 0; i < n; ++i) {

		int i0 = ia[i];
		int i1 = ia[i + 1];

		// ������� ��������� ������� ������������
		for (int k = i0; k < i1; ++k) {

			int j = ja[k]; // ������� j
			int j0 = ia[j]; // i0 ������ j
			int j1 = ia[j + 1]; // i1 ������ j
			sum = 0;
			int ki = i0; // ������ l_ik
			int kj = j0; // ������ u_kj

			while (ki < k && kj < j1) {

				if (ja[ki] == ja[kj]) { // l_ik * u_kj
					sum += al_f[ki] * al_f[kj];

					ki++;
					kj++;
				}
				else { // ���� ��������� �������� i � j ������, ������� ����� �����������
					if (ja[ki] > ja[kj]) kj++;
					else ki++;
				}
			}

			al_f[k] = (al_f[k] - sum) / di_f[j];

		}


		// ������� ������������� ��������
		sum = 0.0;
		for (int k = i0; k < i1; ++k)
			sum += al_f[k] * al_f[k];
		di_f[i] = sqrt(di_f[i] - sum);
	}
}


// ������ ���    L y = F    ==>    y = L^-1 F
vector <real> SLAE::execDirectTraversal(const vector <real> &_F) {

	vector <real> y;
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


// �������� ���    U(sq) x = y    ==>    x = U(sq)^-1 y
vector <real> SLAE::execReverseTraversal(const vector <real> &_y) {

	vector <real> x, y = _y;
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


// ��������� ������� ���������
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


// ���������� ����� � ���������� ������������
real SLAE::calcNormE(vector <real> &x) {
	return sqrt(calcScalarProduct(x, x));
}



// ������� ������������� �������
real SLAE::calcRelativeDiscrepancy() {
	return calcNormE(r) / calcNormE(F);
	//return calcScalarProduct(r,r);
}


// ���������� ���������� ������������
real SLAE::calcScalarProduct(const vector <real> &a, const vector <real> &b) {

	real result = 0;
	for (int i = 0; i < a.size(); ++i)
		result += a[i] * b[i];

	return fabs(result);
}


// �������� - ����������� �����
void SLAE::LOS(std::ofstream& fout) {

	int i;
	x.clear();			// ����� ��������� �����������
	x.resize(n, 0);		// x_0 = (0, 0, ...)

	r = multA(x);		// r_0 = f - A*x_0
	vector <real> xprev = x;
	for (int i = 0; i < n; ++i)
		r[i] = F[i] - r[i];
	z = r;				// z_0 = r_0
	p = multA(z);		// p_0 = A*z_0


	auto begin = std::chrono::steady_clock::now();
	for (i = 0; i < maxiter && calcRelativeDiscrepancy() >= E; ++i) {

		real pp = calcScalarProduct(p, p);
		real alpha = calcScalarProduct(p, r) / pp;
		for (int j = 0; j < n; ++j) {
			x[j] += alpha * z[j];
			r[j] -= alpha * p[j];
		}

		Ftmp = multA(r);
		real beta = -calcScalarProduct(p, Ftmp) / pp;

		for (int j = 0; j < n; ++j) {
			z[j] = r[j] + beta * z[j];
			p[j] = Ftmp[j] + beta * p[j];
		}

		if (x == xprev)
			break;
		xprev = x;
	}


	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	fout << "LOS" << "\t" << i << "\t" << elapsed_ms.count() << "\t" << calcRelativeDiscrepancy() << endl;
}


// �������� - ����������� ����� c �������� ������������ �������������
void SLAE::LOSfactD(std::ofstream& fout) {

	int i;
	// ����� ��������� �����������
	x.clear();
	x.resize(n, 0);		// x_0 = (0, 0, ...)
	vector <real> xprev = x;
	decomposionD();

	r = multA(x);		// r = A*x_0
	for (int i = 0; i < n; ++i)	// r = f - A*x_0
		r[i] = F[i] - r[i];
	r = multD(r);

	z = multD(r);		// z = U^-1 r

	p = multA(z);		// p = A*z
	p = multD(p);		// p = L^-1 A*z


	auto begin = std::chrono::steady_clock::now();
	for (i = 0; i < maxiter && calcRelativeDiscrepancy() >= E; ++i) {

		real pp = calcScalarProduct(p, p);
		real alpha = calcScalarProduct(p, r) / pp;
		for (int j = 0; j < n; ++j) {
			x[j] += alpha * z[j];
			r[j] -= alpha * p[j];
		}

		vector <real> tmp = multD(r);
		tmp = multA(tmp);
		tmp = multD(tmp);
		real beta = -calcScalarProduct(p, tmp) / pp;
		for (int j = 0; j < n; ++j)
			p[j] = tmp[j] + beta * p[j];

		tmp = multD(r);
		for (int j = 0; j < n; ++j)
			z[j] = tmp[j] + beta * z[j];

		if (x == xprev)
			break;
		xprev = x;
	}


	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	fout << "LOS + diag" << "\t" << i << "\t" << elapsed_ms.count() << "\t" << calcRelativeDiscrepancy() << endl;
}





// �������� - ����������� ����� � �������� ������������� LU(sq)
void SLAE::LOSfactLUsq(std::ofstream& fout) {

	int i;
	x.clear();						// ����� ��������� �����������
	x.resize(n, 0);					// x_0 = (0, 0, ...)
	vector <real> xprev = x;
	decomposionLUsq();

	r = multA(x);					// r = A*x_0
	for (int i = 0; i < n; ++i)		// r = f - A*x_0
		r[i] = F[i] - r[i];
	r = execDirectTraversal(r);		// r = L^-1 (f - A*x_0)

	z = execReverseTraversal(r);	// z = U^-1 r

	p = multA(z);					// p = A*z
	p = execDirectTraversal(p);		// p = L^-1 A*z


	auto begin = std::chrono::steady_clock::now();
	for (i = 0; i < maxiter && calcRelativeDiscrepancy() >= E; ++i) {

		real pp = calcScalarProduct(p, p);
		real alpha = calcScalarProduct(p, r) / pp;
		for (int j = 0; j < n; ++j) {
			x[j] += alpha * z[j];
			r[j] -= alpha * p[j];
		}

		vector <real> tmp = execReverseTraversal(r);
		tmp = multA(tmp);
		tmp = execDirectTraversal(tmp);
		real beta = -calcScalarProduct(p, tmp) / pp;

		for (int j = 0; j < n; ++j)
			p[j] = tmp[j] + beta * p[j];

		tmp = execReverseTraversal(r);
		for (int j = 0; j < n; ++j)
			z[j] = tmp[j] + beta * z[j];

		if (x == xprev)
			break;
		xprev = x;
	}


	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	fout << "LOS + LU(sq)" << "\t" << i << "\t" << elapsed_ms.count() << "\t" << calcRelativeDiscrepancy() << endl;
}
