#include "head.h"
#include"matrix.h"
#include "vect.h"


class SLAE : public matrix, public vect {

public:

	void convMatrixToDense();
	void writeDenseMatrixToFile(char* fileName);

	real multLine(vector<real>& line, int i, int mode);
	void mult();

	void Jacobi(real w);
	void GaussSeildel(real w);
	int calcIterative(int mode, real w);
	real findOptimalW(int mode);

	/* 4 ��������� ������ ������ ���������� ��������� ��������� 
	�� ���� ����������� ��� ���������������� ������ � ���������	*/
	vector <real_sum>& sumAtlowerTrianByDiag(vector <real>& line);
	vector <real_sum>& sumAtUpperTrianByDiag(vector <real>& line);
	void BoostedJacobi(real w);
	void BoostedGaussSeildel(real w);
	

protected:
	real calcNormE(vector <real> &x);
	real calcRelativeDiscrepancy();
	int calcCondNumber();
	vector <vector <real>> A;

};