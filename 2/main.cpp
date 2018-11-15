#include "head.h"
#include "slae.h"


void main() {

	std::ofstream fout_res;
	SLAE slae;


	/*slae.generateMatrixWith7Diagonals(10, 3);
	slae.setE(1e-10);
	slae.setMaxiter(200000);
	slae.writeMatrixToFile("A.txt");
	slae.convMatrixToDense();
	slae.writeDenseMatrixToFile("A_dense.txt");

	slae.invertSigns();
	slae.writeMatrixToFile("B.txt");
	slae.convMatrixToDense();
	slae.writeDenseMatrixToFile("B_dense.txt");*/



	fout_res.open("A_Jacobi.txt");
	slae.readMatrixFromFile("A.txt");
	slae.generateVectX(slae.getDimention());
	slae.mult();
	slae.writeFToFile("F_A.txt");
	cout << slae.findOptimalW(1, fout_res) << endl;
	fout_res.close();


	fout_res.open("A_GaussSeidel.txt");
	slae.generateVectX(slae.getDimention());
	slae.mult();
	cout << slae.findOptimalW(2, fout_res) << endl;
	fout_res.close();



	fout_res.open("B_Jacobi.txt");
	slae.readMatrixFromFile("B.txt");
	slae.generateVectX(slae.getDimention());
	slae.mult();
	slae.writeFToFile("F_B.txt");
	cout << slae.findOptimalW(1, fout_res) << endl;
	fout_res.close();


	fout_res.open("B_GaussSeidel.txt");
	slae.generateVectX(slae.getDimention());
	slae.mult();
	cout << slae.findOptimalW(2, fout_res) << endl;
	fout_res.close();
}
