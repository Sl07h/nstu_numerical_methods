#include "slae.h"
#include <iomanip>

int main() {

	
	std::ofstream fout00, fout11, fout12, fout21, fout22;
	
	//std::numeric_limits<real>::digits10 + 1
	
	SLAE slae;
	cout << std::fixed << std::setprecision(4);
	fout00 << std::fixed << std::setprecision(4);
	fout11 << std::fixed << std::setprecision(4);
	fout12 << std::fixed << std::setprecision(4);
	fout21 << std::fixed << std::setprecision(4);
	fout22 << std::fixed << std::setprecision(4);
	
	fout00.open("debugging.txt");
 	fout11.open("table.txt");
	fout12.open("x.txt");
	fout21.open("table_Hilbert.txt");
	fout22.open("x_Hilbert.txt");
	
	fout11 << "Method\tIterations count\tTime\tRelative discrepancy" << endl;
	//string folderName = "A"; 
	//string folderName = "B"; 
	string folderName = "0945"; 
	//string folderName = "4545";
	
	
	slae.readSLAEfromFiles(folderName);
	slae.LOS(fout11);
	slae.writeXToStream(fout12);
	slae.clearAll();
/*
	slae.readSLAEfromFiles(folderName);
	slae.LOSfactD(fout11);
	slae.writeXToStream(fout12);
	slae.clearAll();
	
	slae.readSLAEfromFiles(folderName);
	slae.LOSfactLUsq(fout11);
	slae.writeXToStream(fout12);
	slae.clearAll();*/
	

	/*
	int dimention = 9;
	slae.setE(1e-12);
	slae.setMaxiter(10000);

	slae.createHilbertMatrix(dimention);
	slae.generateVectX(dimention);
	slae.multAAndX();
	slae.LOS(fout21);
	slae.writeXToStream(fout22);
	slae.clearAll();

	slae.createHilbertMatrix(dimention);
	slae.generateVectX(dimention);
	slae.multAAndX();
	slae.LOSfactD(fout21);
	slae.writeXToStream(fout22);
	slae.clearAll();
	
	slae.createHilbertMatrix(dimention);
	slae.generateVectX(dimention);
	slae.multAAndX();
	slae.LOSfactLUsq(fout21);
	slae.writeXToStream(fout22);
	slae.clearAll();
	
	*/

	fout00.close();
	fout11.close();
	fout12.close();
	fout21.close();
	fout22.close();
}