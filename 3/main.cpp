#include "slae.h"
#include <iomanip>


void testSLAE(const string &folderName, bool firstNumberIsOne, bool doWriteHeader) {
	SLAE slae;
	pair <int, real> iterationsCountAndDiscrapancy;

	std::ofstream foutTable, foutX;
	foutTable.open(folderName + "/table.txt");
	foutX.open(folderName + "/x.txt");
	foutTable << std::fixed << std::setprecision(std::numeric_limits<real>::digits10 + 1) << std::scientific;
	foutX << std::fixed << std::setprecision(std::numeric_limits<real>::digits10 + 1);

	if (doWriteHeader)
		foutTable << "Method\tIterations count\tTime\tRelative discrepancy" << endl;


	slae.clearAll();
	slae.readSLAEfromFiles(folderName, firstNumberIsOne);
	auto begin = std::chrono::steady_clock::now();
	iterationsCountAndDiscrapancy = slae.LOS();
	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
	foutTable << "LOS\t" << iterationsCountAndDiscrapancy.first << "\t" << elapsed_ms.count() << "\t" << iterationsCountAndDiscrapancy.second << endl;
	slae.writeXToStream(foutX);



	slae.clearAll();
	slae.readSLAEfromFiles(folderName, firstNumberIsOne);
	begin = std::chrono::steady_clock::now();
	iterationsCountAndDiscrapancy = slae.LOSfactD();
	end = std::chrono::steady_clock::now();
	elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
	foutTable << "LOS + diag\t" << iterationsCountAndDiscrapancy.first << "\t" << elapsed_ms.count() << "\t" << iterationsCountAndDiscrapancy.second << endl;
	slae.writeXToStream(foutX);



	slae.clearAll();
	slae.readSLAEfromFiles(folderName, firstNumberIsOne);
	begin = std::chrono::steady_clock::now();
	iterationsCountAndDiscrapancy = slae.LOSfactLUsq();
	end = std::chrono::steady_clock::now();
	elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
	foutTable << "LOS + LU(sq)\t" << iterationsCountAndDiscrapancy.first << "\t" << elapsed_ms.count() << "\t" << iterationsCountAndDiscrapancy.second << endl;
	slae.writeXToStream(foutX);


	foutTable.close();
	foutX.close();
}




int main() {

	// Сначала нужно создать папки HilbertN, N - размерность матрицы
	SLAE slae;
	slae.createHilbertMatricies(4, 10, 3, "Hilbert");

	testSLAE("A", false, false);
	testSLAE("B", false, false);

	testSLAE("Hilbert4", false, false);
	testSLAE("Hilbert7", false, false);
	testSLAE("Hilbert10", false, false);

	testSLAE("0945", true, false);
	testSLAE("4545", true, false);
}