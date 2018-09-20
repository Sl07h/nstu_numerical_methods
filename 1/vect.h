#include "head.h"


// Vectors F,  y = L\F,  x = L'\y
class vect {

public:
	int readFFromFile(std::ifstream& fin, int size);
	void generateVectX(int size);
	void writexToFile(std::ofstream& fout);


protected:
	vector <real> F;
};