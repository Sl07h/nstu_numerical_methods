#include "head.h"


// Vectors F, x,  y = L\F,  x = L'\y
class vect {

public:
	
	void getVectX(vector <real> &x) { x = F; };
	int readVectFromFile(std::ifstream& fin, int size);
	void generateVectX(int size);
	void writeVectToFile(std::ofstream& fout, char *str);
	void writexCompError(std::ofstream& fout, char *str);
	void writeTableToFile(std::ofstream& fout);
	bool isXcorrect();


protected:
	vector <real> F;
};