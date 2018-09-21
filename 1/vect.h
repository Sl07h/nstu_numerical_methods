#include "head.h"


// Vectors F,  y = L\F,  x = L'\y
class vect {

public:
	
	~vect() { F.clear(); }
	void getVectX(vector <real> &x) { x = F; };
	int readFFromFile(std::ifstream& fin, int size);
	void generateVectX(int size);
	void writexToFile(std::ofstream& fout);
	bool isXcorrect();


protected:
	vector <real> F;
};