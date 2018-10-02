#include "head.h"
#include "profMatrix.h"
#include "vect.h"


// System of linear algebraic equations
class SLAE : public matrix, public vect {

public:

	void writeMatrixtoFile(std::ofstream& fout, char* str);
	
	void execDirectTraversal();
	void execReverseTraversal();

	void convAToDense();
	void convLToDense();
	void mult();
	void createHilbertMatrix(int n);


private:
	vector <vector <double>> A;
};