#include "head.h"
#include "matrix.h"
#include "vect.h"


// System of linear algebraic equations
class SLAE : public matrix, public vect {

public:
	~SLAE() { A.clear(); }
	void execDirectTraversal();
	void execReverseTraversal();

	void convToDense();
	void mult();
	void mult_dense();
	void createHilbertMatrix();


private:
	vector <vector <real>> A;
};