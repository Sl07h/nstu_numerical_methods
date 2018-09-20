#include "head.h"
#include "matrix.h"
#include "vect.h"


// System of linear algebraic equations
class SLAE : public matrix, public vect {

public:
	void execDirectTraversal();
	void execReverseTraversal();

	void convToDense();
	void createHilbertMatrix();


private:
	vector <vector <real>> A;
};