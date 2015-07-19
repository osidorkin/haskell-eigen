#include "eigen-runtime.h"

extern "C" {

enum Decomposition {
	PartialPivLU,
	FullPivLU,
	HouseholderQR,
	ColPivHouseholderQR,
	FullPivHouseholderQR,
	LLT,
	LDLT,
	JacobiSVD,
	// SelfAdjointEigenSolver,
	// ComplexEigenSolver,
	// EigenSolver,
	// GeneralizedSelfAdjointEigenSolver
};

const char* eigen_rank(int, Decomposition d, int*, const void*, int, int);
const char* eigen_kernel(int, Decomposition d, void**, int*, int*, const void*, int, int);
const char* eigen_image(int, Decomposition d, void**, int*, int*, const void*, int, int);

const char* eigen_solve(int, Decomposition d, // Ax=b 
	void*, int, int, // x
	const void*, int, int, // A
	const void*, int, int); // b


const char* eigen_relativeError(int, void* e, 
	const void*, int, int, // x
	const void*, int, int, // A
	const void*, int, int); // b;


} // end extern "C"

