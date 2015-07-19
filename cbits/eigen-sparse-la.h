#include "eigen-runtime.h"

#include "Eigen/Core"
#include "Eigen/SparseCore"
#include "Eigen/OrderingMethods"

#include "Eigen/src/Core/util/DisableStupidWarnings.h"

#include "Eigen/src/misc/Solve.h"
#include "Eigen/src/misc/SparseSolve.h"

#include "Eigen/src/IterativeLinearSolvers/IterativeSolverBase.h"
#include "Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h"
#include "Eigen/src/IterativeLinearSolvers/ConjugateGradient.h"
#include "Eigen/src/IterativeLinearSolvers/BiCGSTAB.h"
// #include "Eigen/src/IterativeLinearSolvers/IncompleteLUT.h" is not MPL2 compliant

#include "Eigen/src/Core/util/ReenableStupidWarnings.h"

#include "Eigen/SparseLU"
#include "Eigen/SparseQR"

using namespace Eigen;

extern "C" {

enum Solver {
	ConjugateGradient,
	BiCGSTAB,
	SparseLU,
	SparseQR
};

const char* eigen_sparse_la_newSolver(int, int, void**);
const char* eigen_sparse_la_freeSolver(int, int, void*);
const char* eigen_sparse_la_factorize(int, int, void*, void*);
const char* eigen_sparse_la_analyzePattern(int, int, void*, void*);
const char* eigen_sparse_la_compute(int, int, void*, void*);
const char* eigen_sparse_la_tolerance(int, int, void*, double*);
const char* eigen_sparse_la_setTolerance(int, int, void*, double);
const char* eigen_sparse_la_maxIterations(int, int, void*, int*);
const char* eigen_sparse_la_setMaxIterations(int, int, void*, int);
const char* eigen_sparse_la_info(int, int, void*, int*);
const char* eigen_sparse_la_error(int, int, void*, double*);
const char* eigen_sparse_la_iterations(int, int, void*, int*);
const char* eigen_sparse_la_solve(int, int, void*, void*, void**);
// const char* eigen_sparse_la_solveWithGuess(int, int, void*, void*, void*, void**);

} // end extern "C"

