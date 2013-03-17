#include "eigen-proxy.h"
#include <Eigen/LU>
#include <Eigen/LeastSquares>
#include <stdio.h>
#include <sstream>

struct eigen_assert_exception : std::exception {
	std::string _what;
	eigen_assert_exception(const std::string& what) : _what(what) {}
	~eigen_assert_exception() throw () {}
	const char* what() const throw () { return _what.c_str(); }
};

void eigen_assert_fail(const char* condition, const char* function, const char* file, int line) {
	std::ostringstream os;
    os << "assertion failed: " << condition << " in function " << function << " at " << file << ":" << line << std::endl;
    throw eigen_assert_exception(os.str());
}

MatrixXd* eigen_create(int rows, int cols) { return new MatrixXd(rows, cols); }
MatrixXd* eigen_clone(const MatrixXd& m) { return new MatrixXd(m); }
void eigen_destroy(const MatrixXd* m) { delete m; }
const char* eigen_set(MatrixXd& dst, int row, int col, double val) { GUARD_START dst(row,col) = val; GUARD_END }
const char* eigen_get(double&r, const MatrixXd& m, int row, int col) { GUARD_START r = m(row,col); GUARD_END }
int eigen_rows(const MatrixXd& m) { return m.rows(); } 
int eigen_cols(const MatrixXd& m) { return m.cols(); }
double* eigen_data(MatrixXd& m) { return m.data(); } 
void eigen_copy(MatrixXd& dst, const MatrixXd& src) { dst = src; }
void eigen_resize(MatrixXd& m, int rows, int cols) { m.resize(rows, cols); }
const char* eigen_add(MatrixXd& ret, const MatrixXd& lhs, const MatrixXd& rhs) { GUARD_START ret = lhs + rhs; GUARD_END }
const char* eigen_sub(MatrixXd& ret, const MatrixXd& lhs, const MatrixXd& rhs) { GUARD_START ret = lhs - rhs; GUARD_END }
const char* eigen_mul(MatrixXd& ret, const MatrixXd& lhs, const MatrixXd& rhs) { GUARD_START ret = lhs * rhs; GUARD_END }
double eigen_norm(const MatrixXd& m) { return m.norm(); }
double eigen_squaredNorm(const MatrixXd& m) { return m.squaredNorm(); }
double eigen_blueNorm(const MatrixXd& m) { return m.squaredNorm(); }
// double mstableNorm(const MatrixXd& m) { return m.stableNorm(); }
double eigen_hypotNorm(const MatrixXd& m) { return m.hypotNorm(); }
double eigen_sum(const MatrixXd& m) { return m.sum(); }
double eigen_prod(const MatrixXd& m) { return m.prod(); }
double eigen_mean(const MatrixXd& m) { return m.mean(); }
double eigen_minCoeff(const MatrixXd& m) { return m.minCoeff(); }
double eigen_maxCoeff(const MatrixXd& m) { return m.maxCoeff(); }
double eigen_trace(const MatrixXd& m) { return m.trace(); }

const char* eigen_inverse(MatrixXd& m) { GUARD_START m = m.inverse().eval(); GUARD_END }
const char* eigen_adjoint(MatrixXd& m) { GUARD_START m.adjointInPlace(); GUARD_END }
const char* eigen_transpose(MatrixXd& m) { GUARD_START m.transposeInPlace(); GUARD_END }
const char* eigen_normalize(MatrixXd& m) { GUARD_START m.normalize(); GUARD_END }
const char* eigen_conjugate(MatrixXd& m) { GUARD_START m = m.conjugate().eval(); GUARD_END }


double eigen_determinant(const MatrixXd& m) { return m.determinant(); }

const char* eigen_solve(Decomposition d, MatrixXd& x, const MatrixXd& A, const MatrixXd& b) {
	GUARD_START
	switch (d) {
		case ::PartialPivLU:
			x = A.partialPivLu().solve(b); 
			break;
		case ::FullPivLU:
			x = A.fullPivLu().solve(b); 
			break;
		case ::HouseholderQR:
			x = A.householderQr().solve(b); 
			break;
		case ::ColPivHouseholderQR:
			x = A.colPivHouseholderQr().solve(b); 
			break;
		case ::FullPivHouseholderQR: 
			x = A.fullPivHouseholderQr().solve(b); 
			break;
		case ::LLT:
			x = A.llt().solve(b); 
			break;
		case ::LDLT:
			x = A.ldlt().solve(b); 
			break;
		case ::JacobiSVD:
			x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b); 
			break;
	}
	GUARD_END
}

const char* eigen_relativeError(double& e, const MatrixXd& x, const MatrixXd& A, const MatrixXd& b) {
	GUARD_START
	e = (A*x - b).norm() / b.norm();
	GUARD_END
}

void eigen_initParallel() {
	initParallel();
}

void eigen_setNbThreads(int n) {
	setNbThreads(n);
}

