#define EIGEN_MPL2_ONLY
#define EIGEN2_SUPPORT
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
void eigen_assert_fail(const char* condition, const char* function, const char* file, int line);
#define eigen_assert(x) do {\
    if (!(x)) eigen_assert_fail(#x, __PRETTY_FUNCTION__, __FILE__, __LINE__);\
} while(false)
#include <Eigen/Core>

using namespace Eigen;

extern "C" {

MatrixXd* eigen_create(int rows, int cols);
void eigen_destroy(const MatrixXd*);
MatrixXd* eigen_clone(const MatrixXd& m);
const char* eigen_set(MatrixXd&,int,int,double);
const char* eigen_get(double&,const MatrixXd&,int,int);
int eigen_rows(const MatrixXd&);
int eigen_cols(const MatrixXd&);
double* eigen_data(MatrixXd&);
void eigen_copy(MatrixXd&, const MatrixXd&);
void eigen_resize(MatrixXd&, int, int);
const char* eigen_add(MatrixXd&, const MatrixXd&, const MatrixXd&);
const char* eigen_sub(MatrixXd&, const MatrixXd&, const MatrixXd&);
const char* eigen_mul(MatrixXd&, const MatrixXd&, const MatrixXd&);
double eigen_norm(const MatrixXd&);
double eigen_squaredNorm(const MatrixXd&);
double eigen_blueNorm(const MatrixXd&);
double eigen_stableNorm(const MatrixXd&);
double eigen_hypotNorm(const MatrixXd&);
double eigen_sum(const MatrixXd&);
double eigen_prod(const MatrixXd&);
double eigen_mean(const MatrixXd&);
double eigen_minCoeff(const MatrixXd&);
double eigen_maxCoeff(const MatrixXd&);
double eigen_trace(const MatrixXd&);
const char* eigen_transpose(MatrixXd&);
const char* eigen_normalize(MatrixXd&);
const char* eigen_inverse(MatrixXd&);
const char* eigen_adjoint(MatrixXd&);
double eigen_determinant(const MatrixXd& src);
void eigen_initParallel();
void eigen_setNbThreads(int);

enum Decomposition {
	PartialPivLU, FullPivLU, HouseholderQR, ColPivHouseholderQR, FullPivHouseholderQR, LLT, LDLT, JacobiSVD
};

const char* eigen_solve(Decomposition d, MatrixXd& x, const MatrixXd& A, const MatrixXd& b); // Ax=b
const char* eigen_relativeError(double& e, const MatrixXd& x, const MatrixXd& A, const MatrixXd& b);


} // end extern "C"

#define GUARD_START try { do {
#define GUARD_END } while(false); return 0; } catch (const std::exception& ex) { return strdup(ex.what()); }
