#define EIGEN_MPL2_ONLY
#define EIGEN2_SUPPORT
#define EIGEN_NO_EIGEN2_DEPRECATED_WARNING
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
void eigen_assert_fail(const char* condition, const char* function, const char* file, int line);
#define eigen_assert(x) do {\
    if (!(x)) eigen_assert_fail(#x, __PRETTY_FUNCTION__, __FILE__, __LINE__);\
} while(false)
#include <Eigen/Core>

using namespace Eigen;

extern "C" {

const char* eigen_add(double*, int, int, const double*, int, int, const double*, int, int);
const char* eigen_sub(double*, int, int, const double*, int, int, const double*, int, int);
const char* eigen_mul(double*, int, int, const double*, int, int, const double*, int, int);
double eigen_norm(const double*, int, int);
double eigen_squaredNorm(const double*, int, int);
double eigen_blueNorm(const double*, int, int);
double eigen_hypotNorm(const double*, int, int);
double eigen_sum(const double*, int, int);
double eigen_prod(const double*, int, int);
double eigen_mean(const double*, int, int);
double eigen_trace(const double*, int, int);
double eigen_determinant(const double*, int, int);
const char* eigen_transpose(double*, int, int, const double*, int, int);
const char* eigen_normalize(double*, int, int);
const char* eigen_inverse(double*, int, int, const double*, int, int);
const char* eigen_conjugate(double*, int, int, const double*, int, int);
const char* eigen_adjoint(double*, int, int, const double*, int, int);
bool eigen_initParallel();
int eigen_getNbThreads();
void eigen_setNbThreads(int);

enum Decomposition {
	PartialPivLU, FullPivLU, HouseholderQR, ColPivHouseholderQR, FullPivHouseholderQR, LLT, LDLT, JacobiSVD
};

const char* eigen_solve(Decomposition d, // Ax=b 
	double*, int, int, // x
	const double*, int, int, // A
	const double*, int, int); // b

const char* eigen_relativeError(double& e, 
	const double*, int, int, // x
	const double*, int, int, // A
	const double*, int, int); // b;


} // end extern "C"

#define GUARD_START try { assert(inited); do {
#define GUARD_END } while(false); return 0; } catch (const eigen_assert_exception& ex) { return strdup(ex.what()); }
