#define EIGEN_MPL2_ONLY
#define EIGEN2_SUPPORT
#define EIGEN_NO_EIGEN2_DEPRECATED_WARNING
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
void eigen_assert_fail(const char* condition, const char* function, const char* file, int line);
#define eigen_assert(x) do {\
    if (!(x)) eigen_assert_fail(#x, __PRETTY_FUNCTION__, __FILE__, __LINE__);\
} while(false)
#include <Eigen/Core>
#include <Eigen/SparseCore>

using namespace Eigen;

extern "C" {

const char* eigen_add(int, void*, int, int, const void*, int, int, const void*, int, int);
const char* eigen_sub(int, void*, int, int, const void*, int, int, const void*, int, int);
const char* eigen_mul(int, void*, int, int, const void*, int, int, const void*, int, int);
const char* eigen_norm(int, void*, const void*, int, int);
const char* eigen_squaredNorm(int, void*, const void*, int, int);
const char* eigen_blueNorm(int, void*, const void*, int, int);
const char* eigen_hypotNorm(int, void*, const void*, int, int);
const char* eigen_sum(int, void*, const void*, int, int);
const char* eigen_prod(int, void*, const void*, int, int);
const char* eigen_mean(int, void*, const void*, int, int);
const char* eigen_trace(int, void*, const void*, int, int);
const char* eigen_determinant(int, void*, const void*, int, int);
const char* eigen_diagonal(int, void*, int, int, const void*, int, int);
const char* eigen_transpose(int, void*, int, int, const void*, int, int);
const char* eigen_normalize(int, void*, int, int);
const char* eigen_random(int, void*, int, int);
const char* eigen_identity(int, void*, int, int);
const char* eigen_inverse(int, void*, int, int, const void*, int, int);
const char* eigen_conjugate(int, void*, int, int, const void*, int, int);
const char* eigen_adjoint(int, void*, int, int, const void*, int, int);
bool eigen_initParallel();
int eigen_getNbThreads();
void eigen_setNbThreads(int);

#define SPARSE_PROP(name) const char* eigen_sparse_##name(int, void*, void*)
SPARSE_PROP(cols);
SPARSE_PROP(rows);
SPARSE_PROP(innerSize);
SPARSE_PROP(outerSize);
SPARSE_PROP(nonZeros);
SPARSE_PROP(isCompressed);
SPARSE_PROP(norm);
SPARSE_PROP(squaredNorm);
SPARSE_PROP(blueNorm);
#undef SPARSE_PROP

#define SPARSE_BINOP(name) const char* eigen_sparse_##name(int, void*, void*, void**);
SPARSE_BINOP(add);
SPARSE_BINOP(sub);
SPARSE_BINOP(mul);
SPARSE_BINOP(prunedRef);
SPARSE_BINOP(scale);
#undef SPARSE_BINOP

#define SPARSE_UNOP(name) const char* eigen_sparse_##name(int, void*, void**);
SPARSE_UNOP(makeCompressed);
SPARSE_UNOP(uncompress);
SPARSE_UNOP(adjont);
SPARSE_UNOP(transponse);
SPARSE_UNOP(pruned);
#undef SPARSE_UNOP

const char* eigen_sparse_fromList(int, int, int, void*, int, void**);
const char* eigen_sparse_toList(int, void*, void*, int);
const char* eigen_sparse_free(int, void*);
const char* eigen_sparse_coeff(int, void*, int, int, void*);
const char* eigen_sparse_block(int, void*, int, int, int, int, void**);
const char* eigen_sparse_fromMatrix(int, void*, int, int, void**);
const char* eigen_sparse_toMatrix(int, void*, void*, int, int);

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

#define GUARD_START try { assert(inited); do {
#define GUARD_END } while(false); return 0; } catch (const std::exception& ex) { return strdup(ex.what()); }
