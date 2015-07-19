#include "eigen-runtime.h"
#include <Eigen/SparseCore>

using namespace Eigen;

extern "C" {

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

} // end extern "C"

