#include "eigen-runtime.h"

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

} // end extern "C"

