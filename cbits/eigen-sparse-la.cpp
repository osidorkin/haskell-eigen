#include "eigen-sparse-la.h"

using namespace Eigen;

typedef SparseMatrix <T0> M0;
typedef SparseMatrix <T1> M1;
typedef SparseMatrix <T2> M2;
typedef SparseMatrix <T3> M3;

template <class T> struct S0 { typedef ConjugateGradient< SparseMatrix <T>, Lower, DiagonalPreconditioner<typename SparseMatrix <T>::Scalar> > type; };
template <class T> struct S1 { typedef ConjugateGradient< SparseMatrix <T>, Lower, IdentityPreconditioner > type; };
template <class T> struct S2 { typedef BiCGSTAB< SparseMatrix <T>, DiagonalPreconditioner<typename SparseMatrix <T>::Scalar> > type; };
template <class T> struct S3 { typedef BiCGSTAB< SparseMatrix <T>, IdentityPreconditioner > type; };
template <class T> struct S4 { typedef SparseLU< SparseMatrix <T>, NaturalOrdering<int> > type; };
template <class T> struct S5 { typedef SparseLU< SparseMatrix <T>, COLAMDOrdering<int> > type; };
template <class T> struct S6 { typedef SparseQR< SparseMatrix <T>, NaturalOrdering<int> > type; };
template <class T> struct S7 { typedef SparseQR< SparseMatrix <T>, COLAMDOrdering<int> > type; };

#define API_ALL(name,args,call) \
extern "C" RET eigen_##name args {\
    GUARD_START\
    switch (code) {\
        case 0: switch(s) {\
            case 0: return name<T0, M0, S0<T0>::type >call;\
            case 1: return name<T0, M0, S1<T0>::type >call;\
            case 2: return name<T0, M0, S2<T0>::type >call;\
            case 3: return name<T0, M0, S3<T0>::type >call;\
            case 4: return name<T0, M0, S4<T0>::type >call;\
            case 5: return name<T0, M0, S5<T0>::type >call;\
            case 6: return name<T0, M0, S6<T0>::type >call;\
            case 7: return name<T0, M0, S7<T0>::type >call;\
        }\
        case 1: switch(s) {\
            case 0: return name<T1, M1 ,S0<T1>::type >call;\
            case 1: return name<T1, M1 ,S1<T1>::type >call;\
            case 2: return name<T1, M1 ,S2<T1>::type >call;\
            case 3: return name<T1, M1 ,S3<T1>::type >call;\
            case 4: return name<T1, M1 ,S4<T1>::type >call;\
            case 5: return name<T1, M1 ,S5<T1>::type >call;\
            case 6: return name<T1, M1 ,S6<T1>::type >call;\
            case 7: return name<T1, M1 ,S7<T1>::type >call;\
        }\
        case 2: switch(s) {\
            case 0: return name<T2, M2 ,S0<T2>::type >call;\
            case 1: return name<T2, M2 ,S1<T2>::type >call;\
            case 2: return name<T2, M2 ,S2<T2>::type >call;\
            case 3: return name<T2, M2 ,S3<T2>::type >call;\
            case 4: return name<T2, M2 ,S4<T2>::type >call;\
            case 5: return name<T2, M2 ,S5<T2>::type >call;\
            case 6: return name<T2, M2 ,S6<T2>::type >call;\
            case 7: return name<T2, M2 ,S7<T2>::type >call;\
        }\
        case 3: switch(s) {\
            case 0: return name<T3, M3 ,S0<T3>::type >call;\
            case 1: return name<T3, M3 ,S1<T3>::type >call;\
            case 2: return name<T3, M3 ,S2<T3>::type >call;\
            case 3: return name<T3, M3 ,S3<T3>::type >call;\
            case 4: return name<T3, M3 ,S4<T3>::type >call;\
            case 5: return name<T3, M3 ,S5<T3>::type >call;\
            case 6: return name<T3, M3 ,S6<T3>::type >call;\
            case 7: return name<T3, M3 ,S7<T3>::type >call;\
        }\
    }\
    GUARD_END\
}

#define API_ITERATIVE(name,args,call) \
extern "C" RET eigen_##name args {\
    GUARD_START\
    switch (code) {\
        case 0: switch(s) {\
            case 0: return name<T0, M0, S0<T0>::type >call;\
            case 1: return name<T0, M0, S1<T0>::type >call;\
            case 2: return name<T0, M0, S2<T0>::type >call;\
            case 3: return name<T0, M0, S3<T0>::type >call;\
        }\
        case 1: switch(s) {\
            case 0: return name<T1, M1, S0<T1>::type >call;\
            case 1: return name<T1, M1, S1<T1>::type >call;\
            case 2: return name<T1, M1, S2<T1>::type >call;\
            case 3: return name<T1, M1, S3<T1>::type >call;\
        }\
        case 2: switch(s) {\
            case 0: return name<T2, M2, S0<T2>::type >call;\
            case 1: return name<T2, M2, S1<T2>::type >call;\
            case 2: return name<T2, M2, S2<T2>::type >call;\
            case 3: return name<T2, M2, S3<T2>::type >call;\
        }\
        case 3: switch(s) {\
            case 0: return name<T3, M3, S0<T3>::type >call;\
            case 1: return name<T3, M3, S1<T3>::type >call;\
            case 2: return name<T3, M3, S2<T3>::type >call;\
            case 3: return name<T3, M3, S3<T3>::type >call;\
        }\
    }\
    return strdup("supported for iterative solver only");\
    GUARD_END\
}

template <class T, class M, class S>
RET sparse_la_newSolver(void** p) {
    *(S**)p = new S;
    return 0;
}
API_ALL(sparse_la_newSolver, (int code, int s, void** p), (p));

template <class T, class M, class S>
RET sparse_la_freeSolver(void* p) {
    delete (S*)p;
    return 0;
}
API_ALL(sparse_la_freeSolver, (int code, int s, void* p), (p));

template <class T, class M, class S>
RET sparse_la_factorize(void* p, void* a) {
    ((S*)p)->factorize(*(M*)a);
    return 0;
}
API_ALL(sparse_la_factorize, (int code, int s, void* p, void* a), (p,a));

template <class T, class M, class S>
RET sparse_la_analyzePattern(void* p, void* a) {
    ((S*)p)->analyzePattern(*(M*)a);
    return 0;
}
API_ALL(sparse_la_analyzePattern, (int code, int s, void* p, void* a), (p,a));

template <class T, class M, class S>
RET sparse_la_compute(void* p, void* a) {
    ((S*)p)->compute(*(M*)a);
    return 0;
}
API_ALL(sparse_la_compute, (int code, int s, void* p, void* a), (p,a));

template <class T, class M, class S>
RET sparse_la_tolerance(void* p, double* x) {
    *x = ((S*)p)->tolerance();
    return 0;
}
API_ITERATIVE(sparse_la_tolerance, (int code, int s, void* p, double* x), (p,x));

template <class T, class M, class S>
RET sparse_la_setTolerance(void* p, double x) {
    ((S*)p)->setTolerance(x);
    return 0;
}
API_ITERATIVE(sparse_la_setTolerance, (int code, int s, void* p, double x), (p,x));

template <class T, class M, class S>
RET sparse_la_maxIterations(void* p, int* x) {
    *x = ((S*)p)->maxIterations();
    return 0;
}
API_ITERATIVE(sparse_la_maxIterations, (int code, int s, void* p, int* x), (p,x));

template <class T, class M, class S>
RET sparse_la_setMaxIterations(void* p, int x) {
    ((S*)p)->setMaxIterations(x);
    return 0;
}
API_ITERATIVE(sparse_la_setMaxIterations, (int code, int s, void* p, int x), (p,x));

template <class T, class M, class S>
RET sparse_la_info(void* p, int* x) {
    *x = ((S*)p)->info();
    return 0;
}
API_ALL(sparse_la_info, (int code, int s, void* p, int* x), (p,x));

template <class T, class M, class S>
RET sparse_la_error(void* p, double* x) {
    *x = ((S*)p)->error();
    return 0;
}
API_ITERATIVE(sparse_la_error, (int code, int s, void* p, double* x), (p,x));

template <class T, class M, class S>
RET sparse_la_iterations(void* p, int* x) {
    *x = ((S*)p)->iterations();
    return 0;
}
API_ITERATIVE(sparse_la_iterations, (int code, int s, void* p, int* x), (p,x));

template <class T, class M, class S>
RET sparse_la_solve(void* p, void* b, void** x) {
    *x = new M(((S*)p)->solve(*(M*)b));
    return 0;
}
API_ALL(sparse_la_solve, (int code, int s, void* p, void* b, void** x), (p,b,x));

// template <class T, class M, class S>
// RET sparse_la_solveWithGuess(void* p, void* b, void* x0, void** x) {
//  *x = new M(((S*)p)->solveWithGuess(Matrix<T,Dynamic,Dynamic>(*(M*)b), Matrix<T,Dynamic,Dynamic>(*(M*)x0)));
//  return 0;
// }
// API_ALL(sparse_la_solveWithGuess, (int code, int s, void* p, void* b, void* x0, void** x), (p,b,x0,x));

#define API_SPARSE_LU(name,args,call) \
extern "C" RET eigen_##name args {\
    GUARD_START\
    switch (code) {\
        case 0: switch(s) {\
            case 4: return name<T0, M0, S4<T0>::type >call;\
            case 5: return name<T0, M0, S5<T0>::type >call;\
        }\
        case 1: switch(s) {\
            case 4: return name<T1, M1, S4<T1>::type >call;\
            case 5: return name<T1, M1, S5<T1>::type >call;\
        }\
        case 2: switch(s) {\
            case 4: return name<T2, M2, S4<T2>::type >call;\
            case 5: return name<T2, M2, S5<T2>::type >call;\
        }\
        case 3: switch(s) {\
            case 4: return name<T3, M3, S4<T3>::type >call;\
            case 5: return name<T3, M3, S5<T3>::type >call;\
        }\
    }\
    return strdup("supported for SparseLU solver only");\
    GUARD_END\
}

#define API_SPARSE_LU_NO_COMPLEX(name,args,call) \
extern "C" RET eigen_##name args {\
    GUARD_START\
    switch (code) {\
        case 0: switch(s) {\
            case 4: return name<T0, M0, S4<T0>::type >call;\
            case 5: return name<T0, M0, S5<T0>::type >call;\
        }\
        case 1: switch(s) {\
            case 4: return name<T1, M1, S4<T1>::type >call;\
            case 5: return name<T1, M1, S5<T1>::type >call;\
        }\
    }\
    return strdup("supported for SparseLU solver and non-complex matrix only");\
    GUARD_END\
}

#define API_SPARSE_QR(name,args,call) \
extern "C" RET eigen_##name args {\
    GUARD_START\
    switch (code) {\
        case 0: switch(s) {\
            case 6: return name<T0, M0, S6<T0>::type >call;\
            case 7: return name<T0, M0, S7<T0>::type >call;\
        }\
        case 1: switch(s) {\
            case 6: return name<T1, M1, S6<T1>::type >call;\
            case 7: return name<T1, M1, S7<T1>::type >call;\
        }\
        case 2: switch(s) {\
            case 6: return name<T2, M2, S6<T2>::type >call;\
            case 7: return name<T2, M2, S7<T2>::type >call;\
        }\
        case 3: switch(s) {\
            case 6: return name<T3, M3, S6<T3>::type >call;\
            case 7: return name<T3, M3, S7<T3>::type >call;\
        }\
    }\
    return strdup("supported for SparseQR solver only");\
    GUARD_END\
}

template <class T, class M, class S>
RET sparse_la_matrixQ(void* p, void** x) {
    // *x = new M(((S*)p)->matrixQ());
    // return 0;
    return strdup("not implemented yet");
}
API_SPARSE_QR(sparse_la_matrixQ, (int code, int s, void* p, void** x), (p,x));

template <class T, class M, class S>
RET sparse_la_matrixR(void* p, void** x) {
    *x = new M(((S*)p)->matrixR());
    return 0;
}
API_SPARSE_QR(sparse_la_matrixR, (int code, int s, void* p, void** x), (p,x));

template <class T, class M, class S>
RET sparse_la_setPivotThreshold(void* p, double x) {
    ((S*)p)->setPivotThreshold(x);
    return 0;
}
API_SPARSE_QR(sparse_la_setPivotThreshold, (int code, int s, void* p, double x), (p,x));

template <class T, class M, class S>
RET sparse_la_rank(void* p, double* x) {
    *x = ((S*)p)->rank();
    return 0;
}
API_SPARSE_QR(sparse_la_rank, (int code, int s, void* p, double* x), (p,x));

template <class T, class M, class S>
RET sparse_la_setSymmetric(void* p, int x) {
    ((S*)p)->isSymmetric(x);
    return 0;
}
API_SPARSE_LU(sparse_la_setSymmetric, (int code, int s, void* p, int x), (p,x));

template <class T, class M, class S>
RET sparse_la_matrixL(void* p, void** x) {
    // *x = new M(((S*)p)->matrixL());
    // return 0;
    return strdup("not implemented yet");
}
API_SPARSE_LU(sparse_la_matrixL, (int code, int s, void* p, void** x), (p,x));

template <class T, class M, class S>
RET sparse_la_matrixU(void* p, void** x) {
    // *x = new M(((S*)p)->matrixU());
    // return 0;
    return strdup("not implemented yet");
}
API_SPARSE_LU(sparse_la_matrixU, (int code, int s, void* p, void** x), (p,x));

template <class T, class M, class S>
RET sparse_la_determinant(void* p, void* x) {
    *(T*)x = ((S*)p)->determinant();
    return 0;
}
API_SPARSE_LU(sparse_la_determinant, (int code, int s, void* p, void* x), (p,x));

template <class T, class M, class S>
RET sparse_la_logAbsDeterminant(void* p, void* x) {
    *(T*)x = ((S*)p)->logAbsDeterminant();
    return 0;
}
API_SPARSE_LU(sparse_la_logAbsDeterminant, (int code, int s, void* p, void* x), (p,x));

template <class T, class M, class S>
RET sparse_la_absDeterminant(void* p, void* x) {
    *(T*)x = ((S*)p)->absDeterminant();
    return 0;
}
API_SPARSE_LU(sparse_la_absDeterminant, (int code, int s, void* p, void* x), (p,x));

template <class T, class M, class S>
RET sparse_la_signDeterminant(void* p, void* x) {
    *(T*)x = ((S*)p)->signDeterminant();
    return 0;
}
API_SPARSE_LU_NO_COMPLEX(sparse_la_signDeterminant, (int code, int s, void* p, void* x), (p,x));
