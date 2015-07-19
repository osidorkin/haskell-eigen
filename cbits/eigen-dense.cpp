#include "eigen-dense.h"
#include <Eigen/Core>
#include <Eigen/LeastSquares>

using namespace Eigen;

template <class T>
Map< Matrix<T,Dynamic,Dynamic> > matrix(void* p, int r, int c) {
    return Map< Matrix<T,Dynamic,Dynamic> >((T*)p, r, c);
}

template <class T>
Map< Matrix<T,Dynamic,Dynamic> > matrix(const void* p, int r, int c) {
    return Map< Matrix<T,Dynamic,Dynamic> >((const T*)p, r, c);
}

#define RET const char*

#define API(name,args,call) \
extern "C" RET eigen_##name args {\
    GUARD_START\
    switch (code) {\
        case 0: return name<T0>call;\
        case 1: return name<T1>call;\
        case 2: return name<T2>call;\
        case 3: return name<T3>call;\
    }\
    GUARD_END\
}


#define BINOP(name,op) \
template <class T>\
RET name(void* p, int r, int c,\
    const void* p1, int r1, int c1,\
    const void* p2, int r2, int c2)\
{\
    matrix<T>(p,r,c) = matrix<T>(p1,r1,c1) op matrix<T>(p2,r2,c2);\
    return 0;\
}\
API(name, (int code,\
    void* p, int r, int c,\
    const void* p1, int r1, int c1,\
    const void* p2, int r2, int c2), (p,r,c,p1,r1,c1,p2,r2,c2));

BINOP(add,+);
BINOP(sub,-);
BINOP(mul,*);

#define PROP(name) \
extern "C" RET eigen_##name(int code, void* q, const void* p, int r, int c) {\
        GUARD_START\
        switch (code) {\
            case 0: *(T0*)q = matrix<T0>(p,r,c).name(); break;\
            case 1: *(T1*)q = matrix<T1>(p,r,c).name(); break;\
            case 2: *(T2*)q = matrix<T2>(p,r,c).name(); break;\
            case 3: *(T3*)q = matrix<T3>(p,r,c).name(); break;\
        }\
        GUARD_END\
    }

PROP(norm);
PROP(squaredNorm);
PROP(blueNorm);
PROP(hypotNorm);
PROP(sum);
PROP(prod);
PROP(mean);
PROP(trace);
PROP(determinant);

#define UNOP(name) \
extern "C" RET eigen_##name(int code, void* p, int r, int c, const void* p1, int r1, int c1) {\
        GUARD_START\
        switch (code) {\
            case 0: matrix<T0>(p,r,c) = matrix<T0>(p1,r1,c1).name(); break;\
            case 1: matrix<T1>(p,r,c) = matrix<T1>(p1,r1,c1).name(); break;\
            case 2: matrix<T2>(p,r,c) = matrix<T2>(p1,r1,c1).name(); break;\
            case 3: matrix<T3>(p,r,c) = matrix<T3>(p1,r1,c1).name(); break;\
        }\
        GUARD_END\
    }

UNOP(inverse);
UNOP(adjoint);
UNOP(conjugate);
UNOP(diagonal);
UNOP(transpose);

extern "C" RET eigen_normalize(int code, void* p, int r, int c)
{
    GUARD_START
    switch (code) {
        case 0: matrix<T0>(p,r,c).normalize(); break;
        case 1: matrix<T1>(p,r,c).normalize(); break;
        case 2: matrix<T2>(p,r,c).normalize(); break;
        case 3: matrix<T3>(p,r,c).normalize(); break;
    }
    GUARD_END
}

extern "C" RET eigen_random(int code, void* p, int r, int c)
{
    GUARD_START
    switch (code) {
        case 0: matrix<T0>(p,r,c) = MatrixXf::Random(r,c); break;
        case 1: matrix<T1>(p,r,c) = MatrixXd::Random(r,c); break;
        case 2: matrix<T2>(p,r,c) = MatrixXcf::Random(r,c); break;
        case 3: matrix<T3>(p,r,c) = MatrixXcd::Random(r,c); break;
    }
    GUARD_END
}

extern "C" RET eigen_identity(int code, void* p, int r, int c)
{
    GUARD_START
    switch (code) {
        case 0: matrix<T0>(p,r,c) = MatrixXf::Identity(r,c); break;
        case 1: matrix<T1>(p,r,c) = MatrixXd::Identity(r,c); break;
        case 2: matrix<T2>(p,r,c) = MatrixXcf::Identity(r,c); break;
        case 3: matrix<T3>(p,r,c) = MatrixXcd::Identity(r,c); break;
    }
    GUARD_END
}

