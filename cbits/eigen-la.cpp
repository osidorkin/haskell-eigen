#include "eigen-la.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/LeastSquares>

using namespace Eigen;

template <class T>
RET rank(Decomposition d, int* v, const void* p, int r, int c) {
    typedef Map< Matrix<T,Dynamic,Dynamic> > MapMatrix;
    MapMatrix A((const T*)p,r,c);
    switch (d) {
        case ::FullPivLU:
            *v = A.fullPivLu().rank();
            break;
        case ::ColPivHouseholderQR:
            *v = A.colPivHouseholderQr().rank();
            break;
        case ::FullPivHouseholderQR:
            *v = A.fullPivHouseholderQr().rank();
            break;
        case ::JacobiSVD:
            *v = A.jacobiSvd(ComputeThinU | ComputeThinV).rank();
            break;
        default:
            return strdup("Selected decomposition doesn't support rank revealing.");
    }
    return 0;
}
API(rank, (int code, Decomposition d, int* v, const void* p, int r, int c), (d,v,p,r,c));

template <class T>
RET kernel(Decomposition d, void** p0, int* r0, int* c0, const void* p1, int r1, int c1) {
    typedef Map< Matrix<T,Dynamic,Dynamic> > MapMatrix;
    if (d != ::FullPivLU)
        return strdup("Selected decomposition doesn't support kernel revealing.");
    MapMatrix A((const T*)p1,r1,c1);
    Matrix<T,Dynamic,Dynamic> B = A.fullPivLu().kernel();
    *r0 = B.rows();
    *c0 = B.cols();
    *p0 = malloc(*r0 * *c0 * sizeof(T));
    MapMatrix((T*)*p0, *r0, *c0) = B;
    return 0;
}
API(kernel, (int code, Decomposition d, void** p0, int* r0, int* c0, const void* p1, int r1, int c1), (d,p0,r0,c0,p1,r1,c1));

template <class T>
RET image(Decomposition d, void** p0, int* r0, int* c0, const void* p1, int r1, int c1) {
    typedef Map< Matrix<T,Dynamic,Dynamic> > MapMatrix;
    if (d != ::FullPivLU)
        return strdup("Selected decomposition doesn't support image revealing.");
    MapMatrix A((const T*)p1,r1,c1);
    Matrix<T,Dynamic,Dynamic> B = A.fullPivLu().image(A);
    *r0 = B.rows();
    *c0 = B.cols();
    *p0 = malloc(*r0 * *c0 * sizeof(T));
    MapMatrix((T*)*p0, *r0, *c0) = B;
    return 0;
}
API(image, (int code, Decomposition d, void** p0, int* r0, int* c0, const void* p1, int r1, int c1), (d,p0,r0,c0,p1,r1,c1));

template <class T>
RET solve(Decomposition d,
    void* px, int rx, int cx,
    const void* pa, int ra, int ca,
    const void* pb, int rb, int cb)
{
    typedef Map< Matrix<T,Dynamic,Dynamic> > MapMatrix;
    MapMatrix x((T*)px, rx, cx);
    MapMatrix A((const T*)pa, ra, ca);
    MapMatrix b((const T*)pb, rb, cb);
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
    return 0;
}
API(solve, (int code, Decomposition d,
    void* px, int rx, int cx,
    const void* pa, int ra, int ca,
    const void* pb, int rb, int cb), (d,px,rx,cx,pa,ra,ca,pb,rb,cb));

template <class T>
RET relativeError(void* e,
    const void* px, int rx, int cx,
    const void* pa, int ra, int ca,
    const void* pb, int rb, int cb)
{
    typedef Map< Matrix<T,Dynamic,Dynamic> > MapMatrix;
    MapMatrix x((const T*)px, rx, cx);
    MapMatrix A((const T*)pa, ra, ca);
    MapMatrix b((const T*)pb, rb, cb);
    *(T*)e = (A*x - b).norm() / b.norm();
    return 0;
}
API(relativeError, (int code, void* e,
    const void* px, int rx, int cx,
    const void* pa, int ra, int ca,
    const void* pb, int rb, int cb), (e,px,rx,cx,pa,ra,ca,pb,rb,cb));

