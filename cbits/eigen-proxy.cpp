#include "eigen-proxy.h"
#include <Eigen/LU>
#include <Eigen/LeastSquares>
#include <stdio.h>
#include <sstream>

static bool inited = eigen_initParallel();

class eigen_assert_exception : public std::exception {
    std::string _what;
public:
    eigen_assert_exception(const std::string& what) : _what(what) {}
    ~eigen_assert_exception() throw() {}
    const char* what() const throw () { return _what.c_str(); }
};

void eigen_assert_fail(const char* condition, const char* function, const char* file, int line) {
    std::ostringstream os;
    os << "assertion failed: " << condition << " in function " << function << " at " << file << ":" << line << std::endl;
    throw eigen_assert_exception(os.str());
}

typedef float T0;
typedef double T1;
typedef std::complex<float> T2;
typedef std::complex<double> T3;

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

extern "C" bool eigen_initParallel() {
    initParallel();
    return true;
}

extern "C" void eigen_setNbThreads(int n) {
    setNbThreads(n);
}

extern "C" int eigen_getNbThreads() {
    return nbThreads();
}

template <class T>
RET sparse_fromList(int rows, int cols, void* data, int size, void** pr) {
    typedef SparseMatrix<T> M;
    typedef Triplet<T> E;
    std::auto_ptr<M> a(new M(rows, cols));
    a->setFromTriplets((E*)data, (E*)data + size);
    *(M**)pr = a.release();
    return 0;
}
API(sparse_fromList, (int code, int rows, int cols, void* data, int size, void** pr), (rows,cols,data,size,pr));

template <class  T>
RET sparse_toList(void* p, void* q, int s) {
    int n = 0;
    typedef SparseMatrix<T> M;
    typedef Triplet<T> E;
    M* m = (M*)p;
    for (int k = 0; k < m->outerSize(); ++k) {
        for (typename M::InnerIterator i(*m, k); i; ++i) {
            if (n >= s)
                return strdup("sparse_toList: buffer overrun detected");
            ((E*)q)[n++] = E(i.row(), i.col(), i.value());
        }
    }
    return n < s ? strdup("sparse_toList: buffer underrun detected") : 0;
}
API(sparse_toList, (int code, void* p, void* q, int s), (p,q,s));

template <class T>
RET sparse_free(void* p) {
    delete (SparseMatrix<T>*)p;
    return 0;
}
API(sparse_free, (int code, void* p), (p));

#define SPARSE_UNOP_INPLACE(name)\
template <class T>\
RET sparse_##name(void* p, void** pr) {\
    typedef SparseMatrix<T> M;\
    std::auto_ptr<M> a(new M(*(M*)p));\
    a->name();\
    *(M**)pr = a.release();\
    return 0;\
}\
API(sparse_##name, (int code, void* p, void** pr), (p, pr));

SPARSE_UNOP_INPLACE(makeCompressed);
SPARSE_UNOP_INPLACE(uncompress);

#define SPARSE_UNOP(name)\
template <class T>\
RET sparse_##name(void* p, void** pr) {\
    typedef SparseMatrix<T> M;\
    *(M**)pr = new M(((M*)p)->name());\
    return 0;\
}\
API(sparse_##name, (int code, void* p, void** pr), (p, pr));

SPARSE_UNOP(adjoint);
SPARSE_UNOP(transpose);

template <class T>
RET sparse_pruned(void* p, void** pr) {
    typedef SparseMatrix<T> M;
    std::auto_ptr<M> a(new M(*(M*)p));
    a->prune(T());
    *(M**)pr = a.release();
    return 0;
}
API(sparse_pruned, (int code, void* p, void** pr), (p, pr));

template <class T>
RET sparse_prunedRef(void* p, void* q, void** pr) {
    typedef SparseMatrix<T> M;
    std::auto_ptr<M> a(new M(*(M*)p));
    a->prune(*(T*)q);
    *(M**)pr = a.release();
    return 0;
}
API(sparse_prunedRef, (int code, void* p, void* q, void** pr), (p, q, pr));

template <class T>
RET sparse_scale(void* p, void* q, void** pr) {
    typedef SparseMatrix<T> M;
    *(M**)pr = new M(*(T*)q * *(M*)p);
    return 0;
}
API(sparse_scale, (int code, void* p, void* q, void** pr), (p, q, pr));

template <class T>
RET sparse_coeff(void* p, int row, int col, void* pr) {
    *(T*)pr = ((SparseMatrix<T>*)p)->coeff(row, col);
    return 0;
}
API(sparse_coeff, (int code, void* p, int row, int col, void* pr), (p, row, col, pr));

#define SPARSE_PROP(name,type)\
template <class T>\
RET sparse_##name(void* p, void* pr) {\
    *(type*)pr = ((SparseMatrix<T>*)p)->name();\
    return 0;\
}\
API(sparse_##name, (int code, void* p, void* pr), (p, pr));

SPARSE_PROP(cols, int);
SPARSE_PROP(rows, int);
SPARSE_PROP(innerSize, int);
SPARSE_PROP(outerSize, int);
SPARSE_PROP(nonZeros, int);
SPARSE_PROP(isCompressed, int);
SPARSE_PROP(norm, T);
SPARSE_PROP(squaredNorm, T);
SPARSE_PROP(blueNorm, T);

#define SPARSE_BINOP(name,op)\
template <class T>\
RET sparse_##name(void* p, void* q, void** pr) {\
    typedef SparseMatrix<T> M;\
    *(M**)pr = new M(*(M*)p op *(M*)q);\
    return 0;\
}\
API(sparse_##name, (int code, void* p, void* q, void** pr), (p, q, pr));

SPARSE_BINOP(add, +);
SPARSE_BINOP(sub, -);
SPARSE_BINOP(mul, *);

template <class T>
RET sparse_block(void* p, int row, int col, int rows, int cols, void** pr) {
    typedef SparseMatrix<T> M;
    *(M**)pr = new M(((M*)p)->block(row,col,rows,cols));
    return 0;
}
API(sparse_block, (int code, void* p, int row, int col, int rows, int cols, void** pr), (p, row, col, rows, cols, pr));

template <class T> bool isZero(T x) { return x == 0; }
template <class T> bool isZero(std::complex<T> x) { return x.real() == 0 && x.imag() == 0; }

template <class T>
RET sparse_fromMatrix(void* p, int rows, int cols, void** pq) {
    typedef SparseMatrix<T> M;
    typedef Map< Matrix<T,Dynamic,Dynamic> > MapMatrix;
    MapMatrix src((const T*)p, rows, cols);
    std::auto_ptr<M> dst(new M(rows, cols));
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            T val = src.coeff(row,col);
            if (!isZero(val))
                dst->insert(row, col) = val;
        }
    }
    *(M**)pq = dst.release();
    return 0;
}
API(sparse_fromMatrix, (int code, void* p, int rows, int cols, void** pq), (p,rows,cols,pq));

template <class T>
RET sparse_toMatrix(void* p, void* q, int rows, int cols) {
    typedef Map< Matrix<T,Dynamic,Dynamic> > MapMatrix;
    MapMatrix((T*)q, rows, cols) = *(SparseMatrix<T>*)p;
    return 0;
}
API(sparse_toMatrix, (int code, void* p, void* q, int rows, int cols), (p,q,rows,cols));
