#include "eigen-sparse.h"
#include <Eigen/LU>
#include <Eigen/LeastSquares>

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
    a->prune(T(1));
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
