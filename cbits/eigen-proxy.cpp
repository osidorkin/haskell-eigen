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

typedef Map< Matrix<double,Dynamic,Dynamic> > MapMatrix;

extern "C" {

const char* eigen_add(
    double* data, int rows, int cols, 
    const double* data1, int rows1, int cols1,
    const double* data2, int rows2, int cols2)
{
    GUARD_START
    MapMatrix(data, rows, cols) = MapMatrix(data1, rows1, cols1) + MapMatrix(data2, rows2, cols2);;
    GUARD_END
}

const char* eigen_sub(
    double* data, int rows, int cols,
    const double* data1, int rows1, int cols1,
    const double* data2, int rows2, int cols2)
{
    GUARD_START
    MapMatrix(data, rows, cols) = MapMatrix(data1, rows1, cols1) - MapMatrix(data2, rows2, cols2);;
    GUARD_END
}

const char* eigen_mul(
    double* data, int rows, int cols,
    const double* data1, int rows1, int cols1,
    const double* data2, int rows2, int cols2)
{
    GUARD_START
    MapMatrix(data, rows, cols) = MapMatrix(data1, rows1, cols1) * MapMatrix(data2, rows2, cols2);;
    GUARD_END
}

double eigen_norm(const double* data, int rows, int cols) { return MapMatrix(data, rows, cols).norm();  }
double eigen_squaredNorm(const double* data, int rows, int cols) { return MapMatrix(data, rows, cols).squaredNorm(); }
double eigen_blueNorm(const double* data, int rows, int cols) { return MapMatrix(data, rows, cols).squaredNorm(); }
double eigen_hypotNorm(const double* data, int rows, int cols) { return MapMatrix(data, rows, cols).hypotNorm(); }
double eigen_sum(const double* data, int rows, int cols) { return MapMatrix(data, rows, cols).sum(); }
double eigen_prod(const double* data, int rows, int cols) { return MapMatrix(data, rows, cols).prod(); }
double eigen_mean(const double* data, int rows, int cols) { return MapMatrix(data, rows, cols).mean(); }
double eigen_trace(const double* data, int rows, int cols) { return MapMatrix(data, rows, cols).trace(); }
double eigen_determinant(const double* data, int rows, int cols) { return MapMatrix(data, rows, cols).determinant(); }

const char* eigen_inverse(double* data, int rows, int cols, const double* data1, int rows1, int cols1)
{
    GUARD_START
    puts("inverse");
    MapMatrix(data, rows, cols) = MapMatrix(data1, rows1, cols1).inverse();
    GUARD_END
}

const char* eigen_adjoint(double* data, int rows, int cols, const double* data1, int rows1, int cols1)
{
    GUARD_START
    puts("adjoint");
    MapMatrix(data, rows, cols) = MapMatrix(data1, rows1, cols1).adjoint();
    GUARD_END
}

const char* eigen_conjugate(double* data, int rows, int cols, const double* data1, int rows1, int cols1)
{
    GUARD_START
    puts("conjugate");
    MapMatrix(data, rows, cols) = MapMatrix(data1, rows1, cols1).conjugate();
    GUARD_END
}

const char* eigen_transpose(double* data, int rows, int cols, const double* data1, int rows1, int cols1)
{
    GUARD_START
    puts("transpose");
    MapMatrix(data, rows, cols) = MapMatrix(data1, rows1, cols1).transpose();
    GUARD_END
}

const char* eigen_normalize(double* data, int rows, int cols)
{
    GUARD_START
    puts("normalize");
    MapMatrix(data, rows, cols).normalize();
    GUARD_END
}


const char* eigen_solve(Decomposition d,
    double* px, int rx, int cx, // x
    const double* pa, int ra, int ca, // A
    const double* pb, int rb, int cb) // b
{
    GUARD_START
    MapMatrix x(px, rx, cx);
    MapMatrix A(pa, ra, ca);
    MapMatrix b(pb, rb, cb);
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
    GUARD_END
}


const char* eigen_relativeError(double& e,
    const double* px, int rx, int cx, // x
    const double* pa, int ra, int ca, // A
    const double* pb, int rb, int cb) // b
{
    GUARD_START
    MapMatrix x(px, rx, cx);
    MapMatrix A(pa, ra, ca);
    MapMatrix b(pb, rb, cb);
    e = (A*x - b).norm() / b.norm();
    GUARD_END
}


bool eigen_initParallel() {
    initParallel();
    return true;
}

void eigen_setNbThreads(int n) {
    setNbThreads(n);
}

int eigen_getNbThreads() {
    return nbThreads();
}

}
