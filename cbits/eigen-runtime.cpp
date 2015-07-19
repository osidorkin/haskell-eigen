#include "eigen-runtime.h"
#include <sstream>
#include <Eigen/Core>

using namespace Eigen;

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

extern "C" {

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

bool eigen_inited = eigen_initParallel();

} // extern "C"
