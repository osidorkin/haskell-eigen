#include <string.h>
#include <stdio.h>
#include <complex>
#include <memory>


#define EIGEN_MPL2_ONLY
#define EIGEN2_SUPPORT
#define EIGEN_NO_EIGEN2_DEPRECATED_WARNING
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
void eigen_assert_fail(const char* condition, const char* function, const char* file, int line);
#define eigen_assert(x) do {\
    if (!(x)) eigen_assert_fail(#x, __PRETTY_FUNCTION__, __FILE__, __LINE__);\
} while(false)

extern "C" {

bool eigen_initParallel();
int eigen_getNbThreads();
void eigen_setNbThreads(int);
extern bool eigen_inited;

} // end extern "C"

#define GUARD_START try { assert(eigen_inited); do {
#define GUARD_END } while(false); return 0; } catch (const std::exception& ex) { return strdup(ex.what()); }

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

typedef float T0;
typedef double T1;
typedef std::complex<float> T2;
typedef std::complex<double> T3;
