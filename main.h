#ifdef DEBUG
#pragma GCC optimize(0)
#else
#pragma GCC optimize(3)
#endif

///10核测试机，8线程main.c约 3500ms，小幅增加线程数出现性能下降（如16,24），32线程后有性能提升到3100ms，
///48线程相较32线程基本没有性能提升，可能跟测试机有关，
///56、64线程基本无提升
#define THREAD 64

#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<string.h>
#include <omp.h>

#define MAX_CHAR_DOUBLE 0x7f
//#define INV_SIGN_BIT 0x7fffffff
//#define MAX_CHAR_INT 0x3f

#define square(a) ((a)*(a))

inline double SumDistance(int k, int n, int dim, double* __restrict coord, double* __restrict pivots);
int combination(int n, int k);