#ifdef DEBUG
#pragma GCC optimize(0)
#else
#pragma GCC optimize(3)
#endif

///10�˲��Ի���8�߳�main.cԼ 3500ms��С�������߳������������½�����16,24����32�̺߳�������������3100ms��
///48�߳����32�̻߳���û���������������ܸ����Ի��йأ�
///56��64�̻߳���������
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