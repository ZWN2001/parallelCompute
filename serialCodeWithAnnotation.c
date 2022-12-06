#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>

// Calculate sum of distance while combining different pivots. Complexity : O( n^2 )
//这部分对应的应该就是给的题目简介的计算方法
double SumDistance(const int k, const int n, const int dim, double *coord, const int *pivots) {
    double *rebuiltCoord = (double *) malloc(sizeof(double) * n * k);
    int i;
    for (i = 0; i < n * k; i++) {
        rebuiltCoord[i] = 0;
    }

    // Rebuild coordinates. New coordinate of one point is its distance to each pivot.
    //重建坐标。一个点的新坐标是它到每个支撑点的距离。
    for (i = 0; i < n; i++) {
        int ki;
        for (ki = 0; ki < k; ki++) {
            double distance = 0;
            int pivoti = pivots[ki];
            int j;
            for (j = 0; j < dim; j++) {
                distance += pow(coord[pivoti * dim + j] - coord[i * dim + j], 2);
            }
            rebuiltCoord[i * k + ki] = sqrt(distance);
        }
    }

    // Calculate the sum of Chebyshev distance with rebuilt coordinates between every points
    //用重建的坐标计算各点之间的切比雪夫距离和
    double chebyshevSum = 0;
    for (i = 0; i < n; i++) {
        int j;
        for (j = 0; j < n; j++) {
            double chebyshev;
            int ki;
            for (ki = 0; ki < k; ki++) {
                double dis = fabs(rebuiltCoord[i * k + ki] - rebuiltCoord[j * k + ki]);
                chebyshev = dis > chebyshev ? dis : chebyshev;
            }
            chebyshevSum += chebyshev;
        }
    }

    free(rebuiltCoord);

    return chebyshevSum;
}

/// Recursive function Combination() : combine pivots and calculate the sum of distance while combining different pivots.
/// ki  : current depth of the recursion 当前递归深度
/// k   : number of pivots  支撑点数量
/// n   : number of points  数据点 数量
/// dim : dimension of metric space  度量空间的维数
/// M   : number of combinations to store  要存储的组合数，也就是说我们要找出多少组组合
/// coord  : coordinates of points 点的坐标
/// pivots : indexes of pivots 支撑点的index
/// maxDistanceSum  : the largest M distance sum
/// maxDisSumPivots : the top M pivots combinations
/// minDistanceSum  : the smallest M distance sum
/// minDisSumPivots : the bottom M pivots combinations
void Combination(int ki, const int k, const int n, const int dim, const int M, double *coord, int *pivots,
                 double *maxDistanceSum, int *maxDisSumPivots, double *minDistanceSum, int *minDisSumPivots) {
    //已经找到 k个支撑点，进行计算与排序
    if (ki == k - 1) {
        int i;
        for (i = pivots[ki - 1] + 1; i < n; i++) {
            pivots[ki] = i;

            // Calculate sum of distance while combining different pivots.
            //在不同支撑点的排列组合之间计算距离
            double distanceSum = SumDistance(k, n, dim, coord, pivots);

            // put data at the end of array
            maxDistanceSum[M] = distanceSum;
            minDistanceSum[M] = distanceSum;
            int kj;
            for (kj = 0; kj < k; kj++) {
                maxDisSumPivots[M * k + kj] = pivots[kj];
            }
            for (kj = 0; kj < k; kj++) {
                minDisSumPivots[M * k + kj] = pivots[kj];
            }
            // sort
            int a;
            for (a = M; a > 0; a--) {
                if (maxDistanceSum[a] > maxDistanceSum[a - 1]) {
                    double temp1 = maxDistanceSum[a];
                    maxDistanceSum[a] = maxDistanceSum[a - 1];
                    maxDistanceSum[a - 1] = temp1;
                    int km;
                    for (km = 0; km < k; km++) {
                        int temp = maxDisSumPivots[a * k + km];
                        maxDisSumPivots[a * k + km] = maxDisSumPivots[(a - 1) * k + km];
                        maxDisSumPivots[(a - 1) * k + km] = temp;
                    }
                }
            }
            for (a = M; a > 0; a--) {
                if (minDistanceSum[a] < minDistanceSum[a - 1]) {
                    double temp1 = minDistanceSum[a];
                    minDistanceSum[a] = minDistanceSum[a - 1];
                    minDistanceSum[a - 1] = temp1;
                    int kn;
                    for (kn = 0; kn < k; kn++) {
                        int temp = minDisSumPivots[a * k + kn];
                        minDisSumPivots[a * k + kn] = minDisSumPivots[(a - 1) * k + kn];
                        minDisSumPivots[(a - 1) * k + kn] = temp;
                    }
                }
            }
        }
        return;
    }

    // Recursively call Combination() to combine pivots 递归对支撑点进行排列组合
    int i;
    //传入的pivots是从下标1开始的，所以前一位值为-1，加一后i是从0开始
    //对所有数据点进行遍历
    for (i = pivots[ki - 1] + 1; i < n; i++) {
        pivots[ki] = i;
        Combination(ki + 1, k, n, dim, M, coord, pivots, maxDistanceSum, maxDisSumPivots, minDistanceSum,
                    minDisSumPivots);

        /** Iteration Log : pivots computed, best pivots, max distance sum, min distance sum pivots, min distance sum
        *** You can delete the logging code. **/
        if (ki == k - 2) {
            int kj;
            for (kj = 0; kj < k; kj++) {
                printf("%d ", pivots[kj]);
            }
            putchar('\t');
            for (kj = 0; kj < k; kj++) {
                printf("%d ", maxDisSumPivots[kj]);
            }
            printf("%lf\t", maxDistanceSum[0]);
            for (kj = 0; kj < k; kj++) {
                printf("%d ", minDisSumPivots[kj]);
            }
            printf("%lf\n", minDistanceSum[0]);
        }
    }
}

int main(int argc, char *argv[]) {

    // 数据集文件名
    char *filename = (char *) "uniformvector-2dim-5h.txt";
    //自定义文件路径，可有可无
    if (argc == 2) {
        filename = argv[1];
    } else if (argc != 1) {
        printf("Usage: ./pivot <filename>\n");
        return -1;
    }
    // M : number of combinations to store
    const int M = 1000;
    // dim : dimension of metric space 度量空间的维数
    int dim;
    // n : number of points
    int n;
    // k : number of pivots 支撑点数量
    int k;

    // Read parameter
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("%s file not found.\n", filename);
        return -1;
    }
    fscanf(file, "%d", &dim);
    fscanf(file, "%d", &n);
    fscanf(file, "%d", &k);
    printf("dim = %d, n = %d, k = %d\n", dim, n, k);

    // Start timing
    struct timeval start;

    // Read Data
    //dim列n行double数据
    double *coord = (double *) malloc(sizeof(double) * dim * n);
    int i;
    for (i = 0; i < n; i++) {
        int j;
        for (j = 0; j < dim; j++) {
            fscanf(file, "%lf", &coord[i * dim + j]);
        }
    }
    fclose(file);
    gettimeofday(&start, NULL);

    // maxDistanceSum : the largest M distance sum
    //最大的M个距离值
    double *maxDistanceSum = (double *) malloc(sizeof(double) * (M + 1));
    for (i = 0; i < M; i++) {
        maxDistanceSum[i] = 0;
    }
    // maxDisSumPivots : the top M pivots combinations
    //M行K列支撑点集合
    int *maxDisSumPivots = (int *) malloc(sizeof(int) * k * (M + 1));
    for (i = 0; i < M; i++) {
        int ki;
        for (ki = 0; ki < k; ki++) {
            maxDisSumPivots[i * k + ki] = 0;
        }
    }
    // minDistanceSum : the smallest M distance sum
    double *minDistanceSum = (double *) malloc(sizeof(double) * (M + 1));
    for (i = 0; i < M; i++) {
        minDistanceSum[i] = __DBL_MAX__;
    }
    // minDisSumPivots : the bottom M pivots combinations
    int *minDisSumPivots = (int *) malloc(sizeof(int) * k * (M + 1));
    for (i = 0; i < M; i++) {
        int ki;
        for (ki = 0; ki < k; ki++) {
            minDisSumPivots[i * k + ki] = 0;
        }
    }

    // temp : indexes of pivots with dummy array head
    int *temp = (int *) malloc(sizeof(int) * (k + 1));
    temp[0] = -1;

    // Main loop. Combine different pivots with recursive function and evaluate them. Complexity : O( n^(k+2) )
    //需要优化的目标
    Combination(0, k, n, dim, M, coord, &temp[1], maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);

    // End timing
    struct timeval end;
    gettimeofday(&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0);

    // Store the result
    FILE *out = fopen("result.txt", "w");
    for (i = 0; i < M; i++) {
        int ki;
        for (ki = 0; ki < k - 1; ki++) {
            fprintf(out, "%d ", maxDisSumPivots[i * k + ki]);
        }
        fprintf(out, "%d\n", maxDisSumPivots[i * k + k - 1]);
    }
    for (i = 0; i < M; i++) {
        int ki;
        for (ki = 0; ki < k - 1; ki++) {
            fprintf(out, "%d ", minDisSumPivots[i * k + ki]);
        }
        fprintf(out, "%d\n", minDisSumPivots[i * k + k - 1]);
    }
    fclose(out);

    // Log
    int ki;
    printf("max : ");
    for (ki = 0; ki < k; ki++) {
        printf("%d ", maxDisSumPivots[ki]);
    }
    printf("%lf\n", maxDistanceSum[0]);
    printf("min : ");
    for (ki = 0; ki < k; ki++) {
        printf("%d ", minDisSumPivots[ki]);
    }
    printf("%lf\n", minDistanceSum[0]);
    // for(i=0; i<M; i++){
    // int ki;
    // for(ki=0; ki<k; ki++){
    // printf("%d\t", maxDisSumPivots[i*k + ki]);
    // }
    // printf("%lf\n", maxDistanceSum[i]);
    // }
    // for(i=0; i<M; i++){
    // int ki;
    // for(ki=0; ki<k; ki++){
    // printf("%d\t", minDisSumPivots[i*k + ki]);
    // }
    // printf("%lf\n", minDistanceSum[i]);
    // }

    return 0;
}
