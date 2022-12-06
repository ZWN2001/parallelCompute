#include "main.h"

double SumDistance(const int k, const int n, const int dim, double* __restrict coord, double* __restrict pivots) {
    double *rebuiltCoord = (double *) malloc(sizeof(double) * n * k);
    memset(rebuiltCoord, 0, sizeof(double) * n * k);
int i,ki,j;
double dis,chebyshev,distance;
#pragma omp parallel for shared(pivots,coord,rebuiltCoord) private(ki,distance,j)
    for ( i = 0; i < n; i++) {
        for ( ki = 0; ki < k; ki++) {
            distance = 0;
            for ( j = 0; j < dim; j++) {
                distance += square(pivots[ki * dim + j] - coord[i * dim + j]);
            }
            rebuiltCoord[i * k + ki] = sqrt(distance);
        }
    }

    double chebyshevSum = 0;
    for ( i = 0; i < n; i++) {
        for ( j = i + 1; j < n; j++) {
             chebyshev = 0;
            for ( ki = 0; ki < k; ki++) {
                 dis = fabs(rebuiltCoord[i * k + ki] - rebuiltCoord[j * k + ki]);
                if(dis > chebyshev){
                    chebyshev = dis;
                }
            }
            chebyshevSum += chebyshev;
        }
    }

    free(rebuiltCoord);
    return chebyshevSum;
}

int main(int argc, char *argv[]) {
    omp_set_num_threads(THREAD);
    // filename : input file namespace
    char *filename = (char *) "uniformvector-2dim-5h.txt";
    if (argc == 2) {
        filename = argv[1];
    } else if (argc != 1) {
        printf("Usage: ./pivot <filename>\n");
        return -1;
    }
    // M : number of combinationKind to store
    int M = 1000;
    // dim : dimension of metric space
    int dim;
    // n : number of points
    int n;
    // k : number of pivots
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
    double *coord = (double *) malloc(sizeof(double) * n * dim );
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < dim; j++) {
            fscanf(file, "%lf", &coord[i * dim + j]);
        }
    }
    fclose(file);
    gettimeofday(&start, NULL);

//========================================================================
    int combinationSize = combination(n,k);
    double *combinationKind = (double *) malloc(sizeof(double) * combinationSize * k * dim);
    int *combinationsIndex = (int *) malloc(sizeof(int) * combinationSize * k);

//·························初始化可能的组合情况················
    int* temp = (int *)malloc(sizeof(int) * (k + 1));
    int base = 0;
    int index_base = 0;
    int w,offset;
    for (int i = 1; i <= k; ++i) {
        temp[i - 1] = i;
    }
    temp[k] = n + 1;
    int u = 0,s;
    while (u < k) {
        offset = 0;
        for (int i = 0; i < k; i++) {
            s = temp[i] - 1;
            for (w = 0; w < dim; ++w) {
                combinationKind[base + offset] = coord[s * dim + w];
                offset++;
            }
            combinationsIndex[index_base + i] = s;
        }
        base += k * dim;
        index_base += k;

        u = 0;
        while (u < k && temp[u] + 1 == temp[u + 1]) {
            temp[u] = u + 1;
            ++u;
        }
        ++temp[u];
    }
//·························初始化可能的组合情况················

    double tempMaxVal[THREAD][M + 1];
    double tempMinVal[THREAD][M + 1];
    int indexMax[THREAD][M][k];
    int indexMin[THREAD][M][k];

    memset(tempMaxVal, 0, sizeof(double) * THREAD * (M + 1));
    memset(tempMinVal, MAX_CHAR_DOUBLE, sizeof(double) * THREAD * (M + 1));
    for (int i = 0; i < THREAD; ++i) {
        tempMaxVal[i][M] = __DBL_MAX__;
        tempMinVal[i][M] = -1;
    }


#pragma omp parallel for default(none) shared(M, combinationSize, k, n, dim, coord, combinationKind, tempMaxVal, tempMinVal, combinationsIndex, indexMax, indexMin)
    for (int i = 0; i < combinationSize; ++i) {
        const int thread_id = omp_in_parallel() ? omp_get_thread_num() : 0;
        double distance = SumDistance(k, n, dim, coord, combinationKind + i * k * dim);
        int j = 0;
        //下面的两个while比较吃性能,但是相对占比不高
        while (distance > tempMaxVal[thread_id][j]) {
            j++;
        }
        if (j != 0) {
            memcpy(tempMaxVal[thread_id], tempMaxVal[thread_id] + 1, sizeof(double) * (j - 1));
            tempMaxVal[thread_id][j - 1] = distance;
            memcpy(indexMin[thread_id][0], indexMin[thread_id][0] + k, sizeof(int) * (j - 1) * k);
            for (int l = 0; l < k; ++l) {
                indexMin[thread_id][j - 1][l] = combinationsIndex[i * k + l];
            }
        }

        j = 0;
        //从大到小的排列
        while (distance < tempMinVal[thread_id][j]) {
            j++;
        }
        if (j != 0) {
            memcpy(tempMinVal[thread_id], tempMinVal[thread_id] + 1, sizeof(double) * (j - 1));
            tempMinVal[thread_id][j - 1] = distance;
            memcpy(indexMax[thread_id][0], indexMax[thread_id][0] + k, sizeof(int) * (j - 1) * k);
            for (int l = 0; l < k; ++l) {
                indexMax[thread_id][j - 1][l] = combinationsIndex[i * k + l];
            }
        }
    }

    int minCount = 0;
    int maxCount = 0;
    int threadIndex[THREAD];
    int finalMin[M][k];
    int finalMax[M][k];

    double maxDistance;
    int maxThreadIndex;
    for (int i = 0; i < THREAD; ++i) {
        threadIndex[i] = M - 1;
    }
    while (maxCount < 1000) {
        maxDistance = -1;
        maxThreadIndex = -1;
        for (int i = 0; i < THREAD; ++i) {
            if (tempMaxVal[i][threadIndex[i]] > maxDistance) {
                maxDistance = tempMaxVal[i][threadIndex[i]];
                maxThreadIndex = i;
            }
        }
        for (int i = 0; i < k; ++i) {
            finalMax[maxCount][i] = indexMin[maxThreadIndex][threadIndex[maxThreadIndex]][i];
        }
        threadIndex[maxThreadIndex]--;
        maxCount++;
    }

    double minDistance;
    int minThreadIndex;
    for (int i = 0; i < THREAD; ++i) {
        threadIndex[i] = M - 1;
    }
    while (minCount < 1000) {
        minDistance = __DBL_MAX__;
        minThreadIndex = -1;
        for (int i = 0; i < THREAD; ++i) {
            if (tempMinVal[i][threadIndex[i]] < minDistance) {
                minDistance = tempMinVal[i][threadIndex[i]];
                minThreadIndex = i;
            }
        }
        for (int i = 0; i < k; ++i) {
            finalMax[maxCount][i] = indexMin[minThreadIndex][threadIndex[minThreadIndex]][i];
        }
        threadIndex[minThreadIndex]--;
        minCount++;
    }
    //==========================================================

    // End timing
    struct timeval end;
    gettimeofday(&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_usec - start.tv_usec) / 1000.0);

    // Store the result
    FILE *out = fopen("result.txt", "w");
    for (int i = 0; i < M; i++) {
        for (int ki = 0; ki < k - 1; ki++) {
            fprintf(out, "%d ", finalMax[i][ki]);
        }
        fprintf(out, "%d\n", finalMax[i][k - 1]);
    }
    for (int i = 0; i < M; i++) {
        for (int ki = 0; ki < k - 1; ki++) {
            fprintf(out, "%d ", finalMin[i][ki]);
        }
        fprintf(out, "%d\n", finalMin[i][k - 1]);
    }
    fclose(out);

    return 0;
}

int combination(int n, int k) {
    unsigned long long res = 1;
    for (int i = n - k; i < n; i++) {
        res *= i + 1;
    }
    for (int i = 0; i < k; ++i) {
        res /= i + 1;
    }
    return res;
}
