# 多核平台上的并行计算

总的来说个人理解就是尽可能均衡地平均负载，压榨每个核的计算能力，减少非必要的非计算耗时。

## 实验题目

<img src="多核平台上的并行计算\3.png" style="zoom:70%;" />

其中，k是给定值（或者说，自定义值）



**切比雪夫距离**（**Chebyshev distance**）或是 $L_∞$ 度量是向量空间中的一种度量，二个点之间的距离定义为其各座标数值差的最大值。以$(x_1,y_1)$ 和 $(x_2,y_2)$ 二点为例，其切比雪夫距离为 $max(|x_2-x_1|,|y_2-y_1|)$。

若将[国际象棋](https://zh.wikipedia.org/wiki/國際象棋)棋盘放在二维直角座标系中，格子的边长定义为1，坐标的x轴及y轴和棋盘方格平行，原点恰落在某一格的中心点，则**王从一个位置走到其他位置需要的步数恰为两个位置的切比雪夫距离**，因此切比雪夫距离也称为**棋盘距离**。例如位置F6和位置E2的切比雪夫距离为4。任何一个不在棋盘边缘的位置，和周围八个位置的切比雪夫距离都是1。

<img src="多核平台上的并行计算\4.png" style="zoom:100%;" />

在平面几何中，若二点p及q的直角坐标系坐标为 ${(x_{1},y_{1})}$ 及 ${(x_{2},y_{2})}$，则切比雪夫距离为 $D_{Chess}=max(|x_2-x_1|,|y_2-y_1|)$

依以上的度量，以任一点为准，和此点切比雪夫距离为r的点会形成一个正方形，其边长为2r，且各边都和坐标轴平行。 

> k=2时切比雪夫距离等价于曼哈顿距离，使用曼哈顿距离进行计算会更快但并不具有普适性。



# 优化思想与代码分析

原有代码的思路：对所有可能的组合情况，通过递归调用`SumDistance()`函数，找出 k 个支撑点，计算数据点基于支撑点的切比雪夫距离，然后进行排序，输出最大和最小的M个数据。

主要的任务就是对递归算法改成循环，一种方法是将递归提前，也就是说，先提前用递归算法计算出所有的组合情况，再针对所有组合情况进行计算，优点是逻辑简单，缺点是只能串行执行，另一种方法是针对特定的 k 直接写出循环来计算出所有排列组合情况，优点是方便并行化且逻辑简单，缺点是需要针对特定的 k 写特定的代码，导致冗余代码较多，最后一种就是针对不同的 k 通过回溯+剪枝生成不同的排列组合，优点是普适性好，代码简洁，缺点是代码逻辑性强，不易理解，不易并行化。

不过由于代码的热点并不集中在这个函数，所以三种方法都可采取，我最终采取最后一种方法。

```c
 int combinationSize = combination(n,k);
    double *combinationKind = (double *) malloc(sizeof(double) * combinationSize * k * dim);
    int *combinationsIndex = (int *) malloc(sizeof(int) * combinationSize * k);

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
```

在获取了所有的排列组合情况后，针对所有的排列组合情况，通过for循环，计算`distance`，通过插入排序对距离计算结果进行临时存储，在对所有情况遍历完成后，针对所有临时数据进行一次处理即可找出符合条件的数据。我们只需要在for循环前进行并行化（如我们使用的OpenMP）即可。

这里有一个优化细节，对于临时数据的存储与写入，如果不对线程id进行区分，那么在写距离数据时就需要进行加减锁操作，我们使用特定线程写特定数据块来避免锁操作，即：

```c
double tempMaxVal[THREAD][M + 1];
double tempMinVal[THREAD][M + 1];
int indexMax[THREAD][M][k];
int indexMin[THREAD][M][k];
```

这样我们避免了线程间对于写数据的竞争关系。

基于此我们写出了第一版并行的代码。

&nbsp;

## 一些细节的优化

### `pointer aliasing`

我们在方法调用中大量使用指针，但其实编译器并不知道这些指针是否会指向同一内存因而优化会相对保守，然而我们可以确定的是我们方法调用中的指针仅用来读取数据而不会进行写操作，所以我们可以使用`__restrict`关键字让编译器放心地进行优化，如下：

```c
double SumDistance(const int k, const int n, const int dim, double* __restrict coord, double* __restrict pivots)
```

用restrict修饰一个指针，意思就是“只要这个指针活着，我保证这个指针独享这片内存，没有‘别人’可以修改这个指针指向的这片内存，所有修改都得通过这个指针来”。由于这个指针的生命周期是已知的，编译器可以放心大胆地把这片内存中前若干字节用寄存器cache起来。

以个人经验而言，我们编写代码时通常会忽略 pointer aliasing 的问题。在性能剖测时，我通过反编译看到很多冗余的读取指令，才会想到加入 restrict 关键字来提升性能。

### 函数`inline `

我们发现`SumDistance`被高频率调用，在调用栈中我们会看到相当多的`call`指令，所以我们可以直接将函数进行内联：

```c
inline double SumDistance(int k, int n, int dim, double* __restrict coord, double* __restrict pivots);
```

### `pow()`or`square`

优化前，我们如下计算`distance`：

```c
distance += pow(coord[pivoti*dim + j] - coord[i*dim + j] ,2);
```

然而`pow()`函数为了达到普适性其实牺牲了一部分性能，我们完全可以使用乘法进行代替：

```c
//main.h
#define square(a) ((a)*(a))

//main.c
distance += square(pivots[ki * dim + j] - coord[i * dim + j]);
```

### 直接使用`pivots`

在原有代码中我们使用的是在`SumDistance` 方法内通过数组下标确定`pivoti`从而在`coord`中找到对应的支撑点数据，然而上面我们缓存了所有可能的排列组合种类，所以我们可以直接通过传递我们缓存好的数据点信息作为函数支撑点参数，从而避免多次数据定位带来的性能损失。

也就是说，在方法调用时，我们传入对应的kind，即我们对支撑点的枚举：

```c
double distance = SumDistance(k, n, dim, coord, combinationKind + i * k * dim);
```

方法内，我们直接使用支撑点数据计算距离：

```c
for ( i = 0; i < n; i++) {
    for ( ki = 0; ki < k; ki++) {
        distance = 0;
        for ( j = 0; j < dim; j++) {
            distance += square(pivots[ki * dim + j] - coord[i * dim + j]);
        }
        rebuiltCoord[i * k + ki] = sqrt(distance);
    }
}
```

&nbsp;

## 不算优化的优化：线程数的取舍

要使用多少线程才能达到最佳的效率呢？

关于如何计算并发线程数，有两种说法。

第一种，《Java Concurrency in Practice》即《java并发编程实践》8.2节 170页

> 对于计算密集型的任务，一个有*Ncpu*个处理器的系统通常通过使用一个*Ncpu* + 1个线程的线程池来获得最优的利用率（计算密集型的线程恰好在某时因为发生一个页错误或者因其他原因而暂停，刚好有一个“额外”的线程，可以确保在这种情况下CPU周期不会中断工作）。对于包含了 I/O和其他阻塞操作的任务，不是所有的线程都会在所有的时间被调度，因此你需要一个更大的池。为了正确地设置线程池的长度，你必须估算出任务花在等待的时间与用来计算的时间的比率；这个估算值不必十分精确，而且可以通过一些监控工具获得。你还可以选择另一种方法来调节线程池的大小，在一个基准负载下，使用 几种不同大小的线程池运行你的应用程序，并观察CPU利用率的水平。
>
> 给定下列定义：
>
> - *Ncpu* = *CPU*的数量
>
> - *Ucpu* = 目标*CPU*的使用率， 0 <= *Ucpu* <= 1
>
> - *W/C* = 等待时间与计算时间的比率
>
> 为保持处理器达到期望的使用率，最优的池的大小等于：
>
>> 　　*Nthreads* = *Ncpu* x *Ucpu* x (1 + *W/C*)
>
> java环境下你可以使用Runtime来获得CPU的数目：
>
> ```
> int N_CPUS = Runtime.getRuntime().availableProcessors();
> ```
>
> 当然，CPU周期并不是唯一你可以使用线程池管理的资源。其他可以约束资源池大小的资源包括：内存、文件句柄、套接字句柄和数据库连接等。计算这些类型资源池的大小约束非常简单：首先累加出每一个任务需要的这些资源的总量，然后除以可用的总量。所得的结果是池大小的上限。
>
> 当任务需要使用池化的资源时，比如数据库连接，那么线程池的长度和资源池的长度会相互影响。如果每一个任务都需要一个数据库连接，那么连接池的大小就限制了线程池的有效大小；类似地，当线程池中的任务是连接池的唯一消费者时，那么线程池的大小反而又会限制了连接池的有效大小。

第二种，《Programming Concurrency on the JVM Mastering》即《Java 虚拟机并发编程》2.1节 12页

> 为了解决上述难题，我们希望至少可以创建处理器核心数那么多个线程。这就保证了有尽可能多地处理器核心可以投入到解决问题的工作中去。通过代码，我们可以很容易地获取到系统可用的处理器核心数。
>
> 所以，应用程序的最小线程数应该等于可用的处理器核数。如果所有的任务都是计算密集型的，则创建处理器可用核心数那么多个线程就可以了。在这种情况下，创建更多的线程对程序性能而言反而是不利的。因为当有多个仟务处于就绪状态时，处理器核心需要在线程间频繁进行上下文切换，而这种切换对程序性能损耗较大。但如果任务都是IO密集型的，那么我们就需要开更多的线程来提高性能。
>
> 当一个任务执行IO操作时，其线程将被阻塞，于是处理器可以立即进行上下文切换以便处理其他就绪线程。如果我们只有处理器可用核心数那么多个线程的话，则即使有待执行的任务也无法处理，因为我们已经拿不出更多的线程供处理器调度了。
>
> 如果任务有50%的时间处于阻塞状态，则程序所需线程数为处理器可用核心数的两倍。 如果任务被阻塞的时间少于50%，即这些任务是计算密集型的，则程序所需线程数将随之减少，但最少也不应低于处理器的核心数。如果任务被阻塞的时间大于执行时间，即该任务是IO密集型的，我们就需要创建比处理器核心数大几倍数量的线程。
>
> 我们可以计算出程序所需线程的总数，总结如下：
>
> > 线程数 = CPU可用核心数/(1 - 阻塞系数），其中阻塞系数的取值在0和1之间。
> 
>  **计算密集型任务的阻塞系数为0，而IO密集型任务的阻塞系数则接近1**。一个完全阻塞的任务是注定要挂掉的，所以我们无须担心阻塞系数会达到1。
> 
>  为了更好地确定程序所需线程数，我们需要知道下面两个关键参数：
> 
>  - 处理器可用核心数；
>  - 任务的阻塞系数；
> 
>  第一个参数很容易确定，我们甚至可以用之前的方法在运行时查到这个值。但确定阻塞系数就稍微困难一些。我们可以先试着猜测，抑或采用一些性能分析工具或java.lang.management API来确定线程花在系统IO操作上的时间与CPU密集任务所耗时间的比值。

说法一和说法二其实是一个公式。

**至此结论为：**

- IO密集型 = 2Ncpu（可以测试后自己控制大小，2Ncpu一般没问题）（常出现于线程中：数据库数据交互、文件上传下载、网络数据传输等等）

- 计算密集型 = Ncpu（常出现于线程中：复杂算法）

当然说法一中还有一种说法：

> 对于计算密集型的任务，一个有*Ncpu*个处理器的系统通常通过使用一个*Ncpu* + 1个线程的线程池来获得最优的利用率（计算密集型的线程恰好在某时因为发生一个页错误或者因其他原因而暂停，刚好有一个“额外”的线程，可以确保在这种情况下CPU周期不会中断工作）。

即，计算密集型 = Ncpu + 1，但是这种做法导致的多一个CPU上下文切换是否值得，这里不考虑。可自己考量。



&nbsp;

## 优化后的代码

```c
///main,c
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
```

```c
///main.h

#ifdef DEBUG
#pragma GCC optimize(0)
#else
#pragma GCC optimize(3)
#endif

///10核测试机，8线程main.c约 3500ms，小幅增加线程数出现性能下降（如16,24），32线程后有性能提升到3100ms，
///48线程相较32线程基本没有性能提升，可能跟测试机有关，
///56、64线程基本无提升
#define THREAD 48

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
```

10核48线程运行结果如下：

<img src="多核平台上的并行计算\45.png" style="zoom:50%;" />

&nbsp;

## 性能探索

但其实我们会发现，代码的执行时间大多数都不是被以上代码所占用，使用VTune或perf方便我们进一步分析：

我在`SumDistance`方法内额外加入了一行`#pragma omp parallel for`，使得这部分代码在调用栈中被独立地进行显示，如下：

<img src="多核平台上的并行计算\42.png" style="zoom:50%;" />

然后使用perf查看了函数执行耗时的占比，如下：

<img src="多核平台上的并行计算\41.png" style="zoom:50%;" />

显然，整个程序的热点其实都集中在这部分代码上。

由于我给`SumDistance`加了`inline`所以该函数并没有出现在栈中，我们取消内联删除上面的并行化，重新编译运行，查看`SumDistance`内的程序热点：

<img src="多核平台上的并行计算\43.png" style="zoom:50%;" />

<img src="多核平台上的并行计算\44.png" style="zoom:50%;" />

可以发现热点是集中在这一部分的。

但对于这部分代码，我们除手动完成向量化之外暂时没有其他优化途径，而且手动向量化在这段代码中也不易应用，所以这是此次实验中的一个不足之处。

&nbsp;

# 一些优化相关技术知识总结

## 流水线

### 概念

<img src="多核平台上的并行计算\5.png" style="zoom:100%;" />

流水线技术的特点：

- **流水过程由多个相联系的子过程组成**，每个过程称为流水线的“级”或“段”，一条流水线的段数，也称为流水线的“深度”或“流水深度”。
  - 比如上例流水线深度为4.
- 每个子过程由专用的功能段实现。
- **各个功能段所需时间应尽量相等**，否则，时间长的功能会造成流水线的“堵塞”和“断流”，**时间最长的段将成为流水线的瓶颈**，这个时间一般为一个时钟周期（拍）或机器周期。
- 流水线需要有“通过时间”（第一个任务流出结果所需的时间），在此之后流水过程才进入稳定工作状态，每一个时钟周期（拍）流出一个结果。
- **流水线技术适合大量重复的时序过程**，只有在输入端能够连续地提供任务，流水线地效率才能充分发挥。
- 流水线的每段后面都要有一个缓冲寄存器，称为锁存器。其作用是在相邻的两段之间传送数据，以提供后面流水段要用到的信息。还有一个作用是隔离各段的工作，避免相邻流水段电路的相互打扰。

### 分类

静态流水线是指在同一时间内，多功能流水线中的各段只能按同一种功能的连接方式工作。当流水线要切换到另一种功能时，必须等前面的任务都流出流水线之后，才能改变连接。如下面的时空图中，必须要等到浮点加法的任务的功能段全部流出流水线之后，该流水线才能够改变功能，让定点乘法的流水线进入。

<img src="多核平台上的并行计算\6.png" style="zoom:100%;" />

动态流水线是指在同一时间内，多功能流水线的各段可以按照不同的方式连接，同时执行多种功能的流水线。它允许在某些段正在实现某些运算时，另一些段却在实现另一种运算。当然，前提是，**任何一个功能段只能参加到一种连接中**。如下图所演示的那样。动态流水线的优点是：更加灵活，能提高各段的使用率，能提高处理速度，但其控制复杂度增加了。

<img src="多核平台上的并行计算\7.png" style="zoom:100%;" />

由于动态流水线的控制很复杂，所以目前大多数都是静态流水线。


###   **流水线性能分析**

- 吞吐率（Throughput Rate）

单位时间内流水线所完成的额任务数或输出结果的数量，包括最大吞吐率和实际吞吐率两种指标。

（1）最大吞吐率 $TP_{max}$ 指的是流水线在连续流动达到稳定状态后所得到的吞吐率。

其取决于流水线中最慢的一段所需要的时间。也就成为了流水线的瓶颈。

$TP_{max}=\dfrac{1}{max\{Δt_i\}}$

<img src="多核平台上的并行计算\8.png" style="zoom:80%;" />

<img src="多核平台上的并行计算\12.png" style="zoom:80%;" />

为了解决瓶颈，可以**将瓶颈段再细分**。

<img src="多核平台上的并行计算\10.png" style="zoom:80%;" />

当瓶颈段不能再细分时，可**设置重复瓶颈段**，使其并行工作。这时最大吞吐率仍然能达到 $\dfrac{1}{Δt_0}$ ，但是各并行段之间的任务分配和同步都比较复杂。

<img src="多核平台上的并行计算\11.png" style="zoom:80%;" />

<img src="多核平台上的并行计算\9.png" style="zoom:80%;" />

对于指令流水线而言，流水线增加了吞吐率但是**并不会真正减少一条指令总的执行时间**。指令吞吐率的提高意味着**即使单条指令执行并没有变快，但是程序总执行时间将会有效减少**。

（2）实际吞吐率。

若各段时间相等，则实际吞吐率为：

$TP=\dfrac{n}{T_{流水}}=\dfrac{n}{mΔt_0+(n−1)Δt_0}=\dfrac{1}{(1+\dfrac{m−1}{n})Δt0}=\dfrac{TP_{max}}{1+\dfrac{m−1}{n}}$

可以看出，实际吞吐率小于最大吞吐率，与段数m、任务数n都有关。

若各段时间不等，实际吞吐率为：

$TP=\dfrac{n}{∑_{i=1}^mΔt_i+(n−1)Δt_i}$

其中 $Δt_i$ 为最慢一段所需时间

#### **加速比（Speedup Ratio）**

流水线的加速比是指**m段流水线的速度与等功能的非流水线的速度之比**。对于连续完成n个任务来讲，如果流水线各段时间相等，则其加速比为：

$S=\dfrac{T}{T_0}=\dfrac{n×m×Δt_0}{mΔt_0+(n−1)Δt_0}=\dfrac{m}{1+\dfrac{m−1}n}$

若各段时间不等，则

$S=\dfrac{n∑_{i=1}^mΔt_i}{∑_{i=1}^mΔt_i+(n−1)Δt_j}$

#### **效率（Efficiency）**

流水线的设备利用率。

若各段时间相等，则各段的效率ei是相等的，即：

$e_0=e_1=e_2=⋯=e_m=\dfrac{nΔt_0}{T_{流水}}$

整个流水线的效率为：

$E=\dfrac{e_1+e_2+⋯+e_m}{m}=\dfrac{me_0}{m}=\dfrac{mnΔt_0}{mT_{流水}}$

从时空图上看，所谓效率即n个任务占用的时空区和m个段总的时空区之比，即**效率含有时间和空间两个方面的因素**。

因为 $E=\dfrac{nΔt_0}{T_{流水}}=\dfrac{n}{m+n−1}=\dfrac{1}{1+\dfrac{m−1}n}$

则 $E=\dfrac{S}{m}$

则效率为实际加速比（S）和最大加速比(m)之比

同时可以得出 $E=TPΔt_0$

即，当 t0 不变时，**流水线的效率和吞吐率成正比**。

若各段时间不等，则：

各任务占用的时空区个段总的时空区$E=\dfrac{n各任务占用的时空区}{m个段总的时空区}=\dfrac{n∑_{i=1}^m}{m[∑_{i=1}^mΔt_i+(n−1)Δt_j]}$

则 $E=TP\dfrac{∑_{i=1}^mΔt_i}m$

对于非线性流水线和多功能流水线，也可仿照上述对线性流水线的性能分析方法，在正确画出时空图的基础上，分析其吞吐率和效率等。

### 一条经典的5段流水线

<img src="多核平台上的并行计算\13.png" style="zoom:80%;" />

- 取指（IF）
  - 以程序计数器（PC）中的内容作为地址，从存储器中取出指令并放入指令寄存器（IR）；
    PC值加4（假设每条指令占4字节），指向顺序的下一条指令。

- 指令译码/读寄存器周期（ID）
  - 对指令进行译码，并用IR中的寄存器地址去访问通用寄存器组，读出所需的操作数；
    对IR中的立即数进行扩展

- 执行/有效地址计算周期（EX）
  - ALU对上一个周期中准备好的操作数进行运算或处理。在这个阶段，不同类型的指令进行的操作不同。
    - load和store指令：ALB把指令中所指定的寄存器的内容与偏移量相加，形成访存有效地址；
    - 寄存器-寄存器 ALU 指令：ALU按照操作码指定的操作对从通用寄存器组中读出的数据进行运算；
    - 寄存器-立即数 ALU 指令：ALU按照操作码指定的操作对从通用寄存器组中读出的操作数和指令中给出的立即数进行运算；
    - 分支指令：ALU把指令中给出的偏移量与PC值相加，形成转移目标的地址。同时，对在前一个周期读出的操作数进行判断，确定分支是否成功。

- 存储器访问/分支完成周期（MEM）
  - load和store指令：load指令根据上一个周期计算出的有效地址从存储器中读出的相应的数据；store把指定的数据写入这个有效地址对应的存储单元。
  - 分支指令：如果分支“成功”，就把前一个周期中计算好的转移目标地址送入PC。分支指令执行完成；否则，就不进行任何操作。

- 写回周期（WB）
  - 把结果写入通用寄存器组。对于ALU运算来说，这个结果来自ALU，而对于load指令来说，这个结果来自存储器。

还要解决流水处理带来的一些问题：

- 保证不会在同意时钟周期要求同一个功能段做两件不同的工作。例如，不能要求ALU既做有效地址计算，又同时做算术运算。RISC指令集比较简洁，所以该要求不难实现；
- 为了避免IF段的访存（取指令）与MEM段的访存（读写数据）发生冲突，必须采用分离的指令存储器和数据存储器，或者仍采用一个公用的存储器，但要采用分离的指令Cache和数据Cache，一般采用后者。
- ID段要对通用寄存器组进行读操作，而WB段要对通用寄存器组进行写操作，为了解决对同一通用寄存器的访问冲突，把写操作安排在时钟周期的前半拍完成，把读操作安排在后半拍完成。

流水线时空图的另外的画法：

<img src="多核平台上的并行计算\14.png" style="zoom:80%;" />

### 相关与流水线冲突

相关（Dependence）是指两条指令之间存在某种依赖关系。考虑两条指令`i`和`j`，`i`在`j`的前面。

- 数据相关

如果下述条件之一成立，则称指令`j`与指令`i`数据相关：

  （1）指令`j`使用指令`i`产生的结果；

  （2）指令`j`与指令`k`数据相关，而指令`k`又与指令`i`数据相关。

  其中第（2）个条件表明，数据相关具有传递性。例如，下面这一段代码存在数据相关：

<img src="多核平台上的并行计算\15.png" style="zoom:80%;" />

其中箭头表示必须保证的执行顺序。它由产生数据的指令指向使用该数据的指令。

当数据的流动是经过寄存器时，相关的检测比较直观和容易，因为统一寄存器在所有指令中的名称都是唯一的。而当数据的流动是经过存储器时，检测就比较复杂了，因为形式上相同的地址其有效地址未必相同，如某条指令中的10(R5) 与另一条指令中的10(R5)可能不同，因为R5的内容发生了变化；而形式不同的地址其有效地址却可能相同。

- 名相关

两条指令使用了相同的寄存器或存储单元的名字，但它们之间并没有数据流动。

（1）反相关：指令`j`所写的名与指令`i`所读的名相同

（2）输出相关：指令`j`所写的名与指令`i`所写的名相同。输出相关的指令的执行顺序不能颠倒，以保证最后的结果是指令`j`写进去的。

名相关的两条指令之间并没有数据流动，只是使用了相同的名而已。因此可以通过改变指令中操作数的名来消除。对于寄存器操作数进行换名称为寄存器换名（Register Renaming）。寄存器换名既可以用编译器静态实现，也可以用硬件动态完成。

考虑下述代码：

<img src="多核平台上的并行计算\16.png" style="zoom:90%;" />

- 控制相关

控制相关是指由分支指令引起的相关，它需要根据分支指令的执行结果来确定后面该执行哪个分支上的指令。

例如：

```text
if p1{
   S1;
}
S;
if p2{
   S2;
}
```

这里的`if p1`和`if p2`编译成目标代码以后都是分支指令。语句`S1`与`p1`控制相关，`S2`与`p2`控制相关，`S`与`p1`和`p2`均无关。

控制相关带来了两个限制：

- 与一条分支指令控制相关的指令不能被移到该分支之前；否则这些指令就不受该分支控制了。
  - 上面的例子，`S1`不能移到`p1`之前；`S2`不能移到`p2`之前；

- 不能把`S`移到`if`语句的`then`部分中。

#### 结构冲突

因硬件资源满足不了指令重叠执行的要求而发生的冲突；

例如，有些流水线处理机只有一个存储器，数据和指令都存放在这个存储器中。在这种情况下，当执行load指令需要存取数时，若又要同时完成其后某条指令的“取指令”，就会发生访存冲突，如下图带阴影的M所示。

<img src="多核平台上的并行计算\17.png" style="zoom:80%;" />

为了消除这种冲突，可以在前一条指令访问存储器时，将流水线停顿一个时钟周期，该停顿周期往往被称为“气泡”。用时空图来表示就是：

<img src="多核平台上的并行计算\18.png" style="zoom:80%;" />

<img src="多核平台上的并行计算\19.png" style="zoom:80%;" />

可以看出，为消除结构冲突引入的停顿将推迟流水线的完成时间，从而影响流水线的性能。**由于这种冲突出现的频度不低，因此一般是采用分别设置独立的指令存储器和数据存储器的方法。或者一个存储器，但是采用两个分离的Cache：指令Cache、数据Cache**。

#### 数据冲突

当指令在流水线中重叠执行时，因需要用到前面指令的结果而发生的冲突；

例如：

<img src="多核平台上的并行计算\20.png" style="zoom:80%;" />

DADD之后的所有指令都要用到DADD指令的计算结果，如图：

<img src="多核平台上的并行计算\21.png" style="zoom:80%;" />

DADD指令在其WB段（第5个时钟周期）才将计算结果写入寄存器R1，但是DSUB指令在其ID段（第3个时钟周期）就要从寄存器R1读取该结果，这就是一个数据冲突。若不采取措施防止这一情况发生，则DSUB指令读到的值就是错误的。XOR指令也受到这种冲突的影响，它在第4个时钟周期从R1读出的值也是错误的。而AND指令则可以正常执行，这是因为它在第5个时钟周期的后半拍才从寄存器读数据的，而DADD指令在第5个时钟周期的前半拍已将结果写入寄存器。

考虑两条指令i ii和j jj，且i ii在j jj之前进入流水线，可能发生的数据冲突有以下几种：

- **写后读冲突**（Read After Write, RAW）：指令j jj用到指令i ii的计算结果，而且在i ii将结果写入寄存器之前就去读该寄存器，因而得到的是旧值。对应于数据相关。

- **写后写冲突**（Wirte After Write, WAW）：指令j jj和指令i ii的结果寄存器相同，而且j jj在i ii写入之前就对该寄存器进行了写入操作，最后在结果寄存器中留下的是i ii写入的值，而不是j jj写入的值。对应于输出相关。
  写后写冲突仅发生在：流水线中不止一个段可以进行写操作；或者指令顺序被重新排列。前面的简单的5段流水线只在WB写了寄存器，所以不会发生写后写冲突。

- **读后写冲突**（Read After Write, RAW）：指令j jj的目的寄存器和指令i ii的源操作数寄存器相同，而且j jj在i ii读取该寄存器之前就先对它进行了写操作，导致i ii读到的值是错误的。这种冲突是由反相关引起的。这在前面所述的简单5段流水线中不会发生，因为该流水线所有的读操作（在ID段）都在写操作（WB段）之前发生。读后写冲突仅发生在这样的情况下：有些指令的写操作提前了，而有些指令的读操作滞后了；或者指令被重新排序了。

可以使用**定向技术**减少数据冲突引起的停顿。

对于前面例子中的数据冲突，简单的处理方法是暂停流水线中DADD之后的所有指令，直到DADD指令将计算结果写入寄存器R1之后，再让DADD之后的指令继续执行，但这种暂停会导致性能下降。

可以建立数据旁路，将数据直接从产生的地方送到使用的地方，如图：

<img src="多核平台上的并行计算\22.png" style="zoom:80%;" />

但并不是所有的数据冲突都可以通过定向技术来解决，如图：

<img src="多核平台上的并行计算\23.png" style="zoom:80%;" />

DADD要使用LD指令的结果，如虚线所示，这显然是无法实现的。解决办法就是停顿。

##### 依靠编译器解决数据冲突

为了减少停顿，对于无法使用定向技术解决的数据冲突，可以通过在编译时让编译器重新组织指令顺序来消除冲突，这种技术称为“指令调度”(Instruction Scheduling)或“流水线调度”(Pipeline Scheduling)。

举个例子，考虑下面的表达式：
$$
A=B+C\\
D=E-F
$$


这两个式子编译后形成的代码为：

<img src="多核平台上的并行计算\24.png" style="zoom:80%;" />

在这个代码序列中，DADD Ra, Rb, Rc 与 LD Rc, C之间存在数据冲突，DSUB Rd, Re, Rf 与 LD Rf, F 之间存在数据冲突。为了保证流水线能正确执行，必须在指令执行过程中插入两个停顿周期（分别在DADD和DSUB执行前）。
将指令顺序调整为：

<img src="多核平台上的并行计算\25.png" style="zoom:80%;" />

再加上相应的数据旁路，就可以消除数据冲突，不必在执行过程中插入任何停顿的周期。

#### 控制冲突

流水线遇到分支指令或其它会改变PC值的指令所引起的冲突。

在流水线中，控制冲突可能会比数据冲突造成更多的性能损失，所以同样需要得到很好的处理。

执行分支指令的结果有两种：一种是分支“成功”，PC值改变为分支转移的目标地址；另一种则是“不成功”或者“失败”，这时PC值保持正常递增，指向顺序的下一条指令。对分支指令“成功”的情况来说，是在条件判定和转移地址计算都完成后，才改变PC值的。对于上面的5段流水线来说，改变PC值是在MEM段进行的。

处理分支指令最简单的方法是“冻结”或者“排空”流水线。即一旦在流水线的译码段ID检测到分支指令，就暂停其后所有指令的执行，直到分支指令到达MEM段、确定是否成功并计算出新的PC值为止。然后，按照新的PC值取指令，如下图所示：

<img src="多核平台上的并行计算\26.png" style="zoom:80%;" />

在这种情况下，分支指令给流水线带来了三个时钟周期的延迟。而且无论分支成功还是失败，都是一样的排空等待。

分支指令在目标代码中出现的频度是不低的，统计结果表明，每三四条指令就有一条是分支指令。所以排空的处理办法所带来的性能损失是相当大的。

另外采取的办法是：**提前判断（或猜测）分支是否成功，尽早计算出分支目标地址**。可以把它们提前到ID段完成（比如BEQZ，只需要判断指定寄存器中的数是否为0，ID段已经读取寄存器的同时就可以做，计算出分支目标地址是将PC+IMM<<2，ID段已经取出IR中的立即数，并且能做符号扩展，当然也能做加法），这样就可以把分支延迟减少到一个时钟周期。

##### 预测分支失败

方法是沿失败的分支继续处理指令，当确定分支失败时，就将分支指令看成一条普通指令，正常流动；当确定分支成功时，就把分支之后取出的指令转换为空操作，并按分支目标地址重新取指令执行。

##### 预测分支成功

当流水线检测到分支指令后，一旦计算出了分支目标地址，就开始从该目标地址取指令执行。

在上面的5段流水线中，由于判断分支是否成功与分支目标地址计算式在同一流水段完成的，所以这种方法对减少该流水段的分支延迟没有任何好处。但在其它的一些流水线处理机中，特别是哪些具有隐含设置条件码或分支条件更复杂（因而更慢）的流水线处理机中，在确定分支是否成功之前，就能得到分支的目标地址。这时采用这种方法便可以减少分支延迟。

注意：这和预测分支失败不是一种方法！

<img src="多核平台上的并行计算\27.png" style="zoom:80%;" />

##### 延迟槽

在分支指令后添加延迟槽：

分支指令

延迟槽

后继指令

具体延迟槽里放的指令是由编译器来选择的。实际上，延迟槽能否带来好处完全取决于编译器能否把有用的指令调度到延迟槽中，这也是一种指令调度技术。常用的调度方法有三种：从前调度、从目标处调度、从失败处调度。

举例，如图：

<img src="多核平台上的并行计算\28.png" style="zoom:80%;" />

a) 从前调度

把位于分支指令之前的一条独立的指令移到延迟槽。

b）从目标处调度

如图，分支指令之前数据存到了R1寄存器，然后R1寄存器被用于了分支指令的判断，所以不能把分支指令之前的这条指令调到延迟槽。

但是可以将分支指令跳转地址的第一条指令调到延迟槽，然后把分支跳转地址修改为第一条指令所在地址的后面，这实际上是预测了分支会成功，如果失败，相当于延迟槽中的指令就无用了。而且这要求分支地址计算在下一条指令进入之前完成。

c）从失败处调度

即将分支指令后面的第一条指令放入延迟槽，正常运行。

### 流水线实践

https://blog.csdn.net/PaddlePaddle/article/details/124600908

https://blog.csdn.net/weixin_38106322/article/details/107017721

&nbsp;

## 分支预测

在结果冒险和数据冒险中，可以发现，所有的流水线停顿的操作都要从**指令执行阶段**开始。流水线的前两个阶段，也就是取取指令（IF）和指令译码（ID）阶段，是不需要停顿的。CPU会在流水线里面去取下一条指令，然后进行译码。

取指令和指令译码不会需要遇到任何瓶颈，这是基于一个假设。这个假设就是，**所有的指令代码都是顺序加载执行的**。不过这个假设，在执行的代码中，一旦遇到`if…else`这样的条件分支，或者for/while循环，就会不成立。

![](多核平台上的并行计算\29.png)

从上面可以看到，在jmp指令发生的时候，CPU可能会跳转去执行其他指令。jmp后的那一条指令是否应该顺序加载执行，在流水线里面进行取指令的时候，我们没法知道。要等jmp指令执行完成，去更新了PC寄存器之后，我们才能知道，是否执行下一条指令，还是跳转到另外一个内存地址，去取别的指令。

**这种为了确保能够取到正确的指令，而不得不进行瞪眼延迟的情况，叫做控制冒险(控制冲突)**。这也是流水线设计里最后一种冒险。

### 分支预测：今天下雨了，明天还会继续下雨么？

在遇到了控制冒险之后，我们的 CPU 具体会怎么应对呢？除了流水线停顿，等待前面的 jmp 指令执行完成之后，再去取最新的指令，还有什么好办法吗？当然有

#### 缩短分支延迟

条件跳转指令其实进行了两种电路操作。

- 第一种，进行条件比较。这个条件比较，需要的输入是，根据指令的opcode，就能确认的条件码寄存器
- 第二种，进行实际的跳转，也就是要把跳转的地址信息写入PC寄存器

无论是 opcode，还是对应的条件码寄存器，还是我们跳转的地址，都是在指令译码（ID）的阶段就能获得的。而对应的条件码比较的电路，只要是简单的逻辑门电路就可以了，并不需要一个完整而复杂的ALU。

所以，我们可以将条件判断、地址跳转，都提前到指令译码阶段进行，而不需要放在指令执行阶段。对应的，我们也在在CPU里面设计对应的旁路，在指令译码阶段，就提供对应的判断条件比较的电路。

这种方式，本质上和前面数据冒险的操作数前推的解决方案类似，就是在硬件电路层面，把一些计算结果更早地反馈到流水线中。这样反馈变得更快了，后面的指令需要等待的时间就变短了。

不过只是改造硬件，并不能彻底解决问题。跳转指令的比较结果，仍然需要在指令执行的时候才能知道。在流水线里，第一条指令进行指令译码的时钟周期里，我们其实就要去取下一条指令了。这个时候，我们其实还没有开始指令执行阶段，自然也就不知道比较的结果。

#### 静态分支预测

所以，这个时候，我们就引入了一个新的解决方案，叫作分支预测（Branch Prediction）技术，也就是说，让CPU猜一猜，条件跳转后执行的指令，应该是“哪一条”。

最简单的分支预测技术，叫做**假装分支不发生**。也就是仍然按照顺序，把指令往下执行。这是一种静态预测技术。就好像猜硬币的时候，你一直猜正面，会有 50% 的正确率。

如果分支预测是正确的，我们自然赚到了。这个意味着，我们节省下来本来需要停顿下来等待的时间。如果分支预测失败了呢？那我们就把后面已经取出指令已经执行的部分，给丢弃掉。这个丢弃的操作，在流水线里面，叫做Zap和Flush。CPU不仅要执行后面的指令，对于这些已经在流水线里面执行到一半的指令，我们还需要做对应的清除操作。比如，清空已经使用的寄存器里面的数据等等，这些清除操作，也有一定的开销。

所以，CPU需要提供对应的丢弃指令的功能，通过控制信号清除掉已经在流水线中执行的指令。只要对应的清除开销不要太大，就是划得来的。

![](多核平台上的并行计算\31.png)

#### 动态分支预测

第三个办法，叫作动态分支预测。

上面的静态预测策略，看起来比较简单，预测的准确率也许有 50%。但是如果运气不好，可能就会特别差。于是，工程师们就开始思考，我们有没有更好的办法呢？比如，**根据之前条件跳转的比较结果来预测**，是不是会更准一点？

我们日常生活里，最经常会遇到的预测就是天气预报。如果没有气象台给你天气预报，你想要猜一猜明天是不是下雨，你会怎么办？

有一个简单的策略，就是完全根据今天的天气来猜。如果今天下雨，我们就预测明天下雨。如果今天天晴，就预测明天也不会下雨。这是一个很符合我们日常生活经验的预测。因为一般下雨天，都是连着下几天，不断地间隔地发生“天晴- 下雨 - 天晴 - 下雨”的情况并不多见。

而同样的策略，我们一样可以放在分支预测上。这种策略，叫做一级分支预测（One Level Branch Prediction），或者叫1 比特饱和计数（1-bit saturating counter）。这个方法，就是**用一个比特，去记录当前分支的比较情况，直接用当前分支的比较情况，来预测下一次分支时候的比较情况**。

只用一天下雨，就预测第二天下雨，这个方法还是有些“草率”，我们可以用更多的信息，而不只是一次的分支信息来进行预测。于是，我们可以引入一个**状态机**来做这个事情。

如果连续发生下雨的情况，我们就认为更有可能下雨。之后如果只有一天放晴了，我们仍然认为会下雨。在连续下雨之后，要连续两天放晴，我们才会认为之后会放晴。

![](多核平台上的并行计算\30.png)

这个状态机里，我们一共有 4 个状态，所以我们需要 2 个比特来记录对应的状态。这样这整个策略，就可以叫作2 比特饱和计数，或者叫双模态预测器

### 为什么循环嵌套的改变会影响性能呢？

```java
public class BranchPrediction {
    public static void main(String args[]) {        
        long start = System.currentTimeMillis();
        for (int i = 0; i < 100; i++) {
            for (int j = 0; j <1000; j ++) {
                for (int k = 0; k < 10000; k++) {
                }
            }
        }
        long end = System.currentTimeMillis();
        System.out.println("Time spent is " + (end - start));
                
        start = System.currentTimeMillis();
        for (int i = 0; i < 10000; i++) {
            for (int j = 0; j <1000; j ++) {
                for (int k = 0; k < 100; k++) {
                }
            }
        }
        end = System.currentTimeMillis();
        System.out.println("Time spent is " + (end - start) + "ms");
    }
}
```

这是一个简单的三重循环，里面没有任何逻辑代码。我们用两种不同的循环顺序各跑一次。第一次，最外重循环循环了100 次，第二重循环 1000 次，最内层的循环了 10000次。第二次，我们把顺序倒过来，最外重循环 10000 次，第二重还是 1000 次，最内层 100 次。

你可以猜一猜，这样两次运行，花费的时间是一样的么？结果应该会让你大吃一惊。我们可以看看对应的命令行输出。

```shell'
Time spent in first loop is 5ms
Time spent in second loop is 15ms
```

同样循环了十亿次，第一段程序只花了 5 毫秒，而第二段程序则花了 15 毫秒，足足多了 2 倍。

这个差异就是因为分支预测。循环本质也是利用cmp和jle这样的先比较后挑战的指令来实现的，每一次循环都有一个 cmp 和 jle 指令。每一个 jle 就意味着，要比较条件码寄存器的状态，决定是顺序执行代码，还是要跳转到另外一个地址。也就是说，在每一次循环的时候，都会有一次“分支”。



<img src="多核平台上的并行计算\32.png" style="zoom:70%;" />

从上图就可以看到，为什么同样空转次数相同的循环代码，第一段代码运行的时间要少得多了。因为第一段代码发生“分支预测”错误的情况比较少，更多的计算机指令，在流水线里顺序执行下去了，而不需要把运行到一半的指令丢弃掉，再去重新加载新的指令执行。



### 另一个例子：铁路分叉道口

```cpp
#include <algorithm>
#include <ctime>
#include <iostream>

int main()
{
    // 随机产生整数，用分区函数填充，以避免出现分桶不均
    const unsigned arraySize = 32768;
    int data[arraySize];

    for (unsigned c = 0; c < arraySize; ++c)
        data[c] = std::rand() % 256;

    // !!! 排序后下面的Loop运行将更快
    std::sort(data, data + arraySize);

    // 测试部分
    clock_t start = clock();
    long long sum = 0;

    for (unsigned i = 0; i < 100000; ++i)
    {
        // 主要计算部分，选一半元素参与计算
        for (unsigned c = 0; c < arraySize; ++c)
        {
            if (data[c] >= 128)
                sum += data[c];
        }
    }

    double elapsedTime = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;

    std::cout << elapsedTime << std::endl;
    std::cout << "sum = " << sum << std::endl;
}

//# 1. 取消std::sort(data, data + arraySize);的注释，即先排序后计算
//10.218
//sum = 312426300000

//# 2. 注释掉std::sort(data, data + arraySize);即不排序，直接计算
//29.6809
//sum = 312426300000
```

让我们回到19世纪，那个远距离无线通信还未普及的年代。你是铁路交叉口的扳道工。当听到火车快来了的时候，你无法猜测它应该朝哪个方向走。于是你叫停了火车，上前去问火车司机该朝哪个方向走，以便你能正确地切换铁轨。

要知道，火车是非常庞大的，切急速行驶时有巨大的惯性。为了完成上述停车-问询-切轨的一系列动作，火车需耗费大量时间减速，停车，重新开启。

既然上述过车非常耗时，那是否有更好的方法？当然有！当火车即将行驶过来前，你可以猜测火车该朝哪个方向走。

- 如果猜对了，它直接通过，继续前行。
- 如果猜错了，车头将停止，倒回去，你将铁轨扳至反方向，火车重新启动，驶过道口。

如果你不幸每次都猜错了，那么火车将耗费大量时间停车-倒回-重启。如果你很幸运，每次都猜对了呢？火车将从不停车，持续前行！

上述比喻可应用于处理器级别的分支跳转指令里：

原程序：

```text
if (data[c] >= 128)
    sum += data[c];   
```

汇编码：

```text
cmp edx, 128
jl SHORT $LN3@main
add rbx, rdx
$LN3@main:
```

让我们回到文章开头的问题。现在假设你是处理器，当看到上述分支时，当你并不能决定该如何往下走，该如何做？只能暂停运行，等待之前的指令运行结束。然后才能继续沿着正确地路径往下走。

要知道，现代编译器是非常复杂的，运行时有着非常长的pipelines， 减速和热启动将耗费巨量的时间。

那么，有没有好的办法可以节省这些状态切换的时间呢？你可以猜测分支的下一步走向！

- 如果猜错了，处理器要flush掉pipelines, 回滚到之前的分支，然后重新热启动，选择另一条路径。
- 如果猜对了，处理器不需要暂停，继续往下执行。

如果每次都猜错了，处理器将耗费大量时间在停止-回滚-热启动这一周期性过程里。如果侥幸每次都猜对了，那么处理器将从不暂停，一直运行至结束。

上述过程就是分支预测(branch prediction)。虽然在现实的道口铁轨切换中，可以通过一个小旗子作为信号来判断火车的走向，但是处理器却无法像火车那样去预知分支的走向--除非最后一次指令运行完毕。

那么处理器该采用怎样的策略来用最小的次数来尽量猜对指令分支的下一步走向呢？答案就是分析历史运行记录： 如果火车过去90%的时间都是走左边的铁轨，本次轨道切换，你就可以猜测方向为左，反之，则为右。如果在某个方向上走过了3次，接下来你也可以猜测火车将继续在这个方向上运行...

换句话说，你试图通过历史记录，识别出一种隐含的模式并尝试在后续铁道切换的抉择中继续应用它。这和处理器的分支预测原理或多或少有点相似。

大多数应用都具有状态良好的(well-behaved)分支，所以现代化的分支预测器一般具有超过90%的命中率。但是面对无法预测的分支，且没有识别出可应用的的模式时，分支预测器就无用武之地了。

关于分支预测期，可参考维基百科相关词条["Branch predictor" article on Wikipedia.](https://link.zhihu.com/?target=https%3A//en.wikipedia.org/wiki/Branch_predictor).

文首导致非排序数组相加耗时显著增加的罪魁祸首便是if逻辑：

```text
if (data[c] >= 128)
    sum += data[c];
```

注意到data数组里的元素是按照0-255的值被均匀存储的(类似均匀的分桶)。数组data有序时，前面一半元素的迭代将不会进入if-statement, 超过一半时，元素迭代将全部进入if-statement.

这样的持续朝同一个方向切换的迭代对分支预测器来说是非常友好的，前半部分元素迭代完之后，后续迭代分支预测器对分支方向的切换预测将全部正确。

简单地分析一下：有序数组的分支预测流程：

```text
T = 分支命中
N = 分支没有命中

data[] = 0, 1, 2, 3, 4, ... 126, 127, 128, 129, 130, ... 250, 251, 252, ...
branch = N  N  N  N  N  ...   N    N    T    T    T  ...   T    T    T  ...

       = NNNNNNNNNNNN ... NNNNNNNTTTTTTTTT ... TTTTTTTTTT  (非常容易预测)
```

无序数组的分支预测流程：

```text
data[] = 226, 185, 125, 158, 198, 144, 217, 79, 202, 118,  14, 150, 177, 182, 133, ...
branch =   T,   T,   N,   T,   T,   T,   T,  N,   T,   N,   N,   T,   T,   T,   N  ...

       = TTNTTTTNTNNTTTN ...   (完全随机--无法预测)
```

在本例中，由于data数组元素填充的特殊性，决定了分支预测器在未排序数组迭代过程中将有50%的错误命中率，因而执行完整个sum操作将会耗时更多。

#### 优化

利用位运算取消分支跳转。基本知识：

```text
|x| >> 31 = 0 # 非负数右移31为一定为0
~(|x| >> 31) = -1 # 0取反为-1

-|x| >> 31 = -1 # 负数右移31为一定为0xffff = -1
~(-|x| >> 31) = 0 # -1取反为0

-1 = 0xffff
-1 & x = x # 以-1为mask和任何数求与，值不变
```

故分支判断可优化为：

```text
int t = (data[c] - 128) >> 31; # statement 1
sum += ~t & data[c]; # statement 2
```

分析：

1. `data[c] < 128`, 则`statement 1`值为: `0xffff = -1`, `statement 2`等号右侧值为: `0 & data[c] == 0`;
2. `data[c] >= 128`, 则`statement 1`值为: 0, statement 2等号右侧值为: ~0 & data[c] == -1 & data[c] == 0xffff & data[c] == data[c];

故上述位运算实现的sum逻辑完全等价于if-statement, 更多的位运算hack操作请参见[bithacks](https://link.zhihu.com/?target=http%3A//graphics.stanford.edu/~seander/bithacks.html).

若想避免移位操作，可以使用如下方式：

```text
int t=-((data[c]>=128)); # generate the mask
sum += ~t & data[c]; # bitwise AND
```



### 总结

应对控制冒险的三个方式。

- 改造我CPU 功能，通过增加对应的电路的方式，来缩短分支带来的延迟
- 静态分支预测：假装分支不发生
- 动态分支预测：根据之前的比较结果预测本次结果
  - 一级分支预测：预测结果和上一次的条件跳转是一致的
  - 双模态预测器：通过一个状态机，多看了一步过去的跳转比较结果。

**将最可能进入的分支放到if中，而不是else中**

Intel处理器有预测分支单元，第一次进入一个分支的时候，由于没有历史信息可供参考，是否跳转取决于Static Predictor（静态预测器）的预测策略。通常静态预测器的预测策略是：向下跳转预测为不跳转，向上跳转预测为跳转。根据此特性，if语句也是需要按照这种方式去编码。

&nbsp;

## 向量化计算(vectorization)

向量化计算(vectorization)，也叫vectorized operation，也叫array programming，说的是一个事情：将多次for循环计算变成一次计算。

<img src="多核平台上的并行计算\33.png" style="zoom:80%;" />

上图中，左侧为vectorization，右侧为寻常的For loop计算。将多次for循环计算变成一次计算完全仰仗于CPU的SIMD指令集，SIMD指令可以在一条cpu指令上处理2、4、8或者更多份的数据。在Intel处理器上，这个称之为SSE以及后来的AVX，在Arm处理上，这个称之为NEON。

因此简单来说，向量化计算就是将一个loop——处理一个array的时候每次处理1个数据共处理N次，转化为vectorization——处理一个array的时候每次同时处理8个数据共处理N/8次。

### **vectorization如何让速度更快？**

我们以x86指令集为例，1997年，x86扩展出了MMX指令集，伴随着80-bit的vector寄存器，首开向量化计算的先河。 之后，x86又扩展出了SSE指令集 (有好几个版本, 从SSE1到SEE4.2)，伴随着128-bit寄存器。而在2011年，Intel发布了Sandy Bridge架构——扩展出了AVX指令集(256-bit寄存器)。在2016年，第一个带有AVX-512寄存器的CPU发布了(512-bit寄存器，可以同时处理16个32-bit的float数)。SSE和AVX各有16个寄存器。SSE的16个寄存器为XMM0-XMM15，AVX的16个寄存器为YMM0-YMM15。XMM registers每个为128 bits，而YMM寄存器每个为256bit（AVX512为512bit）。

SSE有3个数据类型：__m128 , __m128d 和 __m128i，分别代表Float、double (d) 和integer (i)。AVX也有3个数据类型： __m256 , __m256d 和 __m256i，分别代表Float、double (d) 和 integer (i)。

![img](https://pic3.zhimg.com/80/v2-94f7b921e07e7240fdf19601a9cda45a_720w.webp)

Gemfield使用下面一小段C++程序来展示一下AVX带来的计算速度：

```cpp
#include <immintrin.h>
#include <iostream>
#include <chrono>
#include <ctime> 

const int N = 8;
const int loop_num = 100000000;
float gemfield_i[8] = {1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8};
float gemfield_m[8] = {2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9};
float gemfield_a[8] = {11.1,12.2,13.3,14.4,15.5,16.6,17.7,18.8};
float gemfield_o[8] = {0};

__m256 gemfield_v_i = _mm256_set_ps(8.8,7.7,6.6,5.5,4.4,3.3,2.2,1.1);
__m256 gemfield_v_m = _mm256_set_ps(9.9,8.8,7.7,6.6,5.5,4.4,3.3,2.2);
__m256 gemfield_v_a = _mm256_set_ps(18.8,17.7,16.6,15.5,14.4,13.3,12.2,11.1);
__m256 gemfield_v_o = _mm256_set_ps(0,0,0,0,0,0,0,0);


void syszuxMulAndAddV() {
    auto start = std::chrono::system_clock::now();
    for(int j=0; j<loop_num; j++){
        gemfield_v_o += _mm256_fmadd_ps(gemfield_v_i, gemfield_v_m, gemfield_v_a);
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "resultV: ";
    // float* f = (float*)&gemfield_v_o;
    for(int i=0; i<N; i++){
        std::cout<<gemfield_v_o[i]<<" ";
    }
    std::cout<< "\nelapsed time: " << elapsed_seconds.count() << "s\n";
}

void syszuxMulAndAdd(){
    auto start = std::chrono::system_clock::now();
    for(int j=0; j<loop_num; j++){
        for(int i=0; i<N; i++){
            gemfield_o[i] += gemfield_i[i] * gemfield_m[i] + gemfield_a[i];
        }
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "result: ";
    for(int i=0; i<8; i++){
        std::cout<<gemfield_o[i]<<" ";
    }
    std::cout<< "\nelapsed time: " << elapsed_seconds.count() << "s\n";
}

int main() {
    syszuxMulAndAdd();
    syszuxMulAndAddV();
    return 0;
}
```

编译并运行：

```console
#compile civilnet.cpp
gemfield@ThinkPad-X1C:~$ g++ -march=skylake-avx512 civilnet.cpp -o civilnet

#run civilnet
gemfield@ThinkPad-X1C:~$ ./civilnet
result: 2.68435e+08 5.36871e+08 5.36871e+08 1.07374e+09 1.07374e+09 2.14748e+09 2.14748e+09 2.14748e+09 
elapsed time: 2.39723s
resultV: 2.68435e+08 5.36871e+08 5.36871e+08 1.07374e+09 1.07374e+09 2.14748e+09 2.14748e+09 2.14748e+09 
elapsed time: 0.325577s
```

for loop计算消耗了2.39723秒，而vectorization计算消耗了0.325577s，可以看到AVX的计算速度远超for loop，因为AVX使用了下面这样的并行方式：

<img src="多核平台上的并行计算\34.png" style="zoom:80%;" />

### Intel高级向量扩展

#### MMX指令

MultiMedia eXtensions（MMX），MMX指令主要使用的寄存器为MM0 ~ MM7，与浮点运算不能同时进行。MMX指令能一次性地操作1个64-bit的数据、或者两个32-bit的数据、或者4个16-bit的数据、或者8个8-bit的数据。MMX指令集的扩展包括：3DNow!、SSE、AVX

#### SSE指令

Streaming SIMD eXtensions（SSE），SSE指令采用了独立的寄存器组XMM0 ~ XMM7，64位模式下为XMM0 ~ XMM15，并且这些寄存器的长度也增加到了128-bit。SSE指令的升级版包括：SSE2/SSE3/SSSE3/SSE4

#### AVX/AVX2指令

Advanced Vector eXtentions（AVX），AVX对XMM寄存器做了扩展，从原来的128-bit扩展到了256-bit，并从XMM0–XMM7重命名为YMM0–YMM7，仍可通过SSE指令对YMM寄存器的低128位进行操作。新指令使用英特尔所谓的VEX前缀进行编码，这是一个两字节或三字节的前缀，旨在消除当前和将来的x86/x64指令编码的复杂性。AVX2将大多数整数命令扩展为256位，并引入了融合的乘加（FMA）操作。

#### FMA指令

Fused-Multiply-Add（FMA），FMA指令集是128-bit和256-bit的SSE的扩展指令集，以进行乘加运算。共有两种变体：FMA4、FMA3，自2014年以来，从PILEDRIVER架构开始，AMD处理器就支持FMA3；从Haswell处理器和Broadwell处理器开始，英特尔则支持FMA3。

#### AVX512指令

英特尔架构处理器支持旧式和现代指令集，从64位MMX扩展到新的512位指令AVX-512。ZMM的低256-bit与YMM混用。ZMM的前缀为EVEX。与AVX / AVX2相比，AVX-512最显着的新功能是512位矢量寄存器宽度。

#### 其他指令

KNC等其他指令集。

https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html

> 对于实验，我觉得可以考虑AVX/AVX2 

## SIMD汇总 — AVX/AVX2

```text
#include <mmintrin.h>  //MMX
#include <xmmintrin.h> //SSE(include mmintrin.h)
#include <emmintrin.h> //SSE2(include xmmintrin.h)
#include <pmmintrin.h> //SSE3(include emmintrin.h)
#include <tmmintrin.h> //SSSE3(include pmmintrin.h)
#include <smmintrin.h> //SSE4.1(include tmmintrin.h)
#include <nmmintrin.h> //SSE4.2(include smmintrin.h)
#include <wmmintrin.h> //AES(include nmmintrin.h)
#include <immintrin.h> //AVX(include wmmintrin.h)
#include <zmmintrin.h> //AVX512
#include <intrin.h>    //(include immintrin.h)
```

头文件中，下一个头文件包含上一个头文件中内容，例如`xmmintrin.h`为SSE 头文件，此头文件里包含MMX头文件，`emmintrin.h`为SSE2头文件，此头文件里包含SSE头文件。

### SIMD数据类型简介

SIMD数据类型有——（m前面是两个下划线）

`__m64`：64位紧缩整数（MMX）。

`__m128`：128位紧缩单精度（SSE）。

`__m128d`：128位紧缩双精度（SSE2）。

`__m128i`：128位紧缩整数（SSE2）。

`__m256`：256位紧缩单精度（AVX）。

`__m256d`：256位紧缩双精度（AVX）。

`__m256i`：256位紧缩整数（AVX）。

注：紧缩整数包括了8位、16位、32位、64位的带符号和无符号整数。

每一种类型，从2个下划线开头，接一个m，然后是vector的位长度。

如果向量类型是以d结束的，那么向量里面是double类型的数字。如果没有后缀，就代表向量只包含float类型的数字。

整形的向量可以包含各种类型的整形数，例如`char`，`short`，`unsigned long long`。也就是说，`__m256i`可以包含32个char，16个short类型，8个int类型，4个long类型。这些整形数可以是有符号类型也可以是无符号类型。

这些数据类型与寄存器的对应关系为：

- 64位MM寄存器（MM0~MM7）：`__m64`。
- 128位SSE寄存器（XMM0~XMM15）：`__m128`、`__m128d`、`__m128i`。
- 256位AVX寄存器（YMM0~YMM15）：`__m256`、`__m256d`、`__m256i`。

在MMX指令集中，使用的寄存器称作MM0到MM7，实际上借用了浮点处理器的8个寄存器的低64Bit，这样导致了浮点运算速度降低。

SSE指令集推出时，Intel公司在Pentium III CPU中增加了8个128位的SSE指令专用寄存器，称作XMM0到XMM7。这样SSE指令寄存器可以全速运行，保证了与浮点运算的并行性。这些XMM寄存器用于4个单精度浮点数运算的SIMD执行，并可以与MMX整数运算或x87浮点运算混合执行。

2001年在Pentium 4上引入了SSE2技术，进一步扩展了指令集，使得XMM寄存器上可以执行8/16/32位宽的整数SIMD运算或双精度浮点数的SIMD运算。对整型数据的支持使得所有的MMX指令都是多余的了，同时也避免了占用浮点数寄存器。SSE2为了更好地利用高速寄存器，还新增加了几条寄存指令，允许程序员控制已经寄存过的数据。这使得 SIMD技术基本完善。

SSE3指令集扩展的指令包含寄存器的局部位之间的运算，例如高位和低位之间的加减运算；浮点数到整数的转换，以及对超线程技术的支持。

AVX是Intel的SSE延伸架构，把寄存器XMM 128bit提升至YMM 256bit，以增加一倍的运算效率。此架构支持了三运算指令（3-Operand Instructions），减少在编码上需要先复制才能运算的动作。在微码部分使用了LES LDS这两少用的指令作为延伸指令Prefix。AVX的256bit的YMM寄存器分为两个128bit的lanes，AVX指令并不支持跨lanes的操作。其中YMM寄存器的低128位与Intel SSE指令集的128bitXMM寄存器复用。尽管VGX并不要求内存对齐，但是内存对齐有助于提升性能。如对于128-bit访问的16字节对齐和对于256-bit访问的32字节对齐。

AVX虽然已经将支持的SIMD数据宽度增加到了256位，但仅仅增加了对256位的浮点SIMD支持，整点SIMD数据的宽度还停留在128位上，AVX2支持的整点SIMD数据宽度从128位扩展到256位。同时支持了跨lanes操作，加入了增强广播、置换指令支持的数据元素类型、移位操作对各个数据元素可变移位数的支持、跨距访存支持。AVX硬件由16个256bitYMM寄存器（YMM0~YMM15）组成。

**每一代的指令集都是对上一代兼容的，支持上一代的指令，也可以使用上一代的寄存器**，也就是说，AVX2也依然支持128位，64位的操作，也可以使用上一代的寄存器（当然，寄存器的硬件实现可能有区别）。AVX也对部分之前的指令接口进行了重构，所以可以在指令文档中找到几个处于不同代际有着相同功能调用接口却不相同的函数。

另外，**不同代际的指令不要混用，每次状态切换将消耗 50-80 个时钟周期，会拖慢程序的运行速度**。

| MMX    | SSE     | SSE2      | AVX       | AVX2       |            |
| ------ | ------- | --------- | --------- | ---------- | ---------- |
| 寄存器 | MM0-MM7 | XMM0-XMM7 | XMM0-XMM7 | YMM0-YMM15 | YMM0-YMM15 |
| 浮点   |         | 128bit    | 128bit    | 256bit     | 256bit     |
| 整型   | 64bit   |           | 128bit    | 128bit     | 256bit     |



#### 函数命名约定

```c
_mm<bit_width>_<name>_<data_type>
```

- `<bit_width>` 表明了向量的位长度，对于128位的向量，这个参数为空，对于256位的向量，这个参数为256。

- `<name>`描述了内联函数的算术操作。

- `<data_type>` 标识函数主参数的数据类型。
  - `ps` 包含float类型的向量
  - `pd` 包含double类型的向量
  - `epi8/epi16/epi32/epi64` 包含8位/16位/32位/64位的有符号整数
  - `epu8/epu16/epu32/epu64` 包含8位/16位/32位/64位的无符号整数
  - `si128/si256` 未指定的128位或者256位向量
  - `m128/m128i/m128d/m256/m256i/m256d` 当输入向量类型与返回向量的类型不同时，标识输入向量类型

#### 变量命名规范参考

参考匈牙利命名法（Hungarian notation），**在变量名前面增加类型前缀**。类型前缀为3个小写字母，首字母代表寄存器宽度，最后两个字母代表紧缩数据类型。

##### 寄存器宽度（首字母）

`m`：64位MM寄存器。对应 `__m64`

`x`：128位SSE寄存器。对应 `__m128`、`__m128d`、`__m128i`。

`y`：256位AVX寄存器。对应 `__m256`、`__m256d`、`__m256i`。

##### 紧缩数据类型（两个字母）

`mb`：8位数据。用于只知道长度、不知道具体紧缩格式时。（b：Byte）

`mw`：16位数据。（w：Word）

`md`：32位数据。（d：DoubleWord）

`mq`：64位数据。（q：QuadWord）

`mo`：128位数据。（o：OctaWord）

`mh`：256位数据。（h：HexWord）

`ub`：8位无符号整数。

`uw`：16位无符号整数。

`ud`：32位无符号整数。

`uq`：64位无符号整数。

`ib`：8位带符号整数。

`iw`：16位带符号整数。

`id`：32位带符号整数。

`iq`：64位带符号整数。

`fh`：16位浮点数，即半精度浮点数。（h：Half）

`fs`：32位浮点数，即单精度浮点数。（s：Single）

`fd`：64位浮点数，即双精度浮点数。（d：double）

##### 样例

`mub`：64位紧缩字节（64位MMX寄存器，其中存放了8个8位无符号整数）。

`xfs`：128位紧缩单精度（128位SSE寄存器，其中存放了4个单精度浮点数）。

`xid`：128位紧缩带符号字（128位SSE寄存器，其中存放了4个32位带符号整数）。

`yfd`：256位紧缩双精度（256位AVX寄存器，其中存放了4个双精度浮点数）。

`yfh`：256位紧缩半精度（256位AVX寄存器，其中存放了16个半精度浮点数）。

#### 内存对齐

为了方便CPU用指令对内存进行访问，通常**要求某种类型对象的地址必须是某个值K（通常是2、4或8）的倍数**，如果一个变量的内存地址正好位于它长度的整数倍，我们就称他是自然对齐的。不同长度的内存访问会用到不同的汇编指令，这种对齐限制简化了形成处理器和存储器系统之间接口的硬件设计，提高了内存的访问效率。

通常对于各种类型的对齐规则如下：

- 数组 ：按照基本数据类型对齐，第一个对齐了后面的自然也就对齐了。

- union：按其包含的长度最大的数据类型对齐。

- 结构体： 结构体中每个数据类型都要对齐

对于SIMD的内存对齐是指`__m128`等union在内存中存储时的存储方式。然而由于结构体内存对齐的规则略微复杂，我们以结构为例进行说明：

一般情况下，由于内存对齐的原因，存储多种类型数据的结构体所占的内存大小并非元素本身类型大小之和。对于自然对齐而言：

- 对于各**成员变量**来说，存放的起始地址相对于结构的起始地址的**偏移量**必须为**该变量的类型所占用的字节数的倍数**，各成员变量在存放的时候根据在结构中出现的顺序依次申请空间， 同时按照上面的对齐方式调整位置， 空缺的字节自动填充。

- 对于**整个结构体**来说，为了确保结构的大小为结构的字节边界数（即该结构中占用最大的空间的类型的字节数）的倍数，所以在为最后一个成员变量申请空间后，还会根据需要自动填充空缺的字节。

所以一般我们在定义结构体时定义各元素的顺序也会影响实际结构体在存储时的整体大小，**把大小相同或相近的元素放一起，可以减少结构体占用的内存空间**。

除了自然对齐的内存大小，我们也可以设置自己需要的对齐大小，我们称之为对齐系数，如果结构内最大类型的字节数小于对齐系数，结构体内存大小应按最大元素大小对齐，如果最大元素大小超过对齐系数，应按对齐系数大小对齐。

对齐系数大小的设定可以使用下列方法：

`#pragma pack (16)`使用预编译器指令要求对齐。`#pragma pack()`恢复为默认对齐方式。

```cpp
__attribute__ ((aligned (16)))//GCC要求对齐
```

```cpp
__declspec(intrin_type) _CRT_ALIGN(16)//Microsoft Visual C++要求对齐
```

联合的内存对齐方式与结构类似。

SIMD的指令中通常有对内存对齐的要求，例如，SSE中大部分指令要求地址是16bytes对齐的，以`_mm_load_ps`函数来说明，这个函数对应于SSE的loadps指令。

函数原型为：`extern __m128 _mm_load_ps(float const*_A);`

可以看到，它的输入是一个指向float的指针，返回的就是一个`__m128`类型的数据，从函数的角度理解，就是把一个float数组的四个元素依次读取，返回一个组合的`__m128`类型的SSE数据类型，从而可以使用这个返回的结果传递给其它的SSE指令进行运算，比如加法等；从汇编的角度理解，它对应的就是读取内存中连续四个地址的float数据，将其放入SSE的寄存器(XMM)中，从而给其他的指令准备好数据进行计算。其使用示例如下：

```cpp
float input[4] = { 1.0f, 2.0f, 3.0f, 4.0f };
__m128 a = _mm_load_ps(input);	//WARNING
```

这里加载正确的前提是：input这个浮点数阵列都是对齐在16 bytes的边上。否则程序会崩溃或得不到正确结果。如果没有对齐，就需要使用`_mm_loadu_ps`函数，这个函数用于处理没有对齐在16 bytes上的数据，但是其速度会比较慢。

对于上面的例子，如果要将input指定为16 bytes对齐，可以采用的方式是：

```cpp
__declspec(align(16)) float input[4] = {1.0, 2.0, 3.0, 4.0};
```

动态数组（dynamic array）可由`_aligned_malloc`函数为其分配空间：

```cpp
 input = (float*) _aligned_malloc(ARRAY_SIZE * sizeof(float), 16);
```

由`_aligned_malloc`函数分配空间的动态数组可以由`_aligned_free`函数释放其占用的空间：

```cpp
_aligned_free(input);
```

256-bit AVX 指令在内存访问上对内存对齐比128-bit SSE 指令有更高要求。虽然在一个cache-line 之内，Intel 的对齐和非对齐指令已经没有性能差距了，但是**由于AVX 有更长的内存访问宽度（YMM <-> memory），会更频繁地触及cache-line 边界**。所以1）尽量使用对齐内存分配；2）有时候内存对齐不能保证，可以用128-bit（XMM）指令访问内存，然后再组合成256-bit YMM。

##### **通知编译器关于数据对齐**

https://www.intel.cn/content/www/cn/zh/developer/articles/technical/data-alignment-to-assist-vectorization.html

从**矢量对齐和循环并行化**开始时openMP。

&nbsp;

### MMX

[MMX指令集介绍](https://blog.csdn.net/yangjianqiao0/article/details/69388595?spm=1001.2101.3001.6650.6&utm_medium=distribute.pc_relevant.none-task-blog-2%7Edefault%7EBlogCommendFromBaidu%7ERate-6-69388595-blog-1001849.pc_relevant_default&depth_1-utm_source=distribute.pc_relevant.none-task-blog-2%7Edefault%7EBlogCommendFromBaidu%7ERate-6-69388595-blog-1001849.pc_relevant_default&utm_relevant_index=9)

### AVX

https://blog.csdn.net/chen134225/article/details/105935153



### 环境配置

使用软件CPU-Z可以查看CPU支持的指令集。

#### 编译器设置

 我们可以在C/C++使用封装的函数而不是嵌入的汇编代码的方式来调用指令集，这就是Compiler Intrinsics。

 Intrinsics指令是对MMX、SSE等指令集的指令的一种封装，以函数的形式提供，使得程序员更容易编写和使用这些高级指令，在编译的时候，**这些函数会被内联为汇编**，不会产生函数调用的开销。

 除了我们这里使用的intrinsics指令，还有intrinsics函数需要以作区分，这两者既有联系又有区别。编译器指令`#pragma intrinsic()`可以**将一些指定的系统库函数编译为内部函数，从而去掉函数调用参数传递等的开销，这种方式只适用于编译器规定的一部分函数，不是所有函数都能使用，同时会增大生成代码的大小**。

 intrinsics更广泛的使用是指令集的封装，能将函数直接映射到高级指令集，同时隐藏了寄存器分配和调度等，从而使得程序员可以以函数调用的方式来实现汇编能达到的功能，编译器会生成为对应的SSE等指令集汇编。

Intel Intrinsic Guide可以查询到所有的Intrinsic指令、对应的汇编指令以及如何使用等。

 对于VC来说，VC6支持MMX、3DNow!、SSE、SSE2，然后更高版本的VC支持更多的指令集。但是，VC没有提供检测Intrinsic函数集支持性的办法。

而对于GCC来说，它使用`-mmmx`、`-msse`等编译器开关来启用各种指令集，同时定义了对应的 `__MMX__`、`__SSE__`等宏，然后`x86intrin.h`会根据这些宏来声明相应的Intrinsic函数集。`__MMX__`、`__SSE__`等宏可以帮助我们判断Intrinsic函数集是否支持，但这只是GCC的专用功能。

如果使用GCC编译器时，使用intrinsics指令时需要在编写cmake或者makefile文件时加上相关参数，例如使用AVX指令集时添加`-mavx2`参数。

<img src="多核平台上的并行计算\37.png" style="zoom:100%;" />

&nbsp;

## 循环展开

循环展开是通过增加每次迭代计算的元素的数量，减少循环的迭代次数。循环展开只能针**对整形加法和乘法**的性能改进。循环展开是一种牺牲程序的尺寸来加快程序的执行速度的优化方法。可以由程序员完成，也可由编译器自动优化完成。循环展开最常用来降低循环开销，为具有多个功能单元的处理器提供指令级并行。也有利于指令流水线的调度。

循环展开从两个方面改变程序的性能：

- 分支预测失败减少。
- 减少不直接有助于程序结果的操作的数量，如循环索引计算和条件分支。

循环展开对于性能的提升是由帮助的，但这种帮助并不是无限的，**随着展开次数的增多，性能并不会继续增加，相反，循环展开次数过多，会使得程序代码膨胀、代码可读性降低。另外，编译器优化选项`-O1`或`-O2`等，会使得编译器自身会对代码进行优化**，此时手动循环展开并不是一个好的方法。

这里有一篇向量程序中应用循环展开的文献，具体目的是确定最优展开次数，但是并无具体算法，比较难复现。[链接](https://www.jsjkx.com/EN/article/openArticlePDF.jsp?id=16243)

### 循环展开优化程序性能的原因

真正的循环展开，在《计算机体系结构：量化研究方法》这本书里的讲述是因为流水线具有指令调度能力，循环展开有利于流水线的充分调度。

我们可以考虑计算和访存的矛盾。汇编语言中就是取指令（load）和计算（ALU）的矛盾，访存和计算天然是两种独立的资源，两者可以独立运行互不影响。

指令调度的能力，明线是利用CPU流水线的特性，保证CPU满载，所以要找到不重叠的指令减少流水线停顿。在我看来还有一条暗线，利用访存和计算两者可以并行的特性，在计算A的时候把A+1需要的数据先取过来，然后就可以在下一个时钟周期紧锣密鼓的计算A+1。

当循环次数变少，循环体增多的时候，CPU在执行循环体的过程中就有了更大的操作空间进行指令调度，在计算`a[i]`的时候就可以去调度`a[i+1]`相关的指令（取数据）和`a[i-1]`相关的指令（存数据），因此CPU流水线更满载，也就提高了性能。

同样，循环展开带来的坏处也显而易见，除了我在面试中回答的会增加代码体积，降低可读性之外，循环展开会导致在循环体内占用更多的寄存器，因此不能无限展开循环得到越来越好的效果。

### 循环程序并行化的一般方法

https://blog.csdn.net/qq_40765537/article/details/106098554




## 参考

> #### 流水线
>
> https://zhuanlan.zhihu.com/p/297100634
>
> https://blog.csdn.net/weixin_40064300/article/details/124415855
>
> #### 分支预测
>
> https://blog.csdn.net/zhizhengguan/article/details/121269908
>
> https://matt33.com/2020/04/16/cpu-branch-predictor/
>
> #### 向量化
>
> https://zhuanlan.zhihu.com/p/72953129
>
> https://zhuanlan.zhihu.com/p/337756824
>
> https://www.intel.cn/content/www/cn/zh/developer/articles/technical/common-vectorization-tips.html
>
> https://www.intel.cn/content/www/cn/zh/developer/articles/technical/data-alignment-to-assist-vectorization.html
>
> https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#techs=MMX 
>
> https://blog.csdn.net/qq_32916805/article/details/117637192?utm_medium=distribute.pc_relevant.none-task-blog-2~default~baidujs_baidulandingword~default-0-117637192-blog-103905183.pc_relevant_landingrelevant&spm=1001.2101.3001.4242.1&utm_relevant_index=3
>
> AVX编程基础：https://blog.csdn.net/chen134225/article/details/105935153
>
> https://zhuanlan.zhihu.com/p/94649418
>
> #### openMP-SIMD
>
> 编译器优化与SIMD指令集：https://zhuanlan.zhihu.com/p/483476732
>
> 入门知识：https://zhuanlan.zhihu.com/p/32506707
>
> #pragma openmp simd：https://blog.csdn.net/10km/article/details/84579465
>
> #### 循环展开
>
> http://lazybing.github.io/blog/2019/04/17/loop-unroll/
>
> https://zhuanlan.zhihu.com/p/37582101