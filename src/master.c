#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <crts.h>
#include <sys/time.h>
#include "../include/struct.h"

#define TYPE double

#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))

extern void SLAVE_FUN(struct_sptrsv_cross_entry)(void *);
extern void SLAVE_FUN(struct_sptrsv_entry)(void *);
extern void SLAVE_FUN(struct_sptrsv_2skew_entry)(void *);
extern void SLAVE_FUN(benchmark_entry)(void *);

//test
void init_and_rand(TYPE *a, int size)
{
    int i;
    //*a = (TYPE*)malloc(size * sizeof(TYPE));
    for (i = 0; i < size; i++)
    {
        a[i] = (TYPE)rand() / RAND_MAX / 3;
    }
}

void init_and_zero(TYPE *a, int size)
{
    int i;
    //*a = (TYPE*)malloc(size * sizeof(TYPE));
    for (i = 0; i < size; i++)
    {
        a[i] = 0;
    }
}

int check(TYPE *a, TYPE *b, int size)
{
    int i;
    for (i = 0; i < size; i++)
        if (a[i] != b[i]) return i;
    return -1;
}

//MPE.c
int nonzero_analysis_f(int nnz, const int (*pos)[3], int *skew, int *halo, int (*p)[2], int *col_m, int *off_m)
{
    if (nnz <= 0) return -1;

    int i, j;
    *skew = 0;
    *halo = 0;
    int cols = 0;
    for (i = 0; i < nnz; i++)
    {
        int x = pos[i][0],
        y = pos[i][1],
        z = pos[i][2];
        int tmp = x * 100 + y * 10 + z;
        if (tmp >= 0) return -1;
        if (abs(x) > 3 || abs(y) > 3 || abs(z) > 3) return -1;
        if (y > *skew) *skew = y;
        if (abs(z) > *halo) *halo = abs(z);
    }
    
    int flag[11][11];
    for (i = 0; i < 11; i++)
    for (j = 0; j < 11; j++)
        flag[i][j] = 0;

    for (i = 0; i < nnz; i++)
    {
        int x = pos[i][0],
        y = pos[i][1],
        z = pos[i][2];
        if (!flag[x+5][y+5])
        {
            flag[x+5][y+5] = 1;
        }
    }

    for (i = 0; i < 11; i++)
    for (j = 0; j < 11; j++)
    {
        if (flag[i][j])
        {
            if (i < 5 || i == 5 && j < 5)
            {
                p[cols][0] = i - 5;
                p[cols][1] = j - 5;
                cols++;
            }
        }
    }

    for (i = 0; i < nnz; i++)
    {
        if (pos[i][0] == 0 && pos[i][1] == 0)
        {
            col_m[i] = cols;
            off_m[i] = pos[i][2];
        }
        else
            for (j = 0; j < cols; j++)
            {
                if (pos[i][0] == p[j][0] && pos[i][1] == p[j][1])
                {
                    col_m[i] = j;
                    off_m[i] = pos[i][2];
                    break;
                }
            }
    }
    return cols;
}

void show(int nnz, int *col_id, int syc, int *syc_col)
{
    int i;
    printf("col_id\n");
    for (i = 0; i < nnz; i++)
        printf("%d: %d\n", i, col_id[i]);
    printf("syc_col %d\n", syc);
    for (i = 0; i < syc; i++)
        printf("%d: %d\n", i, syc_col[i]);
}

int struct_sptrsv_forward(const TYPE* a, TYPE *x, const TYPE* b, const int dim[3], int nnz, const int (*pos)[3])
{
    int i;
    int col_id[nnz], z_off[nnz];
    int p[nnz][2];
    int skew, halo;
    int cols;
    if ((cols = nonzero_analysis_f(nnz, pos, &skew, &halo, p, col_id, z_off)) < 0)
    {
        printf("nnz error!\n");
        return -1;
    }

    //show(nnz, col_id, syc, syc_col);

    struct_info info;
    info.a = a;
    info.x = x;
    info.b = b;
    info.M = dim[0];
    info.N = dim[1];
    info.L = dim[2];
    info.nnz = nnz;
    info.p = (int*)p;
    info.col_id = col_id;
    info.z_off = z_off;
    info.halo = halo;
    info.cols = cols;

    double read, write;
    //penv_mc_band_init();

    struct timeval st, ed;
    
    gettimeofday(&st, NULL);
    
    if (skew == 0) CRTS_athread_spawn(struct_sptrsv_cross_entry, &info);
    if (skew == 1) CRTS_athread_spawn(struct_sptrsv_entry, &info);
    if (skew == 2) CRTS_athread_spawn(struct_sptrsv_2skew_entry, &info);
    CRTS_athread_join();
    CRTS_athread_halt();

    gettimeofday(&ed, NULL);
    printf("Heterogenous finish! Time: %.6f\n", TIME(st, ed));

    return 0;
}

int struct_sptrsv_sequential_forward(const TYPE* a, TYPE *x, const TYPE* b, const int dim[3], int nnz, const int (*pos)[3])
{
    int i, j, k;
    int ti, tj, tk;
    int u;
    int M = dim[0], N = dim[1], L = dim[2];
    for (i = 0; i < M; i++)
    for (j = 0; j < N; j++)
    for (k = 0; k < L; k++)
    {
        TYPE tmp = b[INDEX3(i, j, k, N, L)];
        int aoff = nnz * INDEX3(i, j, k, N, L);
        for (u = 0; u < nnz; u++)
        {
            ti = i + pos[u][0];
            tj = j + pos[u][1];
            tk = k + pos[u][2];
            if (VALID(ti, tj, tk, M, N, L)) tmp -= a[aoff + u] * x[INDEX3(ti, tj, tk, N, L)];
        }
        x[INDEX3(i, j, k, N, L)] = tmp;
    }
}

int main()
{
    int i, nnz, M, N, L;
    FILE *fp = fopen("../input.struct", "r");
    fscanf(fp, "%d%d%d%d", &M, &N, &L, &nnz);

    int p[nnz][3];

    for (i = 0; i < nnz; i++)
    {
        fscanf(fp, "%d%d%d", &p[i][0], &p[i][1], &p[i][2]);
    }

    fclose(fp);

    int dim[3];
    dim[0] = M;
    dim[1] = N;
    dim[2] = L;

    __attribute__ ((aligned(256))) TYPE *a;
    a = (TYPE*)malloc(M * N * L * nnz * sizeof(TYPE));
    __attribute__ ((aligned(256))) TYPE x1[M * N * L];
    __attribute__ ((aligned(256))) TYPE x2[M * N * L];
    __attribute__ ((aligned(256))) TYPE b[M * N * L];
    
    srand(time(NULL));
    init_and_rand(a, M * N * L * nnz);
    init_and_zero(x1, M * N * L);
    init_and_zero(x2, M * N * L);
    init_and_rand(b, M * N * L);

    CRTS_init();
    struct timeval st, ed;
    double time1, time2;
    gettimeofday(&st, NULL);
    if (struct_sptrsv_forward(a, x1, b, dim, nnz, p)) return -1;
    gettimeofday(&ed, NULL);
    printf("Heterogenous finish! Time: %.6f\n", time1 = TIME(st, ed));
    printf("Gflops: %.2f\n", (double)M * N * L * nnz * 2 / 1024 / 1024 / 1024 / time1);
    
    gettimeofday(&st, NULL);
    struct_sptrsv_sequential_forward(a, x2, b, dim, nnz, p);
    gettimeofday(&ed, NULL);
    printf("Sequential   finish! Time: %.6f\n", time2 = TIME(st, ed));
    printf("Speed up:    %.2f\n", time2 / time1);

    int check_res;
    if ((check_res = check(x1, x2, M * N * L)) == -1)
    {
        printf("Correct!\n");
    }
    else
    {
        int wi, wj, wk, wid;
        wi = check_res / N / L;
        wj = (check_res - wi * N * L) / L;
        wk = check_res - wi * N * L - wj * L;
        printf("Incorrect at %d!\n", check_res);
        printf("index: i %d j %d k %d!\n", wi, wj, wk);
        printf("heterogeneuos: %.6f\n", x1[check_res]);
        printf("sequential:    %.6f\n", x2[check_res]);
    }
}
