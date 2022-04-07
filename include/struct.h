#ifndef STRUCT_SPTRSV_
#define STRUCT_SPTRSV_

#define ROW_P 8
#define COL_P 8

#define INDEX3(i, j, k, N, L) ((i)*(N)*(L) + (j)*(L) + (k))
#define VALID(i, j, k, M, N, L) ((i >= 0) && (i < (M)) && (j >= 0) && (j < (N)) && (k >= 0) && (k < (L)))

#define CAL_ID(row, col) (((row) >= 0 && (col) >= 0 && (row) < ROW_P && (col) < COL_P)? ((row) * ROW_P + (col)) : -1)

#define gridmax(a, b) (((a)>(b))?(a):(b))
#define gridmin(a, b) (((a)<(b))?(a):(b))

#define syc_rcol(sid, syc) (syc - 1 - sid)

#define simd_set_a_nnz(src) \
simd_set_doublev8(*(src), *(src + nnz), *(src + nnz2), *(src + nnz3), \
*(src + nnz4), *(src + nnz5), *(src + nnz6), *(src + nnz7));

typedef struct
{
    const void *a;
    void *x;
    const void *b;
    int M;
    int N;
    int L;
    int nnz;
    int *p;
    int *col_id;
    int *z_off;
    int halo;
    int cols;
} struct_info;

#endif
