#include <stdio.h>
#include <slave.h>
#include <simd.h>
#include <crts.h>
#include "../include/struct.h"
#include "../include/comm.h"

#define PERF

void struct_sptrsv_2skew_f(const TYPE *a, TYPE *x, const TYPE *b, int M, int N, int L, \
int nnz, int halo, int cols, int *p, int *col_m, int *off_m)
{

    #ifdef PERF
        unsigned long time1, time2;
    #endif

    int K = 4096 / sizeof(TYPE);     //vertical block size
    int gra = 64 / sizeof(TYPE);    //granuity of communication, should be a factor of K ***
    int vec = 64 / sizeof(TYPE);
    int nnz2 = nnz * 2, nnz3 = nnz * 3, nnz4 = nnz * 4, \
    nnz5 = nnz * 5, nnz6 = nnz * 6, nnz7 = nnz * 7;

    int oneLEN = 0;
    if (K >= L) oneLEN = 1;

    int overlap = gridmin(K, L);

    //cross partition
    int id, row, col;
    id = CRTS_smng_get_tid();
    row = (CRTS_smng_get_rid()) % ROW_P;
    col = (CRTS_smng_get_cid()) % COL_P;

    //debug
    int pflag = 0;
    if (id == 0) pflag = 1;
    
    int halo2 = 2 * halo;
    int row_num = 8, col_num = 8;

    //DMA reply
    crts_rply_t rply = 0;
    int dmas = 0;

    //buffer
    __attribute((aligned(256))) TYPE x_buf[cols][K + halo2];
    __attribute((aligned(256))) TYPE res_buf[2][K + halo2];
    __attribute((aligned(256))) TYPE a_buf[2][K * nnz];

    int res_id = 0;
    int a_id = 0;

    //column id and vertical offset
    int col_t[nnz];
    int off_t[nnz];
    int pos[2*cols];
    CRTS_dma_iget(pos, p, 2 * cols * sizeof(int), &rply);
    CRTS_dma_iget(col_t, col_m, nnz * sizeof(int), &rply);
    CRTS_dma_iget(off_t, off_m, nnz * sizeof(int), &rply);
    dmas += 3;
    CRTS_dma_wait_value(&rply, dmas);

    //buffer communication
    volatile int done_k[cols]; //neighbour_done_signal
    volatile int recv_k[cols]; //neighbour_receive_signal
    int local_recv_k[cols]; //to avoid sychronization problem
    int local_done_k[cols]; //to avoid sychronization problem

    int done_offset[cols];
    int recv_offset[cols];

    //dma communication
    volatile int dma_done[cols];
    int dma_id[cols];
    int dma_check[cols];
    int dma_dii[cols], dma_djj[cols];

    //whether communication is needed
    int fneed[cols];
    int bneed[cols];

    //communication id
    int bid[cols], fid[cols];

    int i, j, k;
    int u, v;

    int no_cols_nnz = 0;
    u = 0;
    while (u < nnz && col_t[u] != cols) u++;
    no_cols_nnz = u;
    
    int count = -1;
    for (u = 0; u < cols; u++)
    {
        recv_k[u] = done_k[u] = -1;
        dma_done[u] = -1;
        local_recv_k[u] = 0;
        local_done_k[u] = 0;
    }

    // calculate id
    for (u = 0; u < cols; u++)
    {
        int di = pos[u*2], dj = pos[u*2+1];
        int ti = row + di, tj = col + dj + 2 * di;
        if (oneLEN && ti < 0)
        {
            dma_check[u] = 1;
        }
        else
            dma_check[u] = 0;
        if (oneLEN && tj < 0) 
        {
            tj = tj + col_num;
            done_offset[u] = -L;
        }
        else
            done_offset[u] = 0;
        fid[u] = CAL_ID(ti, tj);
        if (!oneLEN && fid[u] == -1) dma_check[u] = 1;
        if (dma_check[u])
        {
            if (ti < 0) dma_dii[u] = -1; else dma_dii[u] = 0;
            if (tj < 0) dma_djj[u] = -1; else dma_djj[u] = 0;
        }

        ti = row - di, tj = col - dj - 2 * di;
        if (oneLEN && ti >= row_num)
        {
            dma_id[u] = CAL_ID(ti - row_num, (tj + col_num) % col_num);
        }
        else
            dma_id[u] = -1;
        if (oneLEN && tj >= col_num) 
        {
            tj = tj - col_num;
            done_k[u] = -L - 1;
        }
        recv_offset[u] = 0;
        bid[u] = CAL_ID(ti, tj);
        if (!oneLEN && bid[u] == -1) dma_id[u] = CAL_ID((ti + row_num) % row_num, (tj + col_num) % col_num);

        //if (id == 7) printf("id %d col_id %d di %d dj %d ti %d tj %d bid %d\n", id, u, di, dj, ti, tj, bid[u]);
    }

    CRTS_ssync_array();

    #ifdef PERF
        time1 = CRTS_stime_cycle();
    #endif

    //M/N_iter iterations have to be done anyway
    int M_iter, N_iter, L_iter;
    M_iter = (M + row_num - 1) / row_num;
    L_iter = (L + K - 1) / K;

    //first column
    int new_col_j = col - 2 * row;
    while (new_col_j < 0) new_col_j += col_num;
    //if (id == 8) printf("new_col_j %d\n", new_col_j);
    if (VALID(row, new_col_j, 0, M, N, 1))
    {
        int nextlen = gridmin(L, K);
        dmas += 2;
        CRTS_dma_iget(a_buf[a_id], a + INDEX3(row, new_col_j, 0, N, L) * nnz, nextlen * nnz * sizeof(TYPE), &rply);
        CRTS_dma_iget(res_buf[res_id] + halo, b + INDEX3(row, new_col_j, 0, N, L), nextlen * sizeof(TYPE), &rply);
    }

    //main loop
    for (int ii = 0; ii < M_iter; ii++)
    {
        int width = gridmin(M - ii * row_num, row_num);
        int right_off;
        right_off = 8 - 2 * (width - 1);
        N_iter = (N - right_off + col_num - 1) / col_num + 1;
    for (int jj = 0; jj < N_iter; jj++)
    {
        i = ii * row_num + row;
        j = jj * col_num + col - 2 * row;

        if (!oneLEN) CRTS_ssync_array();

        for (u = 0; u < cols; u++)
        for (v = 0; v < halo; v++)
            x_buf[u][v] = 0;

        for (v = 0; v < halo; v++)
        {
            res_buf[res_id][v] = 0;
        }
        
    if (VALID(i, j, 0, M, N, 1))
    {
        //*****************main part*****************//
        
        //fneed = need to recv from, bneed = need to send to
        for (u = 0; u < cols; u++)
        {
            int di = pos[u*2], dj = pos[u*2+1];
            int ti = i + di, tj = j + dj;
            if (VALID(ti, tj, 0, M, N, 1)) fneed[u] = 1; else fneed[u] = 0;
            if (fid[u] == -1 || !fneed[u]) local_recv_k[u] += L;

            ti = i - di, tj = j - dj;
            if (VALID(ti, tj, 0, M, N, 1)) bneed[u] = 1; else bneed[u] = 0;
            if (bid[u] != -1 && !bneed[u]) recv_offset[u] -= L;
        }

        int count_max = count + L;

        //vertical block size K
        for (k = 0; k < L; k+=K)
        {
            int len = gridmin(L - k, K);

            //DMA other columns
            int len_halo = gridmin(len + halo, L - k);

            for (u = 0; u < cols; u++)
            {
                if (fid[u] == -1 && fneed[u] == 1)
                {
                    int di = pos[u*2], dj = pos[u*2+1];
                    int ti = i + di, tj = j + dj;
                    if (dma_check[u])
                    {
                        int halo_block_id = ti * N * L + tj * L + k/K;
                        while (dma_done[u] < halo_block_id);
                    }
                    CRTS_dma_iget(x_buf[u] + halo, x + INDEX3(ti, tj, k, N, L), len_halo * sizeof(TYPE), &rply);
                    dmas++;
                }
                else if (fneed[u] == 0)
                {
                    for (v = 0; v < K + halo; v++)
                        x_buf[u][v + halo] = 0;
                }
            }
            CRTS_dma_wait_value(&rply, dmas);

            if (k + len >= L)
            {
                //ti, tj indicate next column position
                int ti = i, tj = j + col_num;
                if (tj >= N)
                {
                    ti = ti + row_num;
                    tj = new_col_j;
                }
                if (VALID(ti, tj, 0, M, N, 1))
                {
                    int nextlen = gridmin(K, L);
                    dmas += 2;
                    CRTS_dma_iget(a_buf[(a_id+1)%2], a + INDEX3(ti, tj, 0, N, L) * nnz, nextlen * nnz * sizeof(TYPE), &rply);
                    CRTS_dma_iget(res_buf[(res_id+1)%2] + halo, b + INDEX3(ti, tj, 0, N, L), nextlen * sizeof(TYPE), &rply);
                }
            }
            else
            {
                int nextlen = gridmin(L - k - len, K);
                dmas += 2;
                CRTS_dma_iget(a_buf[(a_id+1)%2], a + INDEX3(i, j, k + len, N, L) * nnz, nextlen * nnz * sizeof(TYPE), &rply);
                CRTS_dma_iget(res_buf[(res_id+1)%2] + halo, b + INDEX3(i, j, k + len, N, L), nextlen * sizeof(TYPE), &rply);
            }

            int c;

            //calculation loop, vectorized
            int commu_index = halo;

            for (c = 0; c < len; c+=vec)
            {
                count += vec;
                if (c + vec > len) count -= c + vec - len;

                //wait until data is ready
                int flag = 1;

                // if (id == 6) printf("%d bid %d bneed %d\n", id, bid[2], bneed[2]);
                // if (id == 7) printf("%d fid %d fneed %d\n", id, fid[2], fneed[2]);

                while(flag)
                {
                    flag = 0;
                    for (u = 0; u < cols; u++)
                    {
                        if (recv_k[u] + local_recv_k[u] < gridmin(count + halo, count_max))
                            flag = 1;
                        //if (id == 7) printf("%d recv_k[%d] %d\n", id, u, recv_k[u]);
                    }
                }

                //if (i == 5 && j == 14) printf("id %d c %d count %d done_k %d local_done %d\n", id, c, count, done_k[0] + local_done_k[0], local_done_k[0]);
                
                TYPE_V res, atmp, xtmp;
                simd_loadu(res, res_buf[res_id] + c + halo);
                for (u = 0; u < no_cols_nnz; u++)
                {
                    atmp = simd_set_a_nnz(a_buf[a_id] + c * nnz + u);
                    simd_loadu(xtmp, x_buf[col_t[u]] + c + off_t[u] + halo);
                    res = res - atmp * xtmp;
                }
                simd_storeu(res, res_buf[res_id] + c + halo);

                //if (id == 0 && i == 0 && j == 0 && c == 0)
                //    printf("halo %d\n", halo);
                //for (u = 0; u < cols; u++)
                //    printf("id %d i %d j %d c %d u %d recv_k %d done_k %d\n", id, i, j, c, u, recv_k[u] + local_recv_k[u], done_k[u] + local_done_k[u]);

                // printf("id %d i %d j %d recv_k %d %d %d done_k %d %d %d c %d\n", \
                // id, i, j, recv_k[0] + local_recv_k[0], recv_k[1] + local_recv_k[1], \
                // recv_k[2] + local_recv_k[2], done_k[0] + local_done_k[0], \
                // done_k[1] + local_done_k[1], done_k[2] + local_done_k[2], c);

                if (no_cols_nnz < nnz)
                {
                    for (v = 0; v < vec; v++)
                    {
                        TYPE d_xtmp = 0;
                        for (u = no_cols_nnz; u < nnz; u++)
                        {
                            res_buf[res_id][c + v + halo] -= a_buf[a_id][(c + v)* nnz + u] * res_buf[res_id][c + off_t[u] + v + halo];
                        }
                        //x_buf[cols][c + v + halo] = x_buf[cols][c + v + halo] + d_xtmp;
                    }
                }

                for (u = 0; u < cols; u++)
                {
                    if (fid[u] != -1)
                    {
                        int* remote_r = remote_ldm_addr(&done_k[u], fid[u]);
                        (*remote_r) = count + done_offset[u];
                    }
                }

                for (u = cols - 1; u >= 0; u--)
                {
                    if (bid[u] != -1 && bneed[u])
                    {
                        //if (id == 7) printf("u %d %d %d %d\n", u, count + recv_offset[u], done_k[u] + local_done_k[u], overlap - halo - vec);
                        while ((count - (done_k[u] + local_done_k[u])) \
                        >= overlap - halo);

                        rstore8(res_buf[res_id] + commu_index, x_buf[u] + commu_index, bid[u]);

                        int* remote_r = remote_ldm_addr(&recv_k[u], bid[u]);
                        MEM_FENCE;
                        (*remote_r) = count + recv_offset[u];
                    }
                }
                commu_index += gra;

            }

            //dma wb
            CRTS_dma_put(x + INDEX3(i, j, k, N, L), res_buf[res_id] + halo, len * sizeof(TYPE));

            //inform halo dma threads
            for (u = 0; u < cols; u++)
            {
                if (dma_id[u] != -1)
                {
                    int *remote_id = remote_ldm_addr(&dma_done[u], dma_id[u]);
                    (*remote_id) = i * N * L + j * L + k;
                }
            }

            //post
            for (u = 0; u < cols; u++)
            {
                for (v = 0; v < halo; v++)
                {
                    x_buf[u][v] = x_buf[u][v + K];
                }
            }

            for (v = 0; v < halo; v++)
                res_buf[(res_id+1)%2][v] = res_buf[res_id][v + K];

            a_id = (a_id + 1) % 2;
            res_id = (res_id + 1) % 2;
        }
    }
    else
    {
        //if (id == 8) printf("id %d i %d j %d\n", id, i, j);
        count += L;
        for (u = 0; u < cols; u++)
        {
            local_recv_k[u] += L;
            recv_offset[u] -= L;

            if (fid[u] != -1)
            {
                int* remote_r = remote_ldm_addr(&done_k[u], fid[u]);
                (*remote_r) = count + done_offset[u];
            }
        }
    }

    }
    }

    CRTS_ssync_array();

    #ifdef PERF
        time2 = CRTS_stime_cycle();
        if (id == 0)
        {
            double Gflops = 1.0 * M * N * L * 2 * nnz / 1024 / 1024 / 1024;
            double Gmems  = 1.0 * M * N * L * (nnz + 2) * sizeof(TYPE) / 1024 / 1024 / 1024;
            double total_time = (time2 - time1) / 2.25 / 1000 / 1000 / 1000;
            printf("id: %d total time:       %.6f\n", id, total_time);
            printf("id: %d Gflops:           %.6f\n", id, Gflops / total_time);
            printf("id: %d Memory bandwidth: %.6f\n", id, Gmems / total_time);
        }
    #endif

}

void struct_sptrsv_2skew_entry(void *_info)
{
    struct_info info;
    CRTS_dma_get(&info, _info, sizeof(struct_info));

    struct_sptrsv_2skew_f(info.a, info.x, info.b, info.M, info.N, info.L, \
    info.nnz, info.halo, info.cols, info.p, info.col_id, info.z_off);
}