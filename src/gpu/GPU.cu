#include <algorithm>
#include <stdint.h>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <mutex>
#include <vector>
#include <sys/time.h>
#include "GPU.h"


#define BAND_SIZE 256 
#define LOG_BLOCK_SIZE 8
#define LOG_NUM_BLOCKS 8
#define MAX_MARKER_OCC 8
#define HASH_PER_BLOCK MAX_MARKERS_PER_READ
#define MAX_MARKER_OCC_ELEM (1+MAX_MARKER_OCC)/2
#define BYTES_PER_HASH 4*(1+MAX_MARKER_OCC_ELEM)

#define BLOCK_SIZE (1 << LOG_BLOCK_SIZE)
#define NUM_BLOCKS (1 << LOG_NUM_BLOCKS)
#define INVALID_ID 0xffffffff

std::mutex* gpu_lock;

uint32_t** h_reads_pinned;

uint32_t** d_reads;
uint64_t** d_read_pairs;
uint32_t** d_alignments;

uint32_t** d_index_table;
uint32_t** d_hash;
uint32_t** d_marker_h;
uint32_t** d_tb_mem;
uint32_t** d_num_traceback;
uint32_t** d_common_markers;
uint32_t** d_num_common_markers;
uint32_t** d_num_hash_values;

__global__
void find_common_markers (int n, uint32_t* d_num_hash_values, uint32_t* d_reads, uint64_t* d_read_pairs, uint32_t* d_index_table, uint32_t* d_hash, uint32_t* d_marker_h, uint32_t* d_common_markers, uint32_t* d_num_common_markers) {
    int tx = threadIdx.x;
    int bs = blockDim.x;
    int bx = blockIdx.x;
    int gs = gridDim.x;

    //start address of hash table
    uint32_t hs = bx*HASH_PER_BLOCK*(BYTES_PER_HASH/4);
    uint32_t prev_rid1 = INVALID_ID;

    __shared__ uint32_t rid1;
    __shared__ uint32_t s1, e1, s2, e2;
    __shared__ uint32_t h1_arr[BLOCK_SIZE][MAX_MARKER_OCC];
    __shared__ uint32_t prefix[1+BLOCK_SIZE];
    __shared__ uint32_t curr_offset[BLOCK_SIZE];
    __shared__ int has_failed;

    int chunk_size = (n/gs);
    chunk_size += 1;

    int si = bx*chunk_size;
    int ei = si + chunk_size;
    if (ei > n) {
        ei = n;
    }
    
    for (int i = tx; i < HASH_PER_BLOCK; i+=bs) {
        d_hash[hs+i*(1+MAX_MARKER_OCC_ELEM)] = INVALID_ID;
    }

    __syncthreads();

    for (int i = si; i < ei; i++) {
        __syncthreads();
        if (tx == 0) {
            uint64_t v1 = d_read_pairs[2*i];
            uint64_t v2 = d_read_pairs[2*i+1];
            
            s1 = (v1 >> 32);
            e1 = s1 + ((v1 << 32) >> 32);

            s2 = (v2 >> 32);
            e2 = s2 + ((v2 << 32) >> 32);
            
            rid1 = s1;

            prefix[0] = bx*HASH_PER_BLOCK;

            has_failed = 0;
            if ((s1 == e1) || (s2 == e2)) {
                has_failed = 1;
            }
        }

        curr_offset[tx] = 0;
        prefix[tx+1] = 0;

        __syncthreads();

        if (has_failed > 0) {
            if (tx == 0) {
                d_num_common_markers[i] = 0;
            }
            continue;
        }

        if (rid1 != prev_rid1) {
            prev_rid1 = rid1;
            for (uint32_t p = 0; p < e1-s1; p += bs) {
                if (s1+p+tx < e1) {
                    uint32_t idx;
                    uint32_t k, h1;

                    k = d_reads[s1+p+tx];
                    
                    h1 = (~k) + (k << 11);
                    h1 = h1 + (k >> 12);
                    h1 = (h1 + (h1 << 3)) + (h1 << 8);

                    idx = h1 % bs;

                    atomicAdd(&prefix[idx+1], 1);
                }
            }

            __syncthreads();

            if (tx == 0) {
                for (int r = 0; r < BLOCK_SIZE; r++) {
                    prefix[1+r] += prefix[r];
                }
            }

            curr_offset[tx] = 0; 
            __syncthreads();

            for (int p = 0; p < e1-s1; p += bs) {
                if (s1+p+tx < e1) {
                    uint32_t idx;
                    uint32_t k, h1;

                    k = d_reads[s1+p+tx];
                    
                    h1 = (~k) + (k << 11);
                    h1 = h1 + (k >> 12);
                    h1 = (h1 + (h1 << 3)) + (h1 << 8);
                    
                    idx = h1 % bs;

                    uint32_t hst = prefix[idx];
                    uint16_t off = atomicAdd(&curr_offset[idx], 1);

                    uint32_t val = h1 % HASH_PER_BLOCK;
                    val = (val << 16) + (p+tx);

                    d_marker_h[hst+off] = val;
                }
            }

            __syncthreads();

            uint32_t mhs = prefix[tx];
            uint32_t mhe = prefix[tx+1];

            for (uint32_t m = mhs; m < mhe; m++) {
                uint32_t v = d_marker_h[m];
                uint32_t h = (v >> 16);
                uint32_t p = ((v << 16) >> 16);

                uint32_t hidx = hs + h*(1+MAX_MARKER_OCC_ELEM);
                uint32_t rid = d_hash[hidx];

                if (rid != rid1) {
                    d_hash[hidx] = rid1;
                    d_hash[hidx+1] = p+1;
                }
                else {
                    for (int q = 1; q <= MAX_MARKER_OCC_ELEM; q++) {
                        uint32_t val = d_hash[hidx+q];
                        uint32_t l, u;
                        l = ((val << 16) >> 16);
                        u = (val >> 16);
                        if (l == 0) {
                            d_hash[hidx+q] = (p+1);
                            break;
                        }
                        else if (u == 0) {
                            d_hash[hidx+q] = val + ((p+1) << 16);
                            if (q < MAX_MARKER_OCC_ELEM) {
                                d_hash[hidx+q+1] = 0; 
                            }
                            break;
                        }
                    }
                }
            }
        }
        
        __syncthreads();
        
        if (has_failed > 0) {
            if (tx == 0) {
                d_num_common_markers[i] = 0;
            }
            continue;
        }

        if (tx == 0) {
            for (int r = 0; r < bs-1; r++) {
                curr_offset[1+r] += curr_offset[r];
            }
            prefix[0] = i*HASH_PER_BLOCK;
        }

        prefix[1+tx] = 0;
        __syncthreads();

        int n1 = 0;

        for (uint32_t p = 0; p < e2-s2; p += bs) {
            uint32_t h1, rid;
            uint32_t hidx;
            uint32_t l, u;
            uint32_t k;

            rid = INVALID_ID;

            if (s2+p+tx < e2) {
                k = d_reads[s2+p+tx];
                
                h1 = (~k) + (k << 11);
                h1 = h1 + (k >> 12);
                h1 = (h1 + (h1 << 3)) + (h1 << 8);


                h1 = h1 % HASH_PER_BLOCK;

                hidx = hs + (1+MAX_MARKER_OCC_ELEM)*h1;

                rid = d_hash[hidx];
            }

            n1 = 0;
            if (rid == rid1) {
                for (int q = 0; q < MAX_MARKER_OCC_ELEM; q++) {
                    uint32_t val = d_hash[hidx+q+1];
                    l = ((val << 16) >> 16);
                    u = (val >> 16);

                    if (l != 0) {
                        if (k == d_reads[s1+l-1]) {
                            h1_arr[tx][n1] = (p+tx+1) + (l << 16);
                            n1++;
                        }
                        if (u != 0) {
                            if (k == d_reads[s1+u-1]) {
                                h1_arr[tx][n1] = (p+tx+1) + (u << 16);
                                n1++;
                            }
                        }
                        else {
                            break;
                        }
                    }
                    else {
                        break;
                    }
                }
            }

            prefix[1+tx] = n1;
            __syncthreads();

            if (tx == 0) {
                for (int d = 0; d < BLOCK_SIZE; d++) {
                    prefix[1+d] += prefix[d];
                }
            }
            __syncthreads();

            uint32_t addr_s = prefix[tx];
            uint32_t addr_e = prefix[1+tx];

            for (uint32_t addr = addr_s; addr < addr_e; addr++) {
                if (addr-i*HASH_PER_BLOCK < HASH_PER_BLOCK) {
                    d_common_markers[addr] = h1_arr[tx][addr-addr_s];
                }
            }

            __syncthreads();
            
            if (tx == 0) {
                uint32_t addr = prefix[BLOCK_SIZE];
                prefix[0] = addr;
                if (addr-i*HASH_PER_BLOCK < HASH_PER_BLOCK) {
                    d_common_markers[addr] = 0;
                }
            }
            __syncthreads();
        }
        
        __syncthreads();

        int num_common_markers =  prefix[BLOCK_SIZE] - i*HASH_PER_BLOCK;

        if (tx == 0) {
            // TODO: fail if >= HASH_PER_BLOCK?
            d_num_common_markers[i] = (num_common_markers < HASH_PER_BLOCK) ? num_common_markers : HASH_PER_BLOCK;
        }

        __syncthreads();
    }
}

__global__
void find_traceback (int n, size_t maxSkip, uint32_t* d_marker_h, uint32_t* d_common_markers, uint32_t* d_num_common_markers, uint32_t* d_tb_mem, uint32_t* d_alignments, uint32_t* d_num_traceback) {
    int tx = threadIdx.x;
    int bs = blockDim.x;
    int bx = blockIdx.x;
    int gs = gridDim.x;

    __shared__ uint32_t score[BAND_SIZE];
    __shared__ uint32_t score_pos[BAND_SIZE];
    __shared__ int num_common_markers;
    __shared__ bool stop_shared;

    for (int i = bx; i < n; i += gs) {
        uint32_t max_score = 0, max_score_pos = 0;
        uint32_t addr1 = i*MAX_MARKERS_PER_READ;
        uint32_t addr2 = bx*HASH_PER_BLOCK;

//        uint32_t start_addr = 0;
//        if (i > 0) {
//            start_addr = d_num_traceback[i-1];
//        }
//        uint32_t num_tb = d_num_traceback[i] - start_addr;
//        
        if (tx == 0) {
            num_common_markers = d_num_common_markers[i];
//            d_alignments[start_addr] = 0;
            d_alignments[addr1] = 0;
        }
        score[tx] = 0;

        __syncthreads();

        for (int p = 0; p < num_common_markers; p++) {
            uint32_t v = d_common_markers[addr1+p];
            uint32_t l = ((v << 16) >> 16);
            uint32_t u = (v >> 16);

            int ptr = p - tx - 1;

            score[tx] = 1;
            score_pos[tx] = p;

            bool stop = false;
            __syncthreads();

            while (!stop) {
                uint32_t l1, u1;
                if (ptr >= 0) {
                    uint32_t v1 = d_common_markers[addr1+ptr];
                    l1 = ((v1 << 16) >> 16);
                    u1 = (v1 >> 16);
                    if ((l1 < l) && (u1 < u) && (u-u1 < 8) && (l-l1 < 8)) {
                        uint32_t pscore = d_marker_h[addr2+ptr];
                        if (score[tx] < pscore+1) { 
                            score[tx] = pscore+1;
                            score_pos[tx] = ptr;
                        }
                    }
                }
                ptr -= bs;
                if (tx == 0) {
                    if ((ptr < 0) || (l-l1 < 8))  {
                        stop_shared = true;
                    }
                    else {
                        stop_shared = false;
                    }
                }
                __syncthreads();
                stop = stop_shared;
            }

            __syncthreads();

            // parallel reduction (max)
            for(unsigned int s = 1; s < bs; s *= 2) {
                if (tx % (2*s) == 0) {
                    if (score[tx] < score[tx+s]) { 
                        score[tx] = score[tx + s];
                        score_pos[tx] = score_pos[tx + s];
                    }
                }
                __syncthreads();
            }
            
            if (tx == 0) {
                d_marker_h[addr2+p] = score[0];
                d_tb_mem[addr2+p] = score_pos[0];
                if (score[0] > max_score) {
                    max_score = score[0];
                    max_score_pos = score_pos[0];
                }
            }
            __syncthreads();
        }

        __syncthreads();

        if (tx == 0) {
            int num_ptr = 0;

            if (max_score > 0) {
                int curr_pos = max_score_pos;
                int prev_pos = max_score_pos + 1;

                while ((curr_pos >= 0) && (prev_pos > curr_pos)) {
                    prev_pos = curr_pos;
                    d_alignments[addr1+num_ptr] = d_common_markers[addr1+curr_pos];
                    num_ptr++;
                    curr_pos = d_tb_mem[addr2+curr_pos];
                }
            }

            if (num_ptr < HASH_PER_BLOCK) {
                d_alignments[addr1+num_ptr] = 0;
            }
//            d_num_traceback[i] = num_ptr;
        }
        __syncthreads();
    }
}

extern "C" int initializeProcessors () {
    int nDevices;

    cudaGetDeviceCount(&nDevices);
    //    for (int i = 0; i < nDevices; i++) {
    //        cudaDeviceProp prop;
    //        cudaGetDeviceProperties(&prop, i);
    //        printf("Device Number: %d\n", i);
    //        printf("  Device name: %s\n", prop.name);
    //        printf("  Memory Clock Rate (KHz): %d\n",
    //                prop.memoryClockRate);
    //        printf("  Memory Bus Width (bits): %d\n",
    //                prop.memoryBusWidth);
    //        printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
    //                2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    //    }

    gpu_lock = new std::mutex[nDevices];
    
    h_reads_pinned = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));

    d_reads = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_read_pairs = (uint64_t**) malloc(nDevices*sizeof(uint64_t*));
    d_alignments = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));

    d_index_table = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_hash = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_marker_h = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_tb_mem = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_num_traceback = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_common_markers = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_num_common_markers = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_num_hash_values = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));

    cudaError_t err;
    size_t num_bytes;

    for (int k=0; k<nDevices; k++) {
        gpu_lock[k].lock();

        err = cudaSetDevice(k);
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: could not set device %d!\n", k);
            exit(1);
        }
        
        num_bytes = 2*GPU_BATCH_SIZE*MAX_MARKERS_PER_READ*sizeof(uint32_t);
        err = cudaMallocHost(&h_reads_pinned[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "ERROR: cudaMallocHost failed!\n");
            exit(1);
        }

        num_bytes = 2*GPU_BATCH_SIZE*MAX_MARKERS_PER_READ*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_reads[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        num_bytes = 2*GPU_BATCH_SIZE*sizeof(uint64_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_read_pairs[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        num_bytes = GPU_BATCH_SIZE*MAX_MARKERS_PER_READ*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_alignments[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }

        num_bytes = NUM_BLOCKS*MAX_MARKERS_PER_READ*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_index_table[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        num_bytes = NUM_BLOCKS*HASH_PER_BLOCK*BYTES_PER_HASH;
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_hash[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        num_bytes = NUM_BLOCKS*HASH_PER_BLOCK*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_marker_h[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        num_bytes = NUM_BLOCKS*HASH_PER_BLOCK*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_tb_mem[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        num_bytes = GPU_BATCH_SIZE*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_num_traceback[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
//        num_bytes = GPU_BATCH_SIZE*sizeof(uint32_t);
//        if (k==0)
//            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
//        err = cudaMalloc(&d_prefix_num_traceback[k], num_bytes); 
//        if (err != cudaSuccess) {
//            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
//            exit(1);
//        }
        
        num_bytes = GPU_BATCH_SIZE*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_num_common_markers[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        num_bytes = GPU_BATCH_SIZE*HASH_PER_BLOCK*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_common_markers[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        num_bytes = sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_num_hash_values[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        err = cudaMemset(d_num_hash_values[k], 0x0, num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMemset failed!\n");
            exit(1);
        }
        
        gpu_lock[k].unlock();
    }

    return nDevices;
}

void alignBatchGPU (size_t n, size_t deviceId, size_t maxSkip, size_t num_pos, uint32_t* h_reads, uint64_t* h_read_pairs, uint32_t* h_alignments) {
    cudaError_t err;

    bool report_time = false;

    size_t k = deviceId;

    struct timeval t1, t2, t3;
    long useconds, seconds, mseconds;

    std::memcpy(h_reads_pinned[k], h_reads, num_pos*sizeof(uint32_t));

    uint32_t* h_prefix_num_tb = (uint32_t*) malloc(n*sizeof(uint32_t));

    gpu_lock[k].lock();
    
    gettimeofday(&t1, NULL);

    err = cudaSetDevice(k);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: could not set device %zu!\n", k);
        exit(1);
    }

    err = cudaMemcpy(d_reads[k], h_reads_pinned[k], num_pos*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    
    err = cudaMemcpy(d_read_pairs[k], h_read_pairs, 2*n*sizeof(uint64_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }

    find_common_markers <<<NUM_BLOCKS, BLOCK_SIZE>>> (n, d_num_hash_values[k], d_reads[k], d_read_pairs[k], d_index_table[k], d_hash[k], d_marker_h[k], d_common_markers[k], d_num_common_markers[k]);

    err = cudaMemcpy(h_prefix_num_tb, d_num_common_markers[k], n*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }

    thrust::inclusive_scan(h_prefix_num_tb, h_prefix_num_tb+n, h_prefix_num_tb);
    
    err = cudaMemcpy(d_num_traceback[k], h_prefix_num_tb, n*sizeof(uint32_t), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }

    err = cudaMemcpy(h_prefix_num_tb, d_common_markers[k], 320*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: cudaMemcpy failed!\n");
        exit(1);
    }
    
//    for (uint32_t i=0; i<320; i++) {
//        uint32_t v, l, u;
//        v = h_prefix_num_tb[i];
//        u = (v >> 16);
//        l = ((v << 16) >> 16);
//        fprintf(stderr, "%u: %u %u\n", i, u, l);
////        fprintf(stderr, "%u: %u\n", i, v);
//    }

    find_traceback <<<NUM_BLOCKS, BAND_SIZE>>> (n, maxSkip, d_marker_h[k], d_common_markers[k], d_num_common_markers[k], d_tb_mem[k], d_alignments[k], d_num_traceback[k]);

    gettimeofday(&t2, NULL);

    uint32_t total_num_tb = 0;
    if (n > 0) {
        total_num_tb = h_prefix_num_tb[n-1];
    }

    uint32_t* h_tb;
    err =  cudaMallocHost ((uint32_t**) &h_tb, total_num_tb*sizeof(uint32_t));

    //    err = cudaMemcpy(h_tb, d_alignments[k], total_num_tb*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    //    if (err != cudaSuccess) {
    //        fprintf(stderr, "GPU_ERROR: cudaMemcpy failed!\n");
    //        exit(1);
    //    }

    gpu_lock[k].unlock();

    gettimeofday(&t3, NULL);

    for (int i=0; i<n; i++) {
        uint32_t start;
        //uint32_t end;
        uint32_t num_tb;
        start = 0;
        if (i > 0) {
            start = h_prefix_num_tb[i-1];
        }
        //end = h_prefix_num_tb[i];
        //TODO: remove after fixind find_traceback
        num_tb = 0; //end - start;
        std::memcpy(&h_alignments[i*MAX_MARKERS_PER_READ], &h_tb[start], num_tb*sizeof(uint32_t));
        if (num_tb < MAX_MARKERS_PER_READ) {
            h_alignments[i*MAX_MARKERS_PER_READ+num_tb] = 0; 
        }
    }

    if (report_time) {
        useconds = t2.tv_usec - t1.tv_usec;
        seconds = t2.tv_sec - t1.tv_sec;
        mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;

        fprintf(stderr, "Time elapsed (t2-t1): %ld msec \n", mseconds);

        useconds = t3.tv_usec - t1.tv_usec;
        seconds = t3.tv_sec - t1.tv_sec;
        mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;

        fprintf(stderr, "Time elapsed (t3-t1): %ld msec \n", mseconds);
        
        fprintf(stderr, "Num bytes transferred (reads): %zu \n", num_pos*4);
        fprintf(stderr, "Num bytes transferred (traceback): %u \n", total_num_tb*4);
    }

    cudaFree(h_tb);
    free(h_prefix_num_tb);

    return;
}

extern "C" void shutdownProcessors(int nDevices) {
    for (int k=0; k<nDevices; k++) {
        cudaFree(h_reads_pinned[k]);

        cudaFree(d_reads[k]);
        cudaFree(d_read_pairs[k]);
        cudaFree(d_alignments[k]);

        cudaFree(d_index_table[k]);
        cudaFree(d_hash[k]);
        cudaFree(d_marker_h[k]);
        cudaFree(d_tb_mem[k]);
        cudaFree(d_num_traceback[k]);
        cudaFree(d_common_markers[k]);
        cudaFree(d_num_common_markers[k]);
        cudaFree(d_num_hash_values[k]);
    }
        
    free(h_reads_pinned);
    
    free(d_reads);
    free(d_read_pairs);
    free(d_alignments);

    free(d_index_table);
    free(d_hash);
    free(d_marker_h);
    free(d_tb_mem);
    free(d_num_traceback);
    free(d_common_markers);
    free(d_num_common_markers);
    free(d_num_hash_values);

    delete(gpu_lock);

}
