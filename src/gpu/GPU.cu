#include <algorithm>
#include <stdint.h>
#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <mutex>
#include <vector>
#include <sys/time.h>
#include "GPU.h"


#define BAND_SIZE 32 
#define LOG_BLOCK_SIZE 7
#define LOG_NUM_BLOCKS 10

#define BLOCK_SIZE (1 << LOG_BLOCK_SIZE)
#define NUM_BLOCKS (1 << LOG_NUM_BLOCKS)

std::mutex* gpu_lock;

uint32_t num_unique_markers;

uint32_t** d_alignments;
uint32_t** d_score;
uint32_t** d_score_pos;
uint32_t** d_num_traceback;
uint32_t** d_common_markers;
uint32_t** d_num_common_markers;

__global__
void find_common_markers (uint64_t maxMarkerFrequency, uint64_t n, uint32_t num_unique_markers, uint64_t* read_pairs, uint64_t* index_table, uint64_t* rid_marker_pos, uint64_t* sorted_rid_marker_pos, uint32_t* num_common_markers, uint32_t* common_markers)
{
    int tx = threadIdx.x;
    int bs = blockDim.x;
    int bx = blockIdx.x;
    int gs = gridDim.x;

    uint64_t m_mask = ((uint64_t) 1 << 32) - 1;
    uint64_t p_mask = ((uint64_t) 1 << SHASTA_LOG_MAX_MARKERS_PER_READ) - 1;
    
    __shared__ uint32_t prefix[1+BLOCK_SIZE];

    __syncthreads();

    for (int i = bx; i < n; i+=gs) {
        if (tx == 0) {
            prefix[tx] = i*SHASTA_MAX_MARKERS_PER_READ;
            num_common_markers[i] = 0;
        }
        __syncthreads();

        uint64_t v1 = read_pairs[2*i];
        uint64_t v2 = read_pairs[2*i+1];
        uint64_t rid1 = (v1 >> 32);
        uint64_t rid2 = (v2 >> 32);
        uint64_t l1 = ((v1 << 32) >> 32);
        uint64_t l2 = ((v2 << 32) >> 32);

        if ((l1 > 0) && (l2 > 0)) {
            uint64_t s2 = index_table[rid2*num_unique_markers];
            uint64_t e2 = s2+l2;

            for (uint64_t j = s2; j < e2; j += bs) {
                uint64_t idx = tx+j;
                uint64_t marker;
                uint64_t sm1=0, sm2=0, em1=0, em2=0;

                prefix[1+tx] = 0; 

                if (idx < e2) {
                    uint64_t v = rid_marker_pos[idx];
                    marker = ((v >> SHASTA_LOG_MAX_MARKERS_PER_READ) & m_mask);

                    sm1 = index_table[rid1*num_unique_markers+marker];
                    em1 = index_table[rid1*num_unique_markers+marker+1];
                    sm2 = index_table[rid2*num_unique_markers+marker];
                    em2 = index_table[rid2*num_unique_markers+marker+1];

                    if ((em1 - sm1 <= maxMarkerFrequency) && (em2 - sm2 <= maxMarkerFrequency)) {
                        prefix[1+tx] = (em1-sm1);
                    }
                }

                __syncthreads();

                if (tx == 0) {
                    for (int r = 0; r < BLOCK_SIZE; r++) {
                        prefix[1+r] += prefix[r];
                    }
                }

                __syncthreads();

                uint32_t mhs = prefix[tx];
                uint32_t mhe = prefix[1+tx];

                for (uint64_t k1 = 0; k1 < (mhe-mhs); k1++) {
                    if (mhs+k1 < (i+1)*SHASTA_MAX_MARKERS_PER_READ) {
                        uint64_t sv1 = sorted_rid_marker_pos[sm1+k1];
                        uint32_t cm = (sv1 & p_mask) + 1;
                        cm = (cm << 16) + (1+idx-s2);
                        common_markers[mhs+k1] = cm;
                    }
                }

                __syncthreads();

                if (tx == 0) {
                    prefix[tx] = prefix[BLOCK_SIZE];
                }

                __syncthreads();
            }

            if (tx == 0) {
                uint32_t num_common = prefix[tx] - i*SHASTA_MAX_MARKERS_PER_READ;
                if (num_common < SHASTA_MAX_MARKERS_PER_READ) {
                    num_common_markers[i] = num_common;
                }
                else {
                    num_common_markers[i] = 0;
                }
            }
        }

        __syncthreads();
    }
}

__global__
void find_traceback (int n, size_t maxSkip, uint32_t* d_score, uint32_t* d_common_markers, uint32_t* d_num_common_markers, uint32_t* d_score_pos, uint32_t* d_alignments, uint32_t* d_num_traceback) {
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
        uint32_t addr1 = i*SHASTA_MAX_MARKERS_PER_READ;
        uint32_t addr2 = bx*SHASTA_MAX_MARKERS_PER_READ;
        uint32_t addr3 = i*SHASTA_MAX_TB;

        if (tx == 0) {
            num_common_markers = d_num_common_markers[i];
            d_alignments[addr3] = 0;
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
                    if ((l1 < l) && (u1 < u) && (u-u1 <= maxSkip) && (l-l1 <= maxSkip)) {
                        uint32_t pscore = d_score[addr2+ptr];
                        if (score[tx] < pscore+1) { 
                            score[tx] = pscore+1;
                            score_pos[tx] = ptr;
                        }
                    }
                }
                ptr -= bs;
                if (tx == bs-1) {
                    if ((ptr < 0) || (l-l1 > maxSkip))  {
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
                d_score[addr2+p] = score[0];
                d_score_pos[addr2+p] = score_pos[0];
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
                    if (num_ptr < SHASTA_MAX_TB) {
                        d_alignments[addr3+num_ptr] = d_common_markers[addr1+curr_pos];
                    }
                    num_ptr++;
                    curr_pos = d_score_pos[addr2+curr_pos];
                }
            }

            if (num_ptr < SHASTA_MAX_TB) {
                d_alignments[addr3+num_ptr] = 0;
            }
            else {
                d_alignments[addr3] = 0;
            }
            d_num_traceback[i] = num_ptr;
        }
        __syncthreads();
    }
}

extern "C" int shasta_initializeProcessors (size_t numUniqueMarkers) {
    int nDevices;

    num_unique_markers = (uint32_t) numUniqueMarkers;

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
    
    d_alignments = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));

    d_score = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_score_pos = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_num_traceback = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_common_markers = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_num_common_markers = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));

    cudaError_t err;
    size_t num_bytes;

    for (int k=0; k<nDevices; k++) {
        gpu_lock[k].lock();

        err = cudaSetDevice(k);
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: could not set device %d!\n", k);
            exit(1);
        }
        
        num_bytes = SHASTA_GPU_BATCH_SIZE*SHASTA_MAX_TB*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_alignments[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }

        num_bytes = NUM_BLOCKS*SHASTA_MAX_MARKERS_PER_READ*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_score[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        num_bytes = NUM_BLOCKS*SHASTA_MAX_MARKERS_PER_READ*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_score_pos[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        num_bytes = SHASTA_GPU_BATCH_SIZE*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_num_traceback[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        num_bytes = SHASTA_GPU_BATCH_SIZE*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_num_common_markers[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        num_bytes = SHASTA_GPU_BATCH_SIZE*SHASTA_MAX_MARKERS_PER_READ*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_common_markers[k], num_bytes); 
        if (err != cudaSuccess) {
            fprintf(stderr, "GPU_ERROR: cudaMalloc failed!\n");
            exit(1);
        }
        
        gpu_lock[k].unlock();
    }

    return nDevices;
}

extern "C" void shasta_alignBatchGPU (size_t deviceId, size_t maxMarkerFrequency, size_t maxSkip, size_t n, uint64_t num_pos, uint64_t num_reads, uint64_t* batch_rid_marker_pos, uint64_t* batch_rid_markers, uint64_t* batch_read_pairs, uint32_t* h_alignments, uint32_t* h_num_traceback) {
    bool report_time = false;

    size_t k = deviceId;

    struct timeval t1, t2, t3;
    long useconds, seconds, mseconds;

    gpu_lock[k].lock();
    
    cudaError_t err; 

    err = cudaSetDevice(k);
    if (err != cudaSuccess) {
        fprintf(stderr, "GPU_ERROR: could not set device %zu!\n", k);
        exit(1);
    }
    
    gettimeofday(&t1, NULL);

    thrust::device_vector<uint64_t> t_d_rid_marker_pos (batch_rid_marker_pos, batch_rid_marker_pos + num_pos);
    thrust::device_vector<uint64_t> t_d_sorted_rid_marker_pos (batch_rid_marker_pos, batch_rid_marker_pos+num_pos);
    thrust::device_vector<uint64_t> t_d_rid_markers (batch_rid_markers, batch_rid_markers + num_reads*num_unique_markers+1);
    thrust::device_vector<uint64_t> t_d_read_pairs (batch_read_pairs, batch_read_pairs+2*n);
    thrust::device_vector<uint64_t> t_d_index_table (num_reads*num_unique_markers+1);

    thrust::sort(t_d_sorted_rid_marker_pos.begin(), t_d_sorted_rid_marker_pos.end());

    gettimeofday(&t2, NULL);

    thrust::lower_bound(t_d_sorted_rid_marker_pos.begin(),
            t_d_sorted_rid_marker_pos.end(),
            t_d_rid_markers.begin(),
            t_d_rid_markers.end(),
            t_d_index_table.begin());

    uint64_t* d_sorted_rid_marker_pos = thrust::raw_pointer_cast (t_d_sorted_rid_marker_pos.data());
    uint64_t* d_rid_marker_pos = thrust::raw_pointer_cast (t_d_rid_marker_pos.data());
    uint64_t* d_index_table = thrust::raw_pointer_cast (t_d_index_table.data());
    uint64_t* d_read_pairs = thrust::raw_pointer_cast (t_d_read_pairs.data());
    
    find_common_markers <<<NUM_BLOCKS, BLOCK_SIZE>>> (maxMarkerFrequency, n, num_unique_markers, d_read_pairs, d_index_table, d_rid_marker_pos, d_sorted_rid_marker_pos, d_num_common_markers[k], d_common_markers[k]);

    find_traceback <<<NUM_BLOCKS, BAND_SIZE>>>(n, maxSkip, d_score[k], d_common_markers[k], d_num_common_markers[k], d_score_pos[k], d_alignments[k], d_num_traceback[k]);


    err = cudaMemcpy(h_num_traceback, d_num_traceback[k], n*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "1. Error: cudaMemcpy failed!\n");
        exit(1);
    }

    err = cudaMemcpy(h_alignments, d_alignments[k], n*SHASTA_MAX_TB*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        fprintf(stderr, "2. Error: cudaMemcpy failed!\n");
        exit(1);
    }
    gettimeofday(&t3, NULL);
    
    if (report_time) {
        useconds = t2.tv_usec - t1.tv_usec;
        seconds = t2.tv_sec - t1.tv_sec;
        mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
        fprintf(stderr, "Time elapsed (t2-t1): %ld msec \n", mseconds);

        useconds = t3.tv_usec - t1.tv_usec;
        seconds = t3.tv_sec - t1.tv_sec;
        mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
        fprintf(stderr, "Time elapsed (t3-t1): %ld msec \n", mseconds);
    }

    gpu_lock[k].unlock();

    return;
}

extern "C" void shasta_shutdownProcessors(int nDevices) {
    for (int k=0; k<nDevices; k++) {
        cudaFree(d_alignments[k]);

        cudaFree(d_score[k]);
        cudaFree(d_score_pos[k]);
        cudaFree(d_num_traceback[k]);
        cudaFree(d_common_markers[k]);
        cudaFree(d_num_common_markers[k]);
    }
    free(d_alignments);

    free(d_score);
    free(d_score_pos);
    free(d_num_traceback);
    free(d_common_markers);
    free(d_num_common_markers);

    delete(gpu_lock);

}
