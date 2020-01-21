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
#include <condition_variable>
#include "GPU.h"

#include "../stdexcept.hpp"


#define BAND_SIZE 32 
#define LOG_BLOCK_SIZE 7

#define BLOCK_SIZE (1 << LOG_BLOCK_SIZE)

int NUM_BLOCKS;
int NUM_DEVICES;
size_t BATCH_SIZE;

std::mutex mu;
std::condition_variable cv;
std::vector<int> available_gpus;

uint32_t num_unique_markers;

uint32_t** d_alignments;
float** d_score;
uint32_t** d_score_pos;
uint32_t** d_num_traceback;
uint64_t** d_common_markers;
uint32_t** d_num_common_markers;
uint64_t** d_batch_rid_markers;
uint64_t** d_rid_marker_pos;

using namespace shasta;


__global__
void initialize_batch_rid_markers (uint64_t* batch_rid_markers, uint32_t num_unique_markers, uint32_t batch_size) {
    int tx = threadIdx.x;
    int bs = blockDim.x;
    int bx = blockIdx.x;
    int gs = gridDim.x;

    for (uint64_t i=bx; i < 2*batch_size; i+=gs) {
        uint64_t v, val;
        v = (i << (32+SHASTA_LOG_MAX_MARKERS_PER_READ));
        for (uint64_t j=tx; j< num_unique_markers; j+=bs) {
            val = v + (j << SHASTA_LOG_MAX_MARKERS_PER_READ);
            batch_rid_markers[i*num_unique_markers+j] = val;
        }
    }
    if (bx==0) {
        if (tx == 0) {
            uint64_t v = 2*batch_size;
            batch_rid_markers[2*batch_size*num_unique_markers] = (v << (32+SHASTA_LOG_MAX_MARKERS_PER_READ));
        }
    }
}

__global__
void skip_high_frequency_markers (uint64_t maxMarkerFrequency, uint32_t num_unique_markers, uint64_t* index_table, uint64_t* rid_marker_pos, uint64_t* sorted_rid_marker_pos, uint16_t* adjusted_pos) {
    int tx = threadIdx.x;
    int bs = blockDim.x;
    int bx = blockIdx.x;
    
    __shared__ uint64_t s, e;
    
    uint16_t sum = 0;
    uint16_t prev_sum = 0;
    uint64_t m_mask = ((uint64_t) 1 << 32) - 1;

    if (tx == 0) {
        s = index_table[bx*num_unique_markers];
        e = index_table[(bx+1)*num_unique_markers];
    }
    __syncthreads();

    for (uint64_t i = s; i < e; i+=bs) {
        uint64_t idx = i+tx;
        uint64_t v = rid_marker_pos[idx];
        uint64_t marker = ((v >> SHASTA_LOG_MAX_MARKERS_PER_READ) & m_mask);
        uint64_t sm=0, em=0;

        sum = prev_sum;
        if (idx < e) {
            uint64_t v = rid_marker_pos[idx];
            marker = ((v >> SHASTA_LOG_MAX_MARKERS_PER_READ) & m_mask);

            sm = index_table[bx*num_unique_markers+marker];
            em = index_table[bx*num_unique_markers+marker+1];
            if ((em-sm) <= maxMarkerFrequency) {
                sum += 1;
            }
        }

        __syncthreads();

        for (int s = 1; s <= bs; s *= 2) {
            int val = __shfl_up_sync(0xffffffff, sum, s, bs);

            if (tx >= s) {
                sum += val;
            }
        }
            
        __syncthreads();

        if (idx < e) {
            adjusted_pos[idx] = sum;
        }
        
        int prev = __shfl_down_sync(0xffffffff, sum, bs-1, bs);
        if (tx == 0) {
            prev_sum = prev;
        }

        __syncthreads();
    }
}

__global__
void find_common_markers (uint64_t maxMarkerFrequency, uint64_t n, uint32_t num_unique_markers, uint64_t* read_pairs, uint64_t* index_table, uint64_t* rid_marker_pos, uint64_t* sorted_rid_marker_pos, uint16_t* adjusted_pos, uint32_t* num_common_markers, uint64_t* common_markers)
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
            prefix[tx] = i*SHASTA_MAX_COMMON_MARKERS_PER_READ;
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
            uint64_t s1 = index_table[rid1*num_unique_markers];
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
                
                __syncthreads();

                for (uint64_t k1 = 0; k1 < (mhe-mhs); k1++) {
                    if (mhs+k1 < (i+1)*SHASTA_MAX_COMMON_MARKERS_PER_READ) {
                        uint64_t sv1 = sorted_rid_marker_pos[sm1+k1];
                        uint64_t cm = (sv1 & p_mask);
                        uint64_t adj_pos1, adj_pos2;
                        
                        adj_pos1 = adjusted_pos[s1+cm];
                        adj_pos2 = adjusted_pos[idx]; 

                        cm = ((cm+1) << 16) + (1+idx-s2);
                        common_markers[mhs+k1] = cm + (adj_pos1 << 48) + (adj_pos2 << 32);
                    }
                }


                if (tx == 0) {
                    prefix[tx] = prefix[BLOCK_SIZE];
                }

                __syncthreads();
            }

            if (tx == 0) {
                uint32_t num_common = prefix[tx] - i*SHASTA_MAX_COMMON_MARKERS_PER_READ;
                num_common_markers[i] = num_common;
            }
        }

        __syncthreads();
    }
}

__global__
void find_traceback (int n, size_t maxSkip, size_t maxDrift, float* d_score, uint64_t const * __restrict__ d_common_markers, uint32_t const*  __restrict__ d_num_common_markers, uint32_t* d_score_pos, uint32_t* d_alignments, uint32_t* d_num_traceback, bool get_complete_traceback) {
    int tx = threadIdx.x;
    int bs = blockDim.x;
    int bx = blockIdx.x;
    int gs = gridDim.x;

    float score;
    uint32_t score_pos;
    int num_common_markers;

    uint64_t p_mask = ((uint64_t) 1 << SHASTA_LOG_MAX_MARKERS_PER_READ) - 1;

    for (int i = bx; i < n; i += gs) {
        float max_score = 0;
        uint32_t max_score_pos = 0;
        uint32_t addr1 = i*SHASTA_MAX_COMMON_MARKERS_PER_READ;
        uint32_t addr2 = bx*SHASTA_MAX_COMMON_MARKERS_PER_READ;
        uint32_t addr3 = i*SHASTA_MAX_TB;
        if (!get_complete_traceback) {
            addr3 = 2*i;
        }

        num_common_markers = d_num_common_markers[i];
        if (num_common_markers >= SHASTA_MAX_COMMON_MARKERS_PER_READ) {
            num_common_markers = 0;
        }

        __syncthreads();

        for (int p = 0; p < num_common_markers; p++) {
//            uint64_t v =  __ldg(&d_common_markers[addr1+p]);
            uint64_t v =  d_common_markers[addr1+p];
            
            uint64_t l = ((v >> 32) & p_mask);
            uint64_t u = ((v >> 48) & p_mask);
            
            score = 0.01;
            score_pos = p;
            
            int ptr = p + tx;

            bool stop = false;
            __syncthreads();

            while (!stop) {
                ptr -= bs;
                if (ptr >= 0) {
                    uint64_t v1 = d_common_markers[addr1+ptr];
                    uint64_t l1, u1;
                    l1 = ((v1 >> 32) & p_mask);
                    u1 = ((v1 >> 48) & p_mask);
                    float a = l-l1;
                    float b = u-u1;
                    float alpha = fabs(a-b);
                    if ((l1 < l) && (u1 < u) && (u-u1 <= maxSkip) && (l-l1 <= maxSkip) && (alpha <= maxDrift)) {
                        float pscore = d_score[addr2+ptr]+1;
                        if (score < pscore) { 
                            score = pscore;
                            score_pos = ptr;
                        }
                    }

                    if (l > l1+maxSkip)  {
                        stop = true;
                    }
                }
                else {
                    stop = true;
                }
                
                stop = __shfl_sync(0xffffffff, stop, 0);
            }

            __syncthreads();

            // parallel reduction (max)
            for (int s = 1; s <= bs; s *= 2) {
                float val = __shfl_up_sync(0xffffffff, score, s, bs);
                uint32_t val_pos = __shfl_up_sync(0xffffffff, score_pos, s, bs);

                if (tx >= s) {
                    if (val > score) {
                        score = val;
                        score_pos = val_pos; 
                    }
                }
            }
            
            if (tx == bs-1) {
                d_score[addr2+p] = score;
                d_score_pos[addr2+p] = score_pos;
                if (score > max_score) {
                    max_score = score;
                    max_score_pos = p;
                }
            }
            __syncthreads();
        }

        __syncthreads();

        if (tx == bs-1) {
            int num_ptr = 0;
            uint64_t last_common_marker = 0;

            if (max_score > 0) {
                int curr_pos = max_score_pos;
                int prev_pos = max_score_pos + 1;

                while ((curr_pos >= 0) && (prev_pos > curr_pos)) {
                    prev_pos = curr_pos;
                    if (num_ptr < SHASTA_MAX_TB) {
                        if (get_complete_traceback) {
                            d_alignments[addr3+num_ptr] = d_common_markers[addr1+curr_pos];
                        }
                        else {
                            if (num_ptr == 0) {
                                d_alignments[addr3] = d_common_markers[addr1+curr_pos];
                            }
                            else {
                                last_common_marker = d_common_markers[addr1+curr_pos];
                            }
                        }
                    }
                    num_ptr++;
                    curr_pos = d_score_pos[addr2+curr_pos];
                }
            }

            if (num_ptr < SHASTA_MAX_TB) {
                if (get_complete_traceback) {
                    d_alignments[addr3+num_ptr] = 0;
                }
                else {
                    d_alignments[addr3+1] = last_common_marker;
                }
                d_num_traceback[i] = num_ptr;
            }
            else {
                d_alignments[addr3] = 0;
                d_num_traceback[i] = 0;
            }
        }
        __syncthreads();
    }
}

extern "C" std::tuple<int, size_t> shasta_initializeProcessors (size_t numUniqueMarkers) {
    int nDevices;

    num_unique_markers = (uint32_t) numUniqueMarkers;

    cudaError_t err;
    
    err = cudaGetDeviceCount(&nDevices);
    NUM_DEVICES = nDevices;

    if (err != cudaSuccess) {
        throw runtime_error("GPU_ERROR: No GPU device found! Consider running without the --gpu flag.");
    }
    
    size_t device_memory;
    for (int i = 0; i < nDevices; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        device_memory = prop.totalGlobalMem;
        if (device_memory > 0xffffffff) {
            NUM_BLOCKS = (1 << 11);
            BATCH_SIZE = (1 << 12);
        }
        else {
            NUM_BLOCKS = (1 << 8);
            BATCH_SIZE = (1 << 9);
            break;
        }
        //printf("Device Number: %d\n", i);
        //printf("  Device name: %s\n", prop.name);
        //printf("  Memory Clock Rate (KHz): %d\n",
        //prop.memoryClockRate);
        //printf("  Memory Bus Width (bits): %d\n",
        //prop.memoryBusWidth);
        //printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
        //2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    }

    d_alignments = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));

    d_score = (float**) malloc(nDevices*sizeof(float*));
    d_score_pos = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_num_traceback = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_common_markers = (uint64_t**) malloc(nDevices*sizeof(uint64_t*));
    d_num_common_markers = (uint32_t**) malloc(nDevices*sizeof(uint32_t*));
    d_batch_rid_markers = (uint64_t**) malloc(nDevices*sizeof(uint64_t*));
    d_rid_marker_pos = (uint64_t**) malloc(nDevices*sizeof(uint64_t*));

    size_t num_bytes;

    for (int k=0; k<nDevices; k++) {
        
        available_gpus.push_back(k);

        err = cudaSetDevice(k);
        if (err != cudaSuccess) {
            throw runtime_error("GPU_ERROR: could not set device");
        }
        
        num_bytes = BATCH_SIZE*SHASTA_MAX_TB*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_alignments[k], num_bytes); 
        if (err != cudaSuccess) {
            throw runtime_error("GPU_ERROR: cudaMalloc failed!\n");
        }

        num_bytes = NUM_BLOCKS*SHASTA_MAX_COMMON_MARKERS_PER_READ*sizeof(float);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_score[k], num_bytes); 
        if (err != cudaSuccess) {
            throw runtime_error("GPU_ERROR: cudaMalloc failed!\n");
        }
        
        num_bytes = NUM_BLOCKS*SHASTA_MAX_COMMON_MARKERS_PER_READ*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_score_pos[k], num_bytes); 
        if (err != cudaSuccess) {
            throw runtime_error("GPU_ERROR: cudaMalloc failed!\n");
        }
        
        num_bytes = BATCH_SIZE*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_num_traceback[k], num_bytes); 
        if (err != cudaSuccess) {
            throw runtime_error("GPU_ERROR: cudaMalloc failed!\n");
        }
        
        num_bytes = BATCH_SIZE*sizeof(uint32_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_num_common_markers[k], num_bytes); 
        if (err != cudaSuccess) {
            throw runtime_error("GPU_ERROR: cudaMalloc failed!\n");
        }
        
        num_bytes = BATCH_SIZE*SHASTA_MAX_COMMON_MARKERS_PER_READ*sizeof(uint64_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_common_markers[k], num_bytes); 
        if (err != cudaSuccess) {
            throw runtime_error("GPU_ERROR: cudaMalloc failed!\n");
        }

        num_bytes = (1+2*BATCH_SIZE*numUniqueMarkers)*sizeof(uint64_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_batch_rid_markers[k], num_bytes); 
        if (err != cudaSuccess) {
            throw runtime_error("GPU_ERROR: cudaMalloc failed!\n");
        }

        num_bytes = (BATCH_SIZE*SHASTA_MAX_MARKERS_PER_READ)*sizeof(uint64_t);
        if (k==0)
            fprintf(stdout, "\t-Requesting %3.0e bytes on GPU\n", (double)num_bytes);
        err = cudaMalloc(&d_rid_marker_pos[k], num_bytes); 
        if (err != cudaSuccess) {
            throw runtime_error("GPU_ERROR: cudaMalloc failed!\n");
        }

        initialize_batch_rid_markers<<<NUM_BLOCKS, BLOCK_SIZE>>> (d_batch_rid_markers[k], numUniqueMarkers, BATCH_SIZE);  
    }

    return std::make_tuple(nDevices, BATCH_SIZE);
}

extern "C" void shasta_alignBatchGPU (size_t maxMarkerFrequency, size_t maxSkip, size_t maxDrift, size_t n, uint64_t num_pos, uint64_t num_reads, uint64_t* batch_rid_marker_pos, uint64_t* batch_read_pairs, uint32_t* h_alignments, uint32_t* h_num_traceback, bool get_complete_traceback) {
    bool report_time = false;

    int k = -1;

    while (k < 0) {
        std::unique_lock<std::mutex> locker(mu);
        if (available_gpus.empty()) {
            cv.wait(locker, [](){return !available_gpus.empty();});
        }
        k = available_gpus.back();
        available_gpus.pop_back();
        locker.unlock();
    }

    struct timeval t1, t2, t3;
    long useconds, seconds, mseconds;
    
    cudaError_t err; 

    err = cudaSetDevice(k);
    if (err != cudaSuccess) {
        throw runtime_error("GPU_ERROR: could not set device.\n");
    }
    
    gettimeofday(&t1, NULL);

    try {
        err = cudaMemcpy(d_rid_marker_pos[k], batch_rid_marker_pos, num_pos*sizeof(uint64_t), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            throw runtime_error("Error: cudaMemcpy failed!\n");
        }

        thrust::device_ptr<uint64_t> d_batch_rid_markers_ptr (d_batch_rid_markers[k]);
        thrust::device_ptr<uint64_t> d_rid_marker_pos_ptr (d_rid_marker_pos[k]);

        thrust::device_vector<uint64_t> t_d_rid_marker_pos (d_rid_marker_pos_ptr, d_rid_marker_pos_ptr + num_pos);
        thrust::device_vector<uint16_t> t_d_adjusted_pos (num_pos);
        thrust::device_vector<uint64_t> t_d_sorted_rid_marker_pos (num_pos);
        thrust::device_vector<uint64_t> t_d_rid_markers (d_batch_rid_markers_ptr, d_batch_rid_markers_ptr+num_reads*num_unique_markers+1);
        thrust::device_vector<uint64_t> t_d_read_pairs (batch_read_pairs, batch_read_pairs+2*n);
        thrust::device_vector<uint64_t> t_d_index_table (num_reads*num_unique_markers+1);

        thrust::copy (t_d_rid_marker_pos.begin(), t_d_rid_marker_pos.end(), t_d_sorted_rid_marker_pos.begin());
        thrust::sort(t_d_sorted_rid_marker_pos.begin(), t_d_sorted_rid_marker_pos.end());


        thrust::lower_bound(t_d_sorted_rid_marker_pos.begin(),
                t_d_sorted_rid_marker_pos.end(),
                t_d_rid_markers.begin(),
                t_d_rid_markers.end(),
                t_d_index_table.begin());

        gettimeofday(&t2, NULL);

        uint64_t* d_sorted_rid_marker_pos = thrust::raw_pointer_cast (t_d_sorted_rid_marker_pos.data());
        uint16_t* d_adjusted_pos = thrust::raw_pointer_cast (t_d_adjusted_pos.data()); 
        uint64_t* d_rid_marker_pos = thrust::raw_pointer_cast (t_d_rid_marker_pos.data());
        uint64_t* d_index_table = thrust::raw_pointer_cast (t_d_index_table.data());
        uint64_t* d_read_pairs = thrust::raw_pointer_cast (t_d_read_pairs.data());

        skip_high_frequency_markers <<< num_reads, 32>>> (maxMarkerFrequency, num_unique_markers, d_index_table, d_rid_marker_pos, d_sorted_rid_marker_pos, d_adjusted_pos);

        find_common_markers <<<NUM_BLOCKS, BLOCK_SIZE>>> (maxMarkerFrequency, n, num_unique_markers, d_read_pairs, d_index_table, d_rid_marker_pos, d_sorted_rid_marker_pos, d_adjusted_pos, d_num_common_markers[k], d_common_markers[k]);
        
        find_traceback <<<NUM_BLOCKS, BAND_SIZE>>>(n, maxSkip, maxDrift, d_score[k], d_common_markers[k], d_num_common_markers[k], d_score_pos[k], d_alignments[k], d_num_traceback[k], get_complete_traceback);

    }
    catch (std::bad_alloc) {
        throw runtime_error("Insufficient GPU memory. Try on GPU with larger memory or without --gpu option.\n");
    }

    err = cudaMemcpy(h_num_traceback, d_num_traceback[k], n*sizeof(uint32_t), cudaMemcpyDeviceToHost);
    if (err != cudaSuccess) {
        throw runtime_error("Error: cudaMemcpy failed!\n");
    }

    if (get_complete_traceback) {
        err = cudaMemcpy(h_alignments, d_alignments[k], n*SHASTA_MAX_TB*sizeof(uint32_t), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            throw runtime_error("Error: cudaMemcpy failed!\n");
        }
    }
    else {
        err = cudaMemcpy(h_alignments, d_alignments[k], 2*n*sizeof(uint32_t), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            throw runtime_error("Error: cudaMemcpy failed!\n");
        }
    }
    
    {
        std::unique_lock<std::mutex> locker(mu);
        available_gpus.push_back(k);
        locker.unlock();
        cv.notify_one();
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

    return;
}

extern "C" size_t shasta_getGpuBatchSize(){
    return BATCH_SIZE;
}

extern "C" void shasta_shutdownProcessors() {
    cudaDeviceReset();
}
