#include "../timestamp.hpp"

#define MAX_MARKERS_PER_READ (1 << 15)
#define GPU_BATCH_SIZE (1 << 12)
#define MIN_GPU_BATCH_SIZE (1 << 4)

void common_markers (uint64_t* read_pairs, uint32_t offset, int n);
extern "C" int initializeProcessors();
extern "C" void alignBatchGPU (size_t n, size_t deviceId, size_t maxSkip, size_t num_pos, uint32_t* h_reads, uint64_t* h_read_pairs, uint32_t* h_alignments);
extern "C" void shutdownProcessors(int nDevices);
