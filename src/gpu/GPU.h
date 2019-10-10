#include "../timestamp.hpp"

#define SHASTA_MAX_MARKERS_PER_READ (1 << 15)
#define SHASTA_GPU_BATCH_SIZE (1 << 12)
#define SHASTA_MIN_GPU_BATCH_SIZE (1 << 4)

extern "C" int shasta_initializeProcessors();
extern "C" void shasta_alignBatchGPU (size_t n, size_t deviceId, size_t maxSkip, size_t num_pos, uint32_t* h_reads, uint64_t* h_read_pairs, uint32_t* h_alignments);
extern "C" void shasta_shutdownProcessors(int nDevices);
