#include "../timestamp.hpp"

#define SHASTA_LOG_MAX_MARKERS_PER_READ 16
#define SHASTA_LOG_MAX_TB 13
#define SHASTA_LOG_GPU_BATCH_SIZE 11

#define SHASTA_MAX_MARKERS_PER_READ (1 << SHASTA_LOG_MAX_MARKERS_PER_READ)
#define SHASTA_MAX_TB (1 << SHASTA_LOG_MAX_TB)
#define SHASTA_GPU_BATCH_SIZE (1 << SHASTA_LOG_GPU_BATCH_SIZE)

extern "C" int shasta_initializeProcessors(size_t numUniqueMarkers);
extern "C" void shasta_alignBatchGPU (size_t maxMarkerFrequency, size_t maxSkip, size_t n, uint64_t num_pos, uint64_t num_reads, uint64_t* batch_rid_marker_pos, uint64_t* batch_rid_markers, uint64_t* batch_read_pairs, uint32_t* h_alignments, uint32_t* h_num_traceback);
extern "C" void shasta_shutdownProcessors(int nDevices);
