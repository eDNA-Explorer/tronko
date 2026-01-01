---
date: 2026-01-01T12:00:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "GPU Acceleration for tronko-assign via CUDA/OpenCL/Rust"
tags: [research, gpu, cuda, performance, parallelization, tronko-assign]
status: complete
last_updated: 2026-01-01
last_updated_by: Claude
---

# Research: GPU Acceleration for tronko-assign

**Date**: 2026-01-01T12:00:00-08:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

Can we speed up tronko-assign's algorithm via GPU acceleration (CUDA/OpenCL)? Specifically:
- LCA computation is embarrassingly parallel across reads
- Score comparisons across 17K+ trees is a good GPU workload
- What are the options: CUDA basics, thrust library, Rust's wgpu/cuda-rs?

## Summary

**GPU acceleration is highly feasible** for tronko-assign with estimated **5-20x speedup potential** for the core computational loops. The best candidates are:

1. **Score calculation** (`getscore_Arr`) - Log-likelihood lookups are embarrassingly parallel
2. **Wavefront Alignment** - WFA-GPU library exists and could replace WFA2
3. **Cross-tree maximum finding** - Parallel reduction is a textbook GPU operation

**Recommended approach**: Start with CUDA + Thrust for NVIDIA-only deployment, or use cudarc (Rust) if integrating with a Rust port. Existing bioinformatics GPU libraries (WFA-GPU, BEAGLE) provide proven patterns.

---

## Detailed Findings

### 1. Current Algorithm Analysis

#### LCA Computation Pipeline

The LCA computation in tronko-assign follows this flow:

```
Read Input → BWA Alignment → Sequence Alignment → Score Calculation → Vote Collection → LCA Resolution
```

**Key Functions**:

| Function | File:Line | Purpose | GPU Potential |
|----------|-----------|---------|---------------|
| `getscore_Arr()` | `assignment.c:143-210` | Log-likelihood scoring per node | **Excellent** |
| `assignScores_Arr_paired()` | `assignment.c:24-52` | Recursive tree traversal | Moderate |
| `place_paired()` | `placement.c:60` | Score voting and max finding | **Excellent** |
| `getLCA_Arr()` | `tronko-assign.c:98-107` | Depth-based LCA | Moderate |
| `run_bwa()` | `tronko-assign.c:278` | BWA alignment | Separate (BarraCUDA exists) |

#### Score Calculation Inner Loop (`assignment.c:160-198`)

This is the hot path - called for every node × every position:

```c
for (i=0; i<alength; i++){
    if (locQuery[i]=='a' || locQuery[i]=='A')
        score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)];
    else if (locQuery[i]=='c' || locQuery[i]=='C')
        score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 1)];
    // ... G, T
}
```

**GPU suitability**: Each position lookup is independent. This maps perfectly to GPU parallel reduction.

#### Cross-Tree Score Comparison (`placement.c:862-891`)

```c
for (i=0; i<number_of_matches;i++){           // Up to 10 matches
    for (j=...; j<...; j++){                   // Trees per match
        for(k=0; k<2*numspecArr[j]-1; k++){    // All nodes in tree
            if (maximum < nodeScores[i][j][k]){
                maximum = nodeScores[i][j][k];
                // Track best tree/node
            }
        }
    }
}
```

**GPU suitability**: Finding maximum across ~17K trees × thousands of nodes is ideal for parallel reduction. Thrust's `max_element` or custom reduction kernel.

### 2. Current Threading Model Limitations

The existing pthread implementation (`tronko-assign.c:1184-1188`):

- **Coarse-grained**: Work divided at read level, not algorithm level
- **No work stealing**: Each thread processes exactly its assigned range
- **BWA runs single-threaded**: Despite BWA's threading capability, `opt->n_threads=1` is hardcoded
- **Batch-synchronous**: All threads must complete before output

**Memory access patterns**:
- Tree data (`treeArr`, `posteriornc`) is read-only and shared
- Each thread has private result buffers (`nodeScores`, `voteRoot`)

This maps well to GPU: shared read-only data in global/texture memory, per-thread-block results.

### 3. GPU Acceleration Options

#### Option A: CUDA + Thrust (Recommended for NVIDIA)

**Thrust Library** provides high-level parallel algorithms:

```cpp
// Score calculation - transform + reduce
thrust::transform(positions.begin(), positions.end(),
                  scores.begin(),
                  score_functor(posteriornc, query));
float total = thrust::reduce(scores.begin(), scores.end());

// Find maximum across all trees
auto max_iter = thrust::max_element(all_scores.begin(), all_scores.end());
```

**Integration approach**:
1. Keep existing C code in `.c` files
2. Create CUDA kernels in `.cu` files with `extern "C"` wrappers
3. Link with `-lcudart`

**Makefile addition**:
```makefile
NVCC = nvcc
kernels.o: kernels.cu
    $(NVCC) -O3 -c kernels.cu -o kernels.o
tronko-assign: $(OBJS) kernels.o
    $(CC) $(OBJS) kernels.o -lcudart -lpthread -lz -o tronko-assign
```

#### Option B: WFA-GPU (Drop-in for Alignment)

[WFA-GPU](https://github.com/quim0/WFA-GPU) is a CUDA implementation of Wavefront Alignment:

- "Allocates central diagonals of wavefronts in shared memory"
- Batch-based processing via `wfagpu_set_batch_size()`
- Could replace WFA2 in `placement.c`

**Integration**: Batch all sequence alignments per chunk, send to GPU, retrieve results.

#### Option C: cudarc (Rust)

If pursuing a Rust port, [cudarc](https://github.com/coreylowman/cudarc) provides:

- Safe Rust wrappers for CUDA driver API
- Access to cuBLAS, cuDNN, cuRAND
- Runtime kernel compilation via NVRTC

```rust
use cudarc::driver::*;
use cudarc::nvrtc::compile_ptx;

let ptx = compile_ptx(KERNEL_SRC)?;
let dev = CudaDevice::new(0)?;
dev.load_ptx(ptx, "score", &["score_kernel"])?;

let scores_gpu = dev.htod_copy(scores)?;
let f = dev.get_func("score", "score_kernel").unwrap();
unsafe { f.launch(cfg, (&scores_gpu, n)) }?;
```

#### Option D: wgpu (Cross-Platform)

For non-NVIDIA hardware, [wgpu](https://github.com/gfx-rs/wgpu) provides:

- WebGPU-based compute shaders
- Works on Vulkan, Metal, DirectX 12
- No access to cuBLAS etc. - must implement algorithms from scratch

**Trade-off**: Portability vs. performance. CUDA will be faster on NVIDIA hardware.

### 4. Specific GPU Kernel Designs

#### Kernel 1: Batch Score Calculation

```cuda
__global__ void score_kernel(
    float *scores,           // Output: [num_reads × num_nodes]
    const float *posteriornc, // Tree posteriors [num_nodes × num_positions × 4]
    const char *queries,      // Query sequences [num_reads × max_len]
    const int *positions,     // Position mappings [num_reads × max_len]
    int num_reads,
    int num_nodes,
    int max_len
) {
    int read_idx = blockIdx.x;
    int node_idx = threadIdx.x + blockIdx.y * blockDim.x;

    if (read_idx >= num_reads || node_idx >= num_nodes) return;

    float score = 0.0f;
    for (int i = 0; i < max_len; i++) {
        int pos = positions[read_idx * max_len + i];
        char nuc = queries[read_idx * max_len + i];
        int nuc_idx = (nuc == 'A' || nuc == 'a') ? 0 :
                      (nuc == 'C' || nuc == 'c') ? 1 :
                      (nuc == 'G' || nuc == 'g') ? 2 : 3;
        score += posteriornc[node_idx * num_positions * 4 + pos * 4 + nuc_idx];
    }
    scores[read_idx * num_nodes + node_idx] = score;
}
```

**Launch configuration**:
- Grid: `(num_reads, ceil(num_nodes/256))`
- Block: `(256)` threads per block

#### Kernel 2: Parallel Maximum Finding

```cuda
// Using Thrust for simplicity
thrust::device_vector<float> d_scores(all_scores);
thrust::device_vector<int> d_indices(num_scores);
thrust::sequence(d_indices.begin(), d_indices.end());

// Find max and its index
auto max_iter = thrust::max_element(d_scores.begin(), d_scores.end());
int max_idx = max_iter - d_scores.begin();
float max_val = *max_iter;
```

#### Kernel 3: Vote Collection (Threshold-based)

```cuda
__global__ void vote_kernel(
    int *votes,              // Output: [num_nodes] = 0 or 1
    const float *scores,     // Input: [num_nodes]
    float max_score,
    float cinterval,
    int num_nodes
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_nodes) return;

    float score = scores[idx];
    votes[idx] = (score >= max_score - cinterval &&
                  score <= max_score + cinterval) ? 1 : 0;
}
```

### 5. Memory Transfer Optimization

**Key principle**: Minimize host-device transfers.

**Strategy for tronko-assign**:

1. **Tree data** (`posteriornc`): Transfer once at startup, keep on GPU
   - Size: ~17K trees × ~1000 nodes × ~300 positions × 4 nucleotides × 4 bytes = ~80GB
   - **Problem**: May exceed GPU memory. Solution: Process tree batches or use multi-GPU

2. **Query sequences**: Transfer per batch (pinned memory for 2x bandwidth)
   ```c
   cudaMallocHost(&h_queries, batch_size * max_len);  // Pinned
   cudaMemcpyAsync(d_queries, h_queries, size, cudaMemcpyHostToDevice, stream);
   ```

3. **Results**: Transfer back after all kernels complete
   - Only need final LCA indices and scores, not intermediate arrays

**Async streams for overlap**:
```c
cudaStream_t stream1, stream2;
// Stream 1: Transfer batch N, compute batch N-1
// Stream 2: Transfer results batch N-2, write to disk
```

### 6. Existing Bioinformatics GPU Libraries

| Library | Purpose | Relevance |
|---------|---------|-----------|
| [WFA-GPU](https://github.com/quim0/WFA-GPU) | Wavefront Alignment | Direct replacement for WFA2 |
| [BEAGLE](https://beagle-dev.github.io) | Phylogenetic likelihood | Proven GPU likelihood calculations |
| [BarraCUDA](https://seqbarracuda.sourceforge.net/) | BWA on GPU | Could replace BWA alignment step |
| [GASAL2](https://github.com/nahmedraja/GASAL2) | Sequence alignment | 21x faster than CPU Parasail |

**BEAGLE is particularly relevant**: It accelerates phylogenetic likelihood calculations similar to tronko-assign's scoring.

### 7. Performance Estimates

Based on similar bioinformatics GPU implementations:

| Component | Current (CPU) | Expected (GPU) | Speedup |
|-----------|---------------|----------------|---------|
| Score calculation | O(reads × nodes × len) | O(reads × len / 1000) | **10-50x** |
| Maximum finding | O(trees × nodes) | O(log(trees × nodes)) | **100x** |
| Wavefront alignment | O(n²) per pair | O(n²/1000) batched | **5-20x** |
| Memory transfer overhead | N/A | ~10-20% of compute | -10-20% |

**Realistic overall speedup**: **5-20x** depending on batch sizes and tree count.

### 8. Implementation Roadmap

#### Phase 1: Score Calculation Kernel (Highest Impact)
1. Create `kernels.cu` with `extern "C"` wrapper
2. Implement `score_kernel` for batch scoring
3. Modify `assignScores_Arr_paired()` to batch and call GPU
4. Benchmark vs. CPU baseline

#### Phase 2: Parallel Reduction
1. Use Thrust for maximum finding across trees
2. GPU-based vote collection with threshold

#### Phase 3: WFA-GPU Integration
1. Replace WFA2 calls with WFA-GPU batch alignment
2. Requires buffering alignments per batch

#### Phase 4: Full Pipeline
1. Async streams for transfer/compute overlap
2. Multi-GPU support for large tree databases
3. Consider BarraCUDA for BWA step

---

## Code References

- `tronko-assign/assignment.c:143-210` - `getscore_Arr()` hot loop
- `tronko-assign/assignment.c:24-52` - `assignScores_Arr_paired()` tree traversal
- `tronko-assign/placement.c:862-891` - Cross-tree maximum finding
- `tronko-assign/placement.c:911-921` - Vote collection with Cinterval
- `tronko-assign/tronko-assign.c:1184-1188` - pthread creation
- `tronko-assign/global.h:85-93` - Node structure with `posteriornc`
- `tronko-assign/global.h:14-22` - `type_of_PP` (float vs double)

---

## Architecture Insights

1. **Data layout is GPU-friendly**: `posteriornc` as 1D array with `PP_IDX(pos, nuc)` macro already enables coalesced memory access
2. **Tree data is read-only during assignment**: Safe for GPU constant/texture memory
3. **Batch processing already exists**: pthread batch model maps to GPU batch model
4. **Float option exists**: `OPTIMIZE_MEMORY` flag switches to float, which is better for GPU

---

## Historical Context (from thoughts/)

No existing GPU research found. Related performance optimization work:

- `thoughts/shared/research/2025-12-29-rust-port-feasibility.md` - Discusses rayon parallelism, SIMD with block-aligner (5-10x potential)
- `thoughts/shared/research/2025-12-30-rust-port-detailed-implementation.md` - Thread-safe patterns, parallel scoring
- `thoughts/shared/plans/2025-12-30-bulk-allocation-1d-layout.md` - Cache-friendly memory layout (relevant for GPU coalescing)

---

## Related Research

- `thoughts/shared/research/2025-12-29-tronko-assign-data-flow.md` - Complete data flow analysis
- `thoughts/shared/research/2025-12-29-tronko-assign-streaming-architecture.md` - Threading model

---

## Open Questions

1. **GPU memory limits**: 17K trees with full `posteriornc` data may exceed GPU memory. How to partition?
2. **Multi-GPU**: Would NCCL help for very large databases?
3. **CPU/GPU hybrid**: Should BWA stay on CPU while scoring goes to GPU?
4. **OpenCL vs CUDA**: Is cross-vendor support needed, or is NVIDIA-only acceptable?
5. **Rust integration**: If pursuing Rust port, should GPU be native Rust (cudarc) or CUDA with FFI?

---

## Recommendations

1. **Start with CUDA + Thrust** for fastest path to results
2. **Target score calculation first** - highest impact, cleanest isolation
3. **Use pinned memory and async streams** from the start
4. **Consider WFA-GPU** if alignment is a significant bottleneck
5. **Benchmark early** - actual speedup depends on tree sizes and batch characteristics

---

## External Sources

- [NVIDIA CUDA Programming Guide](https://docs.nvidia.com/cuda/cuda-c-programming-guide/)
- [Thrust Documentation](https://nvidia.github.io/cccl/thrust/)
- [WFA-GPU Paper](https://academic.oup.com/bioinformatics/article/39/12/btad701/7425447)
- [WFA-GPU GitHub](https://github.com/quim0/WFA-GPU)
- [BEAGLE Library](https://beagle-dev.github.io)
- [cudarc Documentation](https://docs.rs/cudarc/latest/cudarc/)
- [wgpu Documentation](https://docs.rs/wgpu/latest/wgpu/)
