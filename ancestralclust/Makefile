# the compiler: gcc for C program
# On macOS, use Homebrew gcc (gcc-15, gcc-14, etc.) since Apple clang lacks -fopenmp
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  CC := $(or $(shell which gcc-15 2>/dev/null),$(shell which gcc-14 2>/dev/null),$(shell which gcc-13 2>/dev/null),gcc)
else
  CC ?= gcc
endif

# -w suppresses warnings; -fpermissive downgrades C23 strict errors (GCC 14+)
CFLAGS = -w -fpermissive -std=gnu99
DBGCFLAGS = -g -w -fpermissive -std=gnu99 -fopenmp
LIBS = -lm -pthread -lz
OPENMP = -fopenmp
OPTIMIZATION = -O3

#sources
SOURCES = ancestralclust.c options.c math.c opt.c
NEEDLEMANWUNSCH = needleman_wunsch.c alignment.c alignment_scoring.c
HASHMAP = hashmap.c
KALIGN = kalign/run_kalign.c kalign/tlmisc.c kalign/tldevel.c kalign/parameters.c kalign/rwalign.c kalign/alignment_parameters.c kalign/idata.c kalign/aln_task.c kalign/bisectingKmeans.c kalign/esl_stopwatch.c kalign/aln_run.c kalign/alphabet.c kalign/pick_anchor.c kalign/sequence_distance.c kalign/euclidean_dist.c kalign/aln_mem.c kalign/tlrng.c kalign/aln_setup.c kalign/aln_controller.c kalign/weave_alignment.c kalign/aln_seqseq.c kalign/aln_profileprofile.c kalign/aln_seqprofile.c
WFA2 = WFA2/wavefront_aligner.c WFA2/wavefront_bialigner.c WFA2/mm_allocator.c WFA2/wavefront_align.c WFA2/wavefront_attributes.c WFA2/wavefront_debug.c WFA2/profiler_timer.c WFA2/profiler_counter.c WFA2/cigar.c WFA2/wavefront_compute.c WFA2/wavefront_slab.c WFA2/wavefront.c WFA2/vector.c WFA2/wavefront_components.c WFA2/bitmap.c WFA2/wavefront_backtrace_buffer.c WFA2/wavefront_pcigar.c WFA2/wavefront_heuristic.c WFA2/wavefront_penalties.c WFA2/wavefront_unialign.c WFA2/wavefront_plot.c WFA2/heatmap.c WFA2/wavefront_backtrace.c WFA2/wavefront_extend.c WFA2/wavefront_compute_affine2p.c WFA2/wavefront_backtrace_offload.c WFA2/wavefront_compute_affine.c WFA2/wavefront_compute_linear.c WFA2/wavefront_compute_edit.c WFA2/string_padded.c WFA2/wavefront_bialign.c

# the build target executable:
TARGET = ancestralclust

all: $(TARGET)
$(TARGET): $(TARGET).c
	$(CC) $(OPTIMIZATION) $(CFLAGS) $(OPENMP) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(KALIGN) $(WFA2) $(SOURCES) $(LIBS)
debug: $(TARGET).c
	$(CC) $(DBGCFLAGS) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(KALIGN) $(WFA2) $(SOURCES) $(LIBS)

clean:
	$(RM) $(TARGET)
