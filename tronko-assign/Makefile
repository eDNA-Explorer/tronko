# the compiler: gcc for C program
CC = gcc

#compiler flags:
# -g adds debugging information to the executable file
# -Wall turns on most, but not all, compiler warnings
CFLAGS = -w -pg
DBGCFLAGS = -g -w
STACKPROTECT = -fno-stack-protector
# -lm links the math library
#LIBS = -lm -lpthread -lz
LIBS = -lm -pthread -lz -lzstd -lrt -std=gnu99

# Base optimization
OPTIMIZATION = -O3

# Architecture-specific optimization (use: make NATIVE_ARCH=1)
ifdef NATIVE_ARCH
    ARCH_FLAGS = -march=native -mtune=native
else
    ARCH_FLAGS =
endif

# Fast math optimization (use: make FAST_MATH=1)
ifdef FAST_MATH
    MATH_FLAGS = -ffast-math
else
    MATH_FLAGS =
endif

# OpenMP support for WFA2 (use: make ENABLE_OPENMP=1)
ifdef ENABLE_OPENMP
    OPENMP_FLAGS = -fopenmp
else
    OPENMP_FLAGS =
endif

# Memory optimization flag (use: make OPTIMIZE_MEMORY=1)
ifdef OPTIMIZE_MEMORY
    MEMOPT_FLAGS = -DOPTIMIZE_MEMORY
else
    MEMOPT_FLAGS =
endif

# Parquet support via Carquet (use: make ENABLE_PARQUET=1)
ifdef ENABLE_PARQUET
    PARQUET_FLAGS = -DENABLE_PARQUET -Icarquet/include
    PARQUET_SOURCES = parquet_writer.c
    PARQUET_LIBS = -Lcarquet/build -lcarquet
else
    PARQUET_FLAGS =
    PARQUET_SOURCES =
    PARQUET_LIBS =
endif
#sources
SOURCES = tronko-assign.c sam_parse.c bwamem3.c readreference.c assignment.c allocateMemoryForResults.c allocatetreememory.c getSequenceinRoot.c placement.c options.c printAlignments.c logger.c resource_monitor.c crash_debug.c symbol_resolver.c nodestore.c
NEEDLEMANWUNSCH = needleman_wunsch.c alignment.c alignment_scoring.c
WFA2 = WFA2/profiler_timer.c WFA2/wavefront_backtrace.c WFA2/wavefront_compute.c WFA2/wavefront_pcigar.c WFA2/bitmap.c WFA2/wavefront_compute_edit.c WFA2/wavefront_backtrace_offload.c  WFA2/wavefront_penalties.c WFA2/cigar.c WFA2/string_padded.c WFA2/wavefront_bialign.c WFA2/wavefront_compute_linear.c  WFA2/wavefront_plot.c WFA2/vector.c WFA2/wavefront_bialigner.c WFA2/heatmap.c WFA2/wavefront_debug.c WFA2/wavefront_slab.c WFA2/wavefront_align.c WFA2/wavefront_aligner.c WFA2/wavefront.c WFA2/wavefront_unialign.c WFA2/wavefront_components.c WFA2/wavefront_extend.c WFA2/mm_allocator.c WFA2/wavefront_attributes.c WFA2/wavefront_compute_affine2p.c WFA2/wavefront_heuristic.c WFA2/profiler_counter.c WFA2/wavefront_backtrace_buffer.c WFA2/wavefront_compute_affine.c
HASHMAP = hashmap.c
OBJECTS = (SOURCES: .c = .o)
# the build target executable:
TARGET = tronko-assign

all: carquet_lib $(TARGET)

# Build Carquet library if ENABLE_PARQUET is set
carquet_lib:
ifdef ENABLE_PARQUET
	@if [ ! -f carquet/build/libcarquet.a ]; then \
		echo "Building Carquet library..."; \
		cd carquet && mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$$(nproc); \
	fi
endif

$(TARGET): $(TARGET).c
	$(CC) $(OPTIMIZATION) $(ARCH_FLAGS) $(MATH_FLAGS) $(MEMOPT_FLAGS) $(PARQUET_FLAGS) $(STACKPROTECT) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(WFA2) $(SOURCES) $(PARQUET_SOURCES) $(LIBS) $(PARQUET_LIBS) $(OPENMP_FLAGS)

debug: carquet_lib $(TARGET).c
	$(CC) $(DBGCFLAGS) $(ARCH_FLAGS) $(MATH_FLAGS) $(MEMOPT_FLAGS) $(PARQUET_FLAGS) $(STACKPROTECT) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(WFA2) $(SOURCES) $(PARQUET_SOURCES) $(LIBS) $(PARQUET_LIBS) $(OPENMP_FLAGS)

# Tests (tests/run_tests.sh, described in tests/README.md). The integration tests run the binary
# this Makefile builds and the BWA-MEM3 binary beside it, so build both first.
# The unit test of the node-store scoring path is built with the same optimisation and
# architecture flags as the binary.
TEST_NODESTORE = tests/unit/test_nodestore
$(TEST_NODESTORE): tests/unit/test_nodestore.c nodestore.c nodestore.h assignment.c getSequenceinRoot.c global.h
	$(CC) $(OPTIMIZATION) $(ARCH_FLAGS) $(STACKPROTECT) -std=gnu99 -I. -o $(TEST_NODESTORE) tests/unit/test_nodestore.c nodestore.c assignment.c getSequenceinRoot.c -lm -pthread

test: $(TARGET) bwa-mem3 $(TEST_NODESTORE)
	bash tests/run_tests.sh

# BWA-MEM3, the aligner tronko-assign runs (a separate program, never linked in): the pinned
# release v0.14.0 (commit 5c1d5e3), cloned and built unchanged by scripts/build_bwa_mem3.sh in
# build/bwa-mem3, and copied beside tronko-assign, where tronko-assign looks for it. The shipped
# binary must be built with clang 19 (the default) or GCC 15:
#   make bwa-mem3                                   clang-19 / clang++-19
#   make bwa-mem3 BWA_MEM3_CC=gcc-15 BWA_MEM3_CXX=g++-15
# For development only, an older compiler (e.g. Ubuntu 24.04's gcc 13):
#   make bwa-mem3 BWA_MEM3_CC=gcc BWA_MEM3_CXX=g++ ALLOW_UNSUPPORTED_COMPILER=1
BWA_MEM3_CC ?= clang-19
BWA_MEM3_CXX ?= clang++-19
ALLOW_UNSUPPORTED_COMPILER ?= 0
bwa-mem3:
	CC="$(BWA_MEM3_CC)" CXX="$(BWA_MEM3_CXX)" ALLOW_UNSUPPORTED_COMPILER="$(ALLOW_UNSUPPORTED_COMPILER)" bash scripts/build_bwa_mem3.sh build/bwa-mem3
	cp build/bwa-mem3/bwa-mem3 bwa-mem3

.PHONY: test

clean:
	$(RM) $(TARGET) bwa-mem3

