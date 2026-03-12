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
# Platform-specific libraries
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    PLATFORM_LIBS = -lrt -ldl
else
    PLATFORM_LIBS = -ldl
endif
# Include paths for Homebrew (macOS) or system paths (Linux)
ifeq ($(UNAME_S),Darwin)
    INCLUDE_FLAGS = -I$(shell brew --prefix zstd 2>/dev/null || echo /usr/local)/include
    LDFLAGS_EXTRA = -L$(shell brew --prefix zstd 2>/dev/null || echo /usr/local)/lib
else
    INCLUDE_FLAGS =
    LDFLAGS_EXTRA =
endif
LIBS = -lm -pthread -lz -lzstd $(PLATFORM_LIBS) -std=gnu99
# macOS requires _XOPEN_SOURCE for ucontext.h, plus _DARWIN_C_SOURCE to keep dladdr/Dl_info
ifeq ($(UNAME_S),Darwin)
    XOPEN_FLAGS = -D_XOPEN_SOURCE=600 -D_DARWIN_C_SOURCE
else
    XOPEN_FLAGS =
endif

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
SOURCES = tronko-assign.c readreference.c assignment.c allocateMemoryForResults.c allocatetreememory.c getSequenceinRoot.c placement.c options.c printAlignments.c logger.c resource_monitor.c crash_debug.c symbol_resolver.c minimap2_shim.c minimap2_wrapper.c
NEEDLEMANWUNSCH = needleman_wunsch.c alignment.c alignment_scoring.c
BWA = bwa_source_files/bntseq.c bwa_source_files/bwa.c bwa_source_files/bwamem.c bwa_source_files/bwamem_extra.c bwa_source_files/bwamem_pair.c bwa_source_files/bwashm.c bwa_source_files/bwt.c bwa_source_files/bwt_gen.c bwa_source_files/bwtindex.c bwa_source_files/fastmap.c bwa_source_files/is.c bwa_source_files/kopen.c bwa_source_files/kstring.c bwa_source_files/ksw.c bwa_source_files/kthread.c bwa_source_files/QSufSort.c bwa_source_files/rle.c bwa_source_files/rope.c bwa_source_files/utils.c
WFA2 = WFA2/profiler_timer.c WFA2/wavefront_backtrace.c WFA2/wavefront_compute.c WFA2/wavefront_pcigar.c WFA2/bitmap.c WFA2/wavefront_compute_edit.c WFA2/wavefront_backtrace_offload.c  WFA2/wavefront_penalties.c WFA2/cigar.c WFA2/string_padded.c WFA2/wavefront_bialign.c WFA2/wavefront_compute_linear.c  WFA2/wavefront_plot.c WFA2/vector.c WFA2/wavefront_bialigner.c WFA2/heatmap.c WFA2/wavefront_debug.c WFA2/wavefront_slab.c WFA2/wavefront_align.c WFA2/wavefront_aligner.c WFA2/wavefront.c WFA2/wavefront_unialign.c WFA2/wavefront_components.c WFA2/wavefront_extend.c WFA2/mm_allocator.c WFA2/wavefront_attributes.c WFA2/wavefront_compute_affine2p.c WFA2/wavefront_heuristic.c WFA2/profiler_counter.c WFA2/wavefront_backtrace_buffer.c WFA2/wavefront_compute_affine.c
HASHMAP = hashmap.c
OBJECTS = (SOURCES: .c = .o)
# the build target executable:
TARGET = tronko-assign

# minimap2 library (built as static lib to avoid symbol conflicts with BWA)
MINIMAP2_DIR = minimap2_src
# Rename kseq functions to avoid duplicate symbols with BWA's utils.c
MINIMAP2_COMMON_FLAGS = -DHAVE_KALLOC -I$(MINIMAP2_DIR) -w \
	-Dkseq_init=mm2_kseq_init -Dkseq_read=mm2_kseq_read -Dkseq_destroy=mm2_kseq_destroy
# Core minimap2 sources (excluding main.c, kthread.c, misc.c — conflicts with BWA)
MINIMAP2_CORE_SRCS = $(MINIMAP2_DIR)/kalloc.c $(MINIMAP2_DIR)/bseq.c \
	$(MINIMAP2_DIR)/sketch.c $(MINIMAP2_DIR)/sdust.c $(MINIMAP2_DIR)/options.c \
	$(MINIMAP2_DIR)/index.c $(MINIMAP2_DIR)/lchain.c $(MINIMAP2_DIR)/align.c \
	$(MINIMAP2_DIR)/hit.c $(MINIMAP2_DIR)/seed.c $(MINIMAP2_DIR)/jump.c \
	$(MINIMAP2_DIR)/map.c $(MINIMAP2_DIR)/format.c $(MINIMAP2_DIR)/pe.c \
	$(MINIMAP2_DIR)/esterr.c $(MINIMAP2_DIR)/splitidx.c
MINIMAP2_CORE_OBJS = $(MINIMAP2_CORE_SRCS:.c=.o)
# KSW2 alignment files (need architecture-specific flags)
MINIMAP2_KSW2_SRCS = $(MINIMAP2_DIR)/ksw2_ll_sse.c
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_M),arm64)
    # ARM64 (Apple Silicon / aarch64): use NEON via sse2neon
    MINIMAP2_KSW2_FLAGS = $(MINIMAP2_COMMON_FLAGS) -DKSW_SSE2_ONLY -D__SSE2__ -I$(MINIMAP2_DIR)/sse2neon
    MINIMAP2_KSW2_SRCS += $(MINIMAP2_DIR)/ksw2_extz2_sse.c $(MINIMAP2_DIR)/ksw2_extd2_sse.c $(MINIMAP2_DIR)/ksw2_exts2_sse.c
else
    # x86_64: SSE2 + SSE4.1 dispatch
    MINIMAP2_KSW2_FLAGS = $(MINIMAP2_COMMON_FLAGS) -msse2
    MINIMAP2_KSW2_SRCS += $(MINIMAP2_DIR)/ksw2_extz2_sse.c $(MINIMAP2_DIR)/ksw2_extd2_sse.c $(MINIMAP2_DIR)/ksw2_exts2_sse.c $(MINIMAP2_DIR)/ksw2_dispatch.c
endif
MINIMAP2_KSW2_OBJS = $(MINIMAP2_KSW2_SRCS:.c=.o)
MINIMAP2_ALL_OBJS = $(MINIMAP2_CORE_OBJS) $(MINIMAP2_KSW2_OBJS)

all: carquet_lib libminimap2.a $(TARGET)

# Build Carquet library if ENABLE_PARQUET is set
carquet_lib:
ifdef ENABLE_PARQUET
	@if [ ! -f carquet/build/libcarquet.a ]; then \
		echo "Building Carquet library..."; \
		cd carquet && mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$${CMAKE_PREFIX_PATH:-}" && make -j$$(nproc); \
	fi
endif

# Build minimap2 static library
# Core objects use common flags; KSW2 objects use architecture-specific flags
$(MINIMAP2_CORE_OBJS): %.o: %.c
	$(CC) -c $(OPTIMIZATION) $(MINIMAP2_COMMON_FLAGS) $< -o $@

$(MINIMAP2_KSW2_OBJS): %.o: %.c
	$(CC) -c $(OPTIMIZATION) $(MINIMAP2_KSW2_FLAGS) $< -o $@

libminimap2.a: $(MINIMAP2_ALL_OBJS)
	$(AR) -csru $@ $(MINIMAP2_ALL_OBJS)

$(TARGET): $(TARGET).c libminimap2.a
	$(CC) $(OPTIMIZATION) $(ARCH_FLAGS) $(MATH_FLAGS) $(MEMOPT_FLAGS) $(PARQUET_FLAGS) $(INCLUDE_FLAGS) $(XOPEN_FLAGS) $(STACKPROTECT) -I$(MINIMAP2_DIR) -DHAVE_KALLOC -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(PARQUET_SOURCES) libminimap2.a $(LDFLAGS_EXTRA) $(LIBS) $(PARQUET_LIBS) $(OPENMP_FLAGS)

debug: carquet_lib libminimap2.a $(TARGET).c
	$(CC) $(DBGCFLAGS) $(ARCH_FLAGS) $(MATH_FLAGS) $(MEMOPT_FLAGS) $(PARQUET_FLAGS) $(INCLUDE_FLAGS) $(XOPEN_FLAGS) $(STACKPROTECT) -I$(MINIMAP2_DIR) -DHAVE_KALLOC -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(PARQUET_SOURCES) libminimap2.a $(LDFLAGS_EXTRA) $(LIBS) $(PARQUET_LIBS) $(OPENMP_FLAGS)

clean:
	$(RM) $(TARGET) libminimap2.a $(MINIMAP2_ALL_OBJS)

