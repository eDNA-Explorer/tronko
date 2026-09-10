# the compiler: gcc for C program
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    DEFAULT_CC := $(if $(wildcard /opt/homebrew/bin/gcc-15),/opt/homebrew/bin/gcc-15,gcc)
    SANITIZER_CC := /usr/bin/clang
    LSAN_FLAGS := -fsanitize=address
else
    DEFAULT_CC := gcc
    SANITIZER_CC := $(DEFAULT_CC)
    LSAN_FLAGS := -fsanitize=leak
endif
CC = $(DEFAULT_CC)

#compiler flags:
# -g adds debugging information to the executable file
# -Wall turns on most, but not all, compiler warnings
BASE_CFLAGS = -std=gnu99 -w
BASE_DBGCFLAGS = -std=gnu99 -g -w
CFLAGS ?=
CPPFLAGS ?=
LDFLAGS ?=
STACKPROTECT = -fno-stack-protector
# -lm links the math library
#LIBS = -lm -lpthread -lz
LDLIBS ?= -lm -pthread -lz -lzstd

ifeq ($(UNAME_S),Darwin)
    CPPFLAGS += -D_XOPEN_SOURCE=700 -D_DARWIN_C_SOURCE -D_GNU_SOURCE
    ifneq ($(wildcard /opt/homebrew/include/zstd.h),)
        CPPFLAGS += -I/opt/homebrew/include
        LDFLAGS += -L/opt/homebrew/lib
    endif
else
    CPPFLAGS += -D_GNU_SOURCE
    LDLIBS += -lrt
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
SOURCES = tronko-assign.c bwa_context.c candidate_workspace.c readreference.c assignment.c allocateMemoryForResults.c allocatetreememory.c getSequenceinRoot.c placement.c options.c printAlignments.c logger.c resource_monitor.c crash_debug.c symbol_resolver.c
NEEDLEMANWUNSCH = needleman_wunsch.c alignment.c alignment_scoring.c
BWA = bwa_source_files/bntseq.c bwa_source_files/bwa.c bwa_source_files/bwamem.c bwa_source_files/bwamem_extra.c bwa_source_files/bwamem_pair.c bwa_source_files/bwashm.c bwa_source_files/bwt.c bwa_source_files/bwt_gen.c bwa_source_files/bwtindex.c bwa_source_files/fastmap.c bwa_source_files/is.c bwa_source_files/kopen.c bwa_source_files/kstring.c bwa_source_files/ksw.c bwa_source_files/kthread.c bwa_source_files/QSufSort.c bwa_source_files/rle.c bwa_source_files/rope.c bwa_source_files/utils.c
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

$(TARGET): $(SOURCES) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2)
	$(CC) $(CPPFLAGS) $(OPTIMIZATION) $(BASE_CFLAGS) $(CFLAGS) $(ARCH_FLAGS) $(MATH_FLAGS) $(MEMOPT_FLAGS) $(PARQUET_FLAGS) $(STACKPROTECT) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(PARQUET_SOURCES) $(LDFLAGS) $(LDLIBS) $(PARQUET_LIBS) $(OPENMP_FLAGS)

debug: carquet_lib $(TARGET).c
	$(CC) $(CPPFLAGS) $(BASE_DBGCFLAGS) $(CFLAGS) $(ARCH_FLAGS) $(MATH_FLAGS) $(MEMOPT_FLAGS) $(PARQUET_FLAGS) $(STACKPROTECT) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(PARQUET_SOURCES) $(LDFLAGS) $(LDLIBS) $(PARQUET_LIBS) $(OPENMP_FLAGS)

.PHONY: all clean debug asan lsan tsan

asan:
	$(MAKE) CC=$(SANITIZER_CC) TARGET=tronko-assign-asan OPTIMIZATION=-O1 BASE_CFLAGS=-g STACKPROTECT= \
		CFLAGS="$(CFLAGS) -fno-omit-frame-pointer -fsanitize=address,undefined" \
		LDFLAGS="$(LDFLAGS) -fsanitize=address,undefined"

lsan:
	$(MAKE) CC=$(SANITIZER_CC) TARGET=tronko-assign-lsan OPTIMIZATION=-O1 BASE_CFLAGS=-g STACKPROTECT= \
		CFLAGS="$(CFLAGS) -fno-omit-frame-pointer $(LSAN_FLAGS)" \
		LDFLAGS="$(LDFLAGS) $(LSAN_FLAGS)"

tsan:
	$(MAKE) CC=$(SANITIZER_CC) TARGET=tronko-assign-tsan OPTIMIZATION=-O1 BASE_CFLAGS=-g STACKPROTECT= \
		CFLAGS="$(CFLAGS) -fno-omit-frame-pointer -fsanitize=thread" \
		LDFLAGS="$(LDFLAGS) -fsanitize=thread"

clean:
	$(RM) tronko-assign tronko-assign-asan tronko-assign-lsan tronko-assign-tsan
