# the compiler: gcc for C program
# On macOS, 'gcc' is often Apple clang which lacks -fopenmp.
# Auto-detect real GCC from Homebrew; fall back to 'gcc'.
CC := $(shell for v in 15 14 13 12 11; do \
    if command -v gcc-$$v >/dev/null 2>&1; then echo gcc-$$v; exit; fi; done; echo gcc)

#compiler flags:
# -g adds debugging information to the executable file
# -Wall turns on most, but not all, compiler warnings
CFLAGS = -w -Wno-error=incompatible-pointer-types
DBGCFLAGS = -g -w
# -lm links the math library
#LIBS = -lm -lpthread -lz
LIBS = -lm -pthread -lz -std=gnu99 -lgomp
OPTIMIZATION = -O3
# Without an -march flag GCC targets baseline x86-64, i.e. SSE2. Verified on the
# shipped binary: calculateSPArr contained 50 xmm and 0 ymm instructions, so the
# SP inner loop was running 4 ints/vector instead of 8. Also helps likelihood.c.
# Override with `make ARCH=` (or ARCH=-march=x86-64-v3) for a portable binary --
# -march=native ties the result to this machine's CPU.
ARCH ?= -march=native
# GCC's default -ffp-contract=fast fuses a*b+c into FMA, which changes
# likelihood/posterior values and would break golden-output comparison.
FPFLAGS = -ffp-contract=off
NUMCAT ?= 1
NUMCAT_FLAGS = -DNUMCAT=$(NUMCAT)
#sources
SOURCES = tronko-build.c getclade.c readfasta.c readreference.c allocatetreememory.c math.c likelihood.c opt.c options.c printtree.c
HASHMAP = hashmap.c
OBJECTS = (SOURCES: .c = .o)
# the build target executable:
TARGET = tronko-build

all: $(TARGET)
$(TARGET): $(SOURCES) $(HASHMAP)
	$(CC) $(OPTIMIZATION) $(ARCH) $(FPFLAGS) $(NUMCAT_FLAGS) -fopenmp -w -Wno-error -Wno-incompatible-pointer-types -Wno-int-conversion -Wno-implicit-function-declaration -o $(TARGET) $(HASHMAP) $(SOURCES) $(LIBS)
# debug/asan need -fopenmp too: without it the 18 `#pragma omp threadprivate`
# declarations in global.h/opt.h are silently ignored, so those globals become
# genuinely shared and these builds have different data-sharing semantics than
# the production build -- they report races it does not have, and miss ones it does.
debug: $(SOURCES) $(HASHMAP)
	$(CC) -g -O0 $(NUMCAT_FLAGS) -fopenmp -w -Wno-error -Wno-incompatible-pointer-types -Wno-int-conversion -Wno-implicit-function-declaration -o $(TARGET) $(HASHMAP) $(SOURCES) $(LIBS)

asan: $(SOURCES) $(HASHMAP)
	$(CC) -g -O0 -fsanitize=address -fno-omit-frame-pointer -fopenmp \
	    -Wno-incompatible-pointer-types -Wno-int-conversion \
	    -Wno-implicit-function-declaration \
	    -o $(TARGET) $(HASHMAP) $(SOURCES) $(LIBS)

clean:
	$(RM) $(TARGET)

