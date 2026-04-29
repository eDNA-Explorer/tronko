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
NUMCAT ?= 1
NUMCAT_FLAGS = -DNUMCAT=$(NUMCAT)
#sources
SOURCES = tronko-build.c getclade.c readfasta.c readreference.c allocatetreememory.c math.c likelihood.c opt.c options.c printtree.c
HASHMAP = hashmap.c
OBJECTS = (SOURCES: .c = .o)
# the build target executable:
TARGET = tronko-build

all: $(TARGET)
$(TARGET): $(TARGET).c
	$(CC) $(OPTIMIZATION) $(NUMCAT_FLAGS) -fopenmp -w -Wno-error -Wno-incompatible-pointer-types -Wno-int-conversion -Wno-implicit-function-declaration -o $(TARGET) $(HASHMAP) $(SOURCES) $(LIBS)
debug: $(TARGET).c
	$(CC) -g -O0 $(NUMCAT_FLAGS) -w -Wno-error -Wno-incompatible-pointer-types -Wno-int-conversion -Wno-implicit-function-declaration -o $(TARGET) $(HASHMAP) $(SOURCES) $(LIBS)

asan: $(TARGET).c
	$(CC) -g -O0 -fsanitize=address -fno-omit-frame-pointer \
	    -Wno-incompatible-pointer-types -Wno-int-conversion \
	    -Wno-implicit-function-declaration \
	    -o $(TARGET) $(HASHMAP) $(SOURCES) $(LIBS)

clean:
	$(RM) $(TARGET)

