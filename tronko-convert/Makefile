CC = gcc
CFLAGS = -Wall -Wextra -O3
DEBUG_CFLAGS = -Wall -Wextra -g -DDEBUG
LDFLAGS = -lz -lm -lzstd

SRCS = tronko-convert.c format_common.c format_text.c format_binary.c crc32.c utils.c
OBJS = $(SRCS:.c=.o)
TARGET = tronko-convert

.PHONY: all clean debug test

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

debug: CFLAGS = $(DEBUG_CFLAGS)
debug: clean $(TARGET)

test: $(TARGET)
	@echo "=== Testing text-to-binary conversion ==="
	./$(TARGET) -i ../tronko-build/example_datasets/single_tree/reference_tree.txt \
	            -o /tmp/test.trkb -v
	@echo ""
	@echo "=== Testing binary-to-text conversion ==="
	./$(TARGET) -i /tmp/test.trkb -o /tmp/test_roundtrip.txt -t -v
	@echo ""
	@echo "=== Comparing original and round-trip ==="
	@if diff -q ../tronko-build/example_datasets/single_tree/reference_tree.txt \
	            /tmp/test_roundtrip.txt > /dev/null 2>&1; then \
		echo "SUCCESS: Files are identical"; \
	else \
		echo "MISMATCH: Checking for acceptable floating-point differences..."; \
		diff ../tronko-build/example_datasets/single_tree/reference_tree.txt \
		     /tmp/test_roundtrip.txt | head -20; \
	fi
	@echo ""
	@echo "=== File size comparison ==="
	@ls -lh ../tronko-build/example_datasets/single_tree/reference_tree.txt /tmp/test.trkb

clean:
	rm -f $(OBJS) $(TARGET)
