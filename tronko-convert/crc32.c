#include "crc32.h"
#include <stdio.h>

// Standard CRC-32 (IEEE 802.3 polynomial)
static uint32_t crc32_table[256];
static int table_initialized = 0;

static void init_crc32_table(void) {
    if (table_initialized) return;
    for (uint32_t i = 0; i < 256; i++) {
        uint32_t crc = i;
        for (int j = 0; j < 8; j++) {
            crc = (crc >> 1) ^ ((crc & 1) ? 0xEDB88320 : 0);
        }
        crc32_table[i] = crc;
    }
    table_initialized = 1;
}

uint32_t tronko_crc32(uint32_t crc, const void *buf, size_t len) {
    init_crc32_table();
    const uint8_t *p = (const uint8_t *)buf;
    crc = ~crc;
    while (len--) {
        crc = crc32_table[(crc ^ *p++) & 0xFF] ^ (crc >> 8);
    }
    return ~crc;
}

uint32_t tronko_crc32_file(FILE *fp, long start, long end) {
    init_crc32_table();
    long orig_pos = ftell(fp);
    fseek(fp, start, SEEK_SET);

    uint32_t crc = 0;
    uint8_t buf[8192];
    long remaining = end - start;

    while (remaining > 0) {
        size_t to_read = remaining > (long)sizeof(buf) ? sizeof(buf) : (size_t)remaining;
        size_t n = fread(buf, 1, to_read, fp);
        if (n == 0) break;
        crc = tronko_crc32(crc, buf, n);
        remaining -= n;
    }

    fseek(fp, orig_pos, SEEK_SET);
    return crc;
}
