#ifndef CRC32_H
#define CRC32_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

uint32_t tronko_crc32(uint32_t crc, const void *buf, size_t len);
uint32_t tronko_crc32_file(FILE *fp, long start, long end);

#endif
