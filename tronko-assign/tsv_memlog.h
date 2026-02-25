/*
 * tsv_memlog.h - TSV memory logging for tronko-assign
 * Provides machine-parseable TSV output for memory analysis
 */

#ifndef TSV_MEMLOG_H
#define TSV_MEMLOG_H

#include <stdio.h>
#include <stdarg.h>
#include "resource_monitor.h"

// Write a TSV log entry
static inline void tsv_memlog_write(FILE *f, const char *phase, const char *extra_fmt, ...) {
    if (!f) return;

    resource_stats_t stats;
    get_resource_stats(&stats);

    fprintf(f, "%.3f\t%s\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\t",
            stats.wall_time_sec, phase,
            stats.memory_rss_kb / 1024.0,
            stats.memory_vm_size_kb / 1024.0,
            stats.memory_vm_rss_peak_kb / 1024.0,
            stats.user_time_sec,
            stats.system_time_sec);

    if (extra_fmt && extra_fmt[0] != '\0') {
        va_list args;
        va_start(args, extra_fmt);
        vfprintf(f, extra_fmt, args);
        va_end(args);
    }

    fprintf(f, "\n");
    fflush(f);
}

#define TSV_LOG(file, phase, ...) \
    do { if (file) tsv_memlog_write(file, phase, __VA_ARGS__); } while(0)

#define TSV_LOG_SIMPLE(file, phase) \
    do { if (file) tsv_memlog_write(file, phase, ""); } while(0)

// BWA-specific logging macro for bounds tracking
#define TSV_LOG_BWA(fp, phase, leaf_iter, max_matches, concordant, discordant, dropped) \
    do { \
        if (fp) { \
            resource_stats_t __stats; \
            if (get_resource_stats(&__stats) == 0) { \
                fprintf(fp, "%.3f\t%s\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\tleaf_iter=%d,max=%d,conc=%d,disc=%d,dropped=%d\n", \
                    __stats.wall_time_sec, (phase), \
                    __stats.memory_rss_kb / 1024.0, \
                    __stats.memory_vm_size_kb / 1024.0, \
                    __stats.memory_vm_rss_peak_kb / 1024.0, \
                    __stats.user_time_sec, __stats.system_time_sec, \
                    (leaf_iter), (max_matches), (concordant), (discordant), (dropped)); \
                fflush(fp); \
            } \
        } \
    } while(0)

#endif /* TSV_MEMLOG_H */
