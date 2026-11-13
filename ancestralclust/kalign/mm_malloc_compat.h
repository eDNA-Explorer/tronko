/* mm_malloc compatibility shim for non-x86 platforms (ARM/Apple Silicon).
 * On x86, mm_malloc.h is provided by the compiler (GCC/Clang).
 * On ARM, we implement _mm_malloc/_mm_free using posix_memalign/free. */
#ifndef MM_MALLOC_COMPAT_H
#define MM_MALLOC_COMPAT_H

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
  #include <mm_malloc.h>
#else
  #include <stdlib.h>
  static inline void* _mm_malloc(size_t size, size_t align) {
    void* ptr = NULL;
    posix_memalign(&ptr, align, size);
    return ptr;
  }
  static inline void _mm_free(void* ptr) {
    free(ptr);
  }
#endif

#endif /* MM_MALLOC_COMPAT_H */
