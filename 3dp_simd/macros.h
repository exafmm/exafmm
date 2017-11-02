#ifndef macros_h
#define macros_h

// Detect SIMD Byte length of architecture
#if __MIC__ | __AVX512F__
const int SIMD_BYTES = 64;                                      //!< SIMD byte length of MIC and AVX512
#elif __AVX__ | __bgq__
const int SIMD_BYTES = 32;                                      //!< SIMD byte length of AVX and BG/Q
#elif __SSE__ | __sparc_v9__ | _SX
const int SIMD_BYTES = 16;                                      //!< SIMD byte length of SSE, FX, SX
#else
#error no SIMD
#endif

#ifndef __CUDACC__
#define __host__
#define __device__
#define __forceinline__
#endif

#if _SX
#define __attribute__(x)
#endif

// Bluegene/Q and K computer don't have single precision arithmetic
#if __bgq__ | __sparc_v9__
#ifdef EXAFMM_SINGLE
#error Please use double precision for BG/Q, FX10, FX100
#endif
#endif

// Suppress Intel compiler warnings
#if __INTEL_COMPILER
#pragma warning disable 68 111
#endif

#endif
