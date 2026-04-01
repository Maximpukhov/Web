#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <emmintrin.h>

#include <smmintrin.h>
// Has _mm_mullo_epi32

// TODO: WHEN GNUC IS DEFINED????

#ifdef __GNUC__
  #define likely(x) __builtin_expect(!!(x),1)
  #define unlikely(x) __builtin_expect(!!(x),0)
#else
  #define likely(x)       (x)
  #define unlikely(x)     (x)
#endif

struct coord { 
  int x,y,z;
} ;


typedef union {
  coord c;
  __m128i vect;
} lcoord __attribute((aligned(64)));

typedef union{
  __m64 m64[2];
  __m128i m128;
} conv128_64 __attribute((aligned(64)));
  
  
static inline __m128i mul4x32bit(__m128i a, __m128i b)
{
#ifdef __SSE4_1__ 
  return _mm_mullo_epi32(a,b);
#else
__m128i tmp1 = _mm_mul_epu32(a,b); /* mul 2,0*/    
__m128i tmp2 = _mm_mul_epu32( _mm_srli_si128(a,4), _mm_srli_si128(b,4)); /* mul 3,1 */    
  return _mm_unpacklo_epi32( _mm_shuffle_epi32(tmp1, _MM_SHUFFLE (0,0,2,0)), _mm_shuffle_epi32(tmp2, _MM_SHUFFLE (0,0,2,0))); /* shuffle results to [63..0] and pack */
#endif
}


#endif