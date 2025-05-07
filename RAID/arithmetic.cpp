/*********************************************************
 * arithmetic.cpp  - implementation for arithmetical functions
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/

#include "arithmetic.h"

#include "misc.h"

#include <emmintrin.h>
#include <immintrin.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>

using namespace std;

/** Get the pointer aligned to ARITHMETIC_ALIGNMENT boundary
 */
unsigned char* AlignedMalloc(size_t Size) {
#ifdef _WIN32
  return (unsigned char*) _aligned_malloc(Size, ARITHMETIC_ALIGNMENT);
#else
  void* pResult;
  if (posix_memalign(&pResult, ARITHMETIC_ALIGNMENT, Size)) {
    return 0;
  } else {
    return (unsigned char*) pResult;
  }
#endif
};

/**Deallocate the aligned block of memory
 */

void AlignedFree(void* ptr) {
#ifdef _WIN32
  _aligned_free(ptr);
#else
  free(ptr);
#endif
};

/**XOR arrays A and B, storing the result in A
Depending on alignment of the arrays, different implementations are used*/
void XOR(unsigned char* pA, const unsigned char* pB, unsigned Size) {
  LOCKEDADD(opXOR, Size);
#ifdef AVX
  for (unsigned i = 0; i < Size; i += 32) {
    __m256* a = (__m256*) (pA + i);
    const __m256* b = (const __m256*) (pB + i);
    *a = _mm256_xor_ps(*a, *b);
  };
#else
  if (size_t(pB) & 0xF) {
    // B is not aligned
    if (!(size_t(pA) & 0xF)) {
      // A is aligned
      for (unsigned i = 0; i < Size; i += 16) {
        __m128i B = _mm_lddqu_si128((const __m128i*) (pB + i));
        *(__m128i*) (pA + i) = _mm_xor_si128(*(__m128i*) (pA + i), B);
      }
    } else {
      // A is not aligned
      for (unsigned i = 0; i < Size; i += 16) {
        __m128i B = _mm_lddqu_si128((const __m128i*) (pB + i));
        __m128i A = _mm_loadu_si128((__m128i*) (pA + i));
        A = _mm_xor_si128(A, B);
        _mm_storeu_si128((__m128i*) (pA + i), A);
      }
    }
  } else {
    // B is aligned
    if (!(size_t(pA) & 0xF)) {
      // everything is aligned
      for (unsigned i = 0; i < Size; i += 16) {
        *(__m128i*) (pA + i) = _mm_xor_si128(*(__m128i*) (pA + i), *(__m128i*) (pB + i));
      }
    } else {
      // A is not aligned
      for (unsigned i = 0; i < Size; i += 16) {
        __m128i B = _mm_lddqu_si128((const __m128i*) (pB + i));
        __m128i A = _mm_loadu_si128((__m128i*) (pA + i));
        A = _mm_xor_si128(A, *(const __m128i*) (pB + i));
        _mm_storeu_si128((__m128i*) (pA + i), A);
      }
    }
  };

#endif
};

/**XOR arrays A and B, storing the result in C
Depending on alignment of the arrays, different implementations are used*/
void XOR(const unsigned char* pA, const unsigned char* pB, unsigned char* pC, unsigned Size) {
  LOCKEDADD(opXOR, Size);
#ifdef AVX
  for (unsigned i = 0; i < Size; i += 32) {
    __m256* a = (__m256*) (pA + i);
    const __m256* b = (const __m256*) (pB + i);
    *a = _mm256_xor_ps(*a, *b);
  };
#else
  if (size_t(pB) & 0xF) {
    // B is not aligned
    if (!(size_t(pA) & 0xF)) {
      // A is aligned
      if (!(size_t(pC) & 0xF)) {
        // C is aligned
        for (unsigned i = 0; i < Size; i += 16) {
          __m128i B = _mm_lddqu_si128((const __m128i*) (pB + i));
          *(__m128i*) (pC + i) = _mm_xor_si128(*(__m128i*) (pA + i), B);
        }
      } else {
        // C is not aligned
        for (unsigned i = 0; i < Size; i += 16) {
          __m128i B = _mm_lddqu_si128((const __m128i*) (pB + i));
          __m128i C = _mm_xor_si128(*(__m128i*) (pA + i), B);
          _mm_storeu_si128((__m128i*) (pC + i), C);
        };
      }
    } else
      // A is not aligned
      if (!(size_t(pC) & 0xF)) {
        // C is aligned
        for (unsigned i = 0; i < Size; i += 16) {
          __m128i B = _mm_lddqu_si128((const __m128i*) (pB + i));
          __m128i A = _mm_loadu_si128((__m128i*) (pA + i));
          *(__m128i*) (pC + i) = _mm_xor_si128(A, B);
        }
      } else {
        // C is not aligned
        for (unsigned i = 0; i < Size; i += 16) {
          __m128i B = _mm_lddqu_si128((const __m128i*) (pB + i));
          __m128i A = _mm_loadu_si128((__m128i*) (pA + i));
          __m128i C = _mm_xor_si128(A, B);
          _mm_storeu_si128((__m128i*) (pC + i), C);
        }
      };
  } else {
    // B is aligned
    if (!(size_t(pA) & 0xF)) {
      // A is aligned
      if (!(size_t(pC) & 0xF)) {
        // everything is aligned
        for (unsigned i = 0; i < Size; i += 16) {
          *(__m128i*) (pC + i) = _mm_xor_si128(*(__m128i*) (pA + i), *(__m128i*) (pB + i));
        };
      } else {
        // everything but C is aligned
        for (unsigned i = 0; i < Size; i += 16) {
          __m128i C = _mm_xor_si128(*(__m128i*) (pA + i), *(__m128i*) (pB + i));
          _mm_storeu_si128((__m128i*) (pC + i), C);
        };
      };
    } else {
      // A is not aligned
      if (!(size_t(pC) & 0xF)) {
        // C is aligned
        for (unsigned i = 0; i < Size; i += 16) {
          __m128i B = _mm_lddqu_si128((const __m128i*) (pB + i));
          __m128i A = _mm_loadu_si128((__m128i*) (pA + i));
          *(__m128i*) (pC + i) = _mm_xor_si128(A, *(const __m128i*) (pB + i));
        }
      } else {
        // C is not aligned
        for (unsigned i = 0; i < Size; i += 16) {
          __m128i B = _mm_lddqu_si128((const __m128i*) (pB + i));
          __m128i A = _mm_loadu_si128((__m128i*) (pA + i));
          __m128i C = _mm_xor_si128(A, B);
          _mm_storeu_si128((__m128i*) (pC + i), C);
        }
      };
    };
  };

#endif
};

/**XOR arrays A, B and C, storing the result in D
Everything must be aligned*/
void XOR(const unsigned char* pA, const unsigned char* pB, const unsigned char* pC, unsigned char* pD, unsigned Size) {
  LOCKEDADD(opXOR, 2 * Size);
  for (unsigned i = 0; i < Size; i += 16) {
    __m128i U = _mm_xor_si128(*(__m128i*) (pA + i), *(__m128i*) (pB + i));
    *(__m128i*) (pD + i) = _mm_xor_si128(U, *(__m128i*) (pC + i));
  };
};

/**XOR arrays A and B, and XOR the result to  C
Everything must be aligned*/
void XORXOR(const unsigned char* pA, const unsigned char* pB, unsigned char* pC, unsigned Size) {
  LOCKEDADD(opXOR, 2 * Size);
#ifdef AVX
  for (unsigned i = 0; i < Size; i += 32) {
    __m256* a = (__m256*) (pA + i);
    const __m256* b = (const __m256*) (pB + i);
    *a = _mm256_xor_ps(*a, *b);
  };
#else
  for (unsigned i = 0; i < Size; i += 16) {
    __m128i Temp = _mm_xor_si128(*(__m128i*) (pA + i), *(__m128i*) (pB + i));
    *(__m128i*) (pC + i) = _mm_xor_si128(*(__m128i*) (pC + i), Temp);
  };

#endif
};

unsigned* GF = 0;
int* LogTable = 0;
MultiplyHelper* pHelper = 0;
int FieldSize_1 = 0;
unsigned Extension = 0;

// generator polynomials for Galois fields
const unsigned GFGenerators[] = {
    0,
    0,
    7,
    0xB,
    0x13,
    0x25,
    0x43,
    0x83,
    0x11D /*, 0x211,02011,04005,010123,020033,21993,39065,85245*/
};

unsigned char multBy2(unsigned char v, unsigned w) {
  if (!v) {
    return 0;
  }
  unsigned char c = LogTable[v] + 1;
  if (c >= (1 << w) - 1) {
    c -= (1 << w) - 1;
  }
  return GF[c + 1];
}

GFValue singleMult(GFValue v1, GFValue v2, unsigned w) {
  if ((v1 == 0) || (v2 == 0)) {
    return 0;
  }
  int c = LogTable[v1] + LogTable[v2];
  if (c >= (1 << w) - 1) {
    c -= (1 << w) - 1;
  }
  return GF[c + 1];
}

/// construct GF(2^m) tables. The primitive polynomial is taken from a table
///@return false if m is too large
void InitGF(unsigned m) {
  if (GF) {
    throw Exception("GF(%d) is already initialized", FieldSize_1 + 1);
  }
  if (m > sizeof(GFGenerators) / sizeof(GFGenerators[0]) - 1) {
    throw Exception("Don't know the primitive  polynomial for GF(2^%d)", m);
  }
  if (sizeof(GFValue) * 8 < m) {
    throw Exception("GFValue is too small", m);
  }
  if (m < 2) {
    throw Exception("Multiplication tables are not needed for GF(2)");
  }

  unsigned GenPoly = GFGenerators[m];
  FieldSize_1 = (1 << m) - 1;
  Extension = m;
  GF = new unsigned[2 * FieldSize_1 + 1];
  LogTable = new int[FieldSize_1 + 1];
  GF[0] = 0;
  LogTable[0] = -1;
  GF[1] = 1;
  LogTable[1] = 0;
  for (int i = 2; i <= FieldSize_1; i++) {
    // compute alpha^i
    GF[i] = GF[i - 1] << 1;
    if ((GF[i] >> m) & 1) {
      GF[i] ^= GenPoly;
    }
    LogTable[GF[i]] = i - 1;
  };
  // extend the table to avoid modulo reduction
  memcpy(GF + FieldSize_1 + 1, GF + 1, sizeof(GF[0]) * FieldSize_1);
  // construct multiplication helpers
  if (m <= 8) {
    pHelper = new MultiplyHelper[FieldSize_1];
    // construct helper for each \alpha^x
    for (int x = 0; x < FieldSize_1; x++) {
      MultiplyHelper& H = pHelper[x];
      H.LookupLow[0] = H.LookupHigh[0] = 0;
      for (GFValue y = 1; y < min(16, FieldSize_1 + 1); y++) {
        int L = LogTable[y] + x;
        H.LookupLow[y] = GF[1 + L];
      };
      for (GFValue y = 1; y < min(16, (FieldSize_1 + 1) >> 4); y++) {
        int L = LogTable[y << 4] + x;
        H.LookupHigh[y] = GF[1 + L];
      };
    };
  }
};

volatile unsigned long long MultiplyCount = 0;

const __m128i Mask0F = _mm_set1_epi8(0x0F);
const __m128i* pMask0F = &Mask0F;

/**
Multiply each value in pSrc by \alpha^x and store the result to pDest.
The implementation is based on the identity (y0+\alpha^4*y1)*x=y0*x+(\alpha^4*y1)*x,
where y0 and y1 are the half-bytes of y=y0+\alpha^4*y1
The products y0*x, \alpha^4*y1*x are stored in helper tables

*/
void Multiply(
    int x, /// scale factor
    const GFValue* pSrc, /// source data block
    GFValue* pDest, /// destination block to be updated
    unsigned Size /// size of the block
) {
  if (x < 0) {
    // nothing to do
    return;
  }
  const MultiplyHelper& H = pHelper[x];
  unsigned PackedCount = Size * sizeof(GFValue) / sizeof(__m128i);
  const __m128i* pmSrc = (const __m128i*) pSrc;
  __m128i* pmDest = (__m128i*) pDest;

  LOCKEDADD(opGFMul, Size);

  if (size_t(pSrc) & 0xF) {
    // source data is not aligned
    if (size_t(pDest) & 0xF) {
      // destination is not aligned
      for (unsigned i = 0; i < PackedCount; i++) {
        // load the source  value
        __m128i A = _mm_lddqu_si128(pmSrc + i);
        // extract the first half byte
        __m128i A0 = _mm_and_si128(A, Mask0F);
        // extract the second half byte
        __m128i A1 = _mm_srli_epi16(A, 4);
        A1 = _mm_and_si128(A1, Mask0F);
        // compute the product of half-bytes by x
        A0 = _mm_shuffle_epi8(H.Lookup0, A0);
        A1 = _mm_shuffle_epi8(H.Lookup1, A1);
        // add them up
        A0 = _mm_xor_si128(A0, A1);
        // write it back
        _mm_storeu_si128(pmDest + i, A0);
      }
    } else {
      // destination is aligned
      for (unsigned i = 0; i < PackedCount; i++) {
        __m128i A = _mm_lddqu_si128(pmSrc + i);
        // extract the first half byte
        __m128i A0 = _mm_and_si128(A, Mask0F);
        // extract the second half byte
        __m128i A1 = _mm_srli_epi16(A, 4);
        A1 = _mm_and_si128(A1, Mask0F);
        // compute the product of half-bytes by x
        A0 = _mm_shuffle_epi8(H.Lookup0, A0);
        A1 = _mm_shuffle_epi8(H.Lookup1, A1);
        // add them up
        pmDest[i] = _mm_xor_si128(A0, A1);
      };
    }
  } else {
    // source data is aligned
    if (size_t(pDest) & 0xF) {
      // destination is not aligned
      for (unsigned i = 0; i < PackedCount; i++) {
        // extract the first half byte of source
        __m128i A0 = _mm_and_si128(pmSrc[i], Mask0F);
        // extract the second half byte of source
        __m128i A1 = _mm_srli_epi16(pmSrc[i], 4);
        A1 = _mm_and_si128(A1, Mask0F);
        // compute the product of half-bytes by x
        A0 = _mm_shuffle_epi8(H.Lookup0, A0);
        A1 = _mm_shuffle_epi8(H.Lookup1, A1);
        // add them up
        A0 = _mm_xor_si128(A0, A1);
        // write it back
        _mm_storeu_si128(pmDest + i, A0);
      }
    } else {
      // everything is aligned
      for (unsigned i = 0; i < PackedCount; i++) {
        // extract the first half byte
        __m128i A0 = _mm_and_si128(pmSrc[i], Mask0F);
        // extract the second half byte
        __m128i A1 = _mm_srli_epi16(pmSrc[i], 4);
        A1 = _mm_and_si128(A1, Mask0F);
        // compute the product of half-bytes by x
        A0 = _mm_shuffle_epi8(H.Lookup0, A0);
        A1 = _mm_shuffle_epi8(H.Lookup1, A1);
        // add them up
        pmDest[i] = _mm_xor_si128(A0, A1);
      };
    }
  };
};

volatile unsigned long long MultiplyAddCount = 0;

/**
Multiply each value in pSrc by \alpha^x and add the result to pDest.
The implementation is based on the identity (y0+\alpha^4*y1)*x=y0*x+(\alpha^4*y1)*x,
where y0 and y1 are the half-bytes of y=y0+\alpha^4*y1
The products y0*x, \alpha^4*y1*x are stored in helper tables

*/
void MultiplyAdd(
    int x, /// scale factor
    const GFValue* pSrc, /// source data block
    GFValue* pDest, /// destination block to be updated
    unsigned Size /// block size
) {
  if (x < 0) {
    // nothing to do
    return;
  }
  if (x == 0) {
    // no multiplication
    XOR(pDest, pSrc, Size);
    return;
  };

  const MultiplyHelper& H = pHelper[x];
  unsigned PackedCount = Size * sizeof(GFValue) / sizeof(__m128i);
  const __m128i* pmSrc = (const __m128i*) pSrc;
  __m128i* pmDest = (__m128i*) pDest;

  LOCKEDADD(opGFMulAdd, Size);

  if (size_t(pSrc) & 0xF) {
    // source data is not aligned
    if (size_t(pDest) & 0xF) {
      // destination is not aligned
      for (unsigned i = 0; i < PackedCount; i++) {
        // load the source  value
        __m128i A = _mm_lddqu_si128(pmSrc + i);
        // load the destination value
        __m128i B = _mm_lddqu_si128(pmDest + i);
        // extract the first half byte
        __m128i A0 = _mm_and_si128(A, Mask0F);
        // extract the second half byte
        __m128i A1 = _mm_srli_epi16(A, 4);
        A1 = _mm_and_si128(A1, Mask0F);
        // compute the product of half-bytes by x
        A0 = _mm_shuffle_epi8(H.Lookup0, A0);
        A1 = _mm_shuffle_epi8(H.Lookup1, A1);
        // add them up
        A0 = _mm_xor_si128(A0, A1);
        // add them to the destination
        B = _mm_xor_si128(B, A0);
        // write it back
        _mm_storeu_si128(pmDest + i, B);
      }
    } else {
      // destination is aligned
      for (unsigned i = 0; i < PackedCount; i++) {
        __m128i A = _mm_lddqu_si128(pmSrc + i);
        // extract the first half byte
        __m128i A0 = _mm_and_si128(A, Mask0F);
        // extract the second half byte
        __m128i A1 = _mm_srli_epi16(A, 4);
        A1 = _mm_and_si128(A1, Mask0F);
        // compute the product of half-bytes by x
        A0 = _mm_shuffle_epi8(H.Lookup0, A0);
        A1 = _mm_shuffle_epi8(H.Lookup1, A1);
        // add them up
        A0 = _mm_xor_si128(A0, A1);
        // add them to the destination
        pmDest[i] = _mm_xor_si128(pmDest[i], A0);
      };
    }
  } else {
    // source data is aligned
    if (size_t(pDest) & 0xF) {
      // destination is not aligned
      for (unsigned i = 0; i < PackedCount; i++) {
        // load the destination value
        __m128i B = _mm_lddqu_si128(pmDest + i);
        // extract the first half byte of source
        __m128i A0 = _mm_and_si128(pmSrc[i], Mask0F);
        // extract the second half byte of source
        __m128i A1 = _mm_srli_epi16(pmSrc[i], 4);
        A1 = _mm_and_si128(A1, Mask0F);
        // compute the product of half-bytes by x
        A0 = _mm_shuffle_epi8(H.Lookup0, A0);
        A1 = _mm_shuffle_epi8(H.Lookup1, A1);
        // add them up
        A0 = _mm_xor_si128(A0, A1);
        // add them to the destination
        B = _mm_xor_si128(B, A0);
        // write it back
        _mm_storeu_si128(pmDest + i, B);
      }
    } else {
      // everything is aligned
      for (unsigned i = 0; i < PackedCount; i++) {
        // extract the first half byte
        __m128i A0 = _mm_and_si128(pmSrc[i], Mask0F);
        // extract the second half byte
        __m128i A1 = _mm_srli_epi16(pmSrc[i], 4);
        A1 = _mm_and_si128(A1, Mask0F);
        // compute the product of half-bytes by x
        A0 = _mm_shuffle_epi8(H.Lookup0, A0);
        A1 = _mm_shuffle_epi8(H.Lookup1, A1);
        // add them up
        A0 = _mm_xor_si128(A0, A1);
        // add them to the destination
        pmDest[i] = _mm_xor_si128(pmDest[i], A0);
      };
    }
  };
};

/**
multiply each value in pSrc by x and add to the results values from pCorrection
pSrc[i]=(pSrc[i]*x)^pCorrection[i]

The implementation is based on the identity (y0+\alpha^4*y1)*x=y0*x+(\alpha^4*y1)*x,
where y0 and y1 are the half-bytes of y=y0+\alpha^4*y1
The products y0*x, \alpha^4*y1*x are stored in helper tables

*/
void AddMultiply(
    int x, /// scale factor
    GFValue* pSrc, /// data block to be updated
    const GFValue* pCorrection, /// correction
    unsigned Size /// block size
) {
  if (x < 0) {
    // nothing to do
    return;
  }
  const MultiplyHelper& H = pHelper[x];
  unsigned PackedCount = Size * sizeof(GFValue) / sizeof(__m128i);
  __m128i* pmSrc = (__m128i*) pSrc;
  const __m128i* pmCorrection = (const __m128i*) pCorrection;

  LOCKEDADD(opGFMulAdd, Size);

  if (size_t(pSrc) & 0xF) {
    // source data is not aligned
    if (size_t(pCorrection) & 0xF) {
      // destination is not aligned
      for (unsigned i = 0; i < PackedCount; i++) {
        // load the source  value
        __m128i A = _mm_lddqu_si128(pmSrc + i);
        // load the destination value
        __m128i B = _mm_lddqu_si128(pmCorrection + i);
        // extract the first half byte
        __m128i A0 = _mm_and_si128(A, Mask0F);
        // extract the second half byte
        __m128i A1 = _mm_srli_epi16(A, 4);
        A1 = _mm_and_si128(A1, Mask0F);
        // compute the product of half-bytes by x
        A0 = _mm_shuffle_epi8(H.Lookup0, A0);
        A1 = _mm_shuffle_epi8(H.Lookup1, A1);
        // add them up
        A0 = _mm_xor_si128(A0, A1);
        // add them to the destination
        B = _mm_xor_si128(B, A0);
        // write it back
        _mm_storeu_si128(pmSrc + i, B);
      }
    } else {
      // destination is aligned
      for (unsigned i = 0; i < PackedCount; i++) {
        __m128i A = _mm_lddqu_si128(pmSrc + i);
        // extract the first half byte
        __m128i A0 = _mm_and_si128(A, Mask0F);
        // extract the second half byte
        __m128i A1 = _mm_srli_epi16(A, 4);
        A1 = _mm_and_si128(A1, Mask0F);
        // compute the product of half-bytes by x
        A0 = _mm_shuffle_epi8(H.Lookup0, A0);
        A1 = _mm_shuffle_epi8(H.Lookup1, A1);
        // add them up
        A0 = _mm_xor_si128(A0, A1);
        // add them to the destination
        pmSrc[i] = _mm_xor_si128(pmCorrection[i], A0);
      };
    }
  } else {
    // source data is aligned
    if (size_t(pCorrection) & 0xF) {
      // destination is not aligned
      for (unsigned i = 0; i < PackedCount; i++) {
        // load the destination value
        __m128i B = _mm_lddqu_si128(pmCorrection + i);
        // extract the first half byte of source
        __m128i A0 = _mm_and_si128(pmSrc[i], Mask0F);
        // extract the second half byte of source
        __m128i A1 = _mm_srli_epi16(pmSrc[i], 4);
        A1 = _mm_and_si128(A1, Mask0F);
        // compute the product of half-bytes by x
        A0 = _mm_shuffle_epi8(H.Lookup0, A0);
        A1 = _mm_shuffle_epi8(H.Lookup1, A1);
        // add them up
        A0 = _mm_xor_si128(A0, A1);
        // add them to the destination
        B = _mm_xor_si128(B, A0);
        // write it back
        _mm_storeu_si128(pmSrc + i, B);
      }
    } else {
      // everything is aligned
      for (unsigned i = 0; i < PackedCount; i++) {
        // extract the first half byte
        __m128i A0 = _mm_and_si128(pmSrc[i], Mask0F);
        // extract the second half byte
        __m128i A1 = _mm_srli_epi16(pmSrc[i], 4);
        A1 = _mm_and_si128(A1, Mask0F);
        // compute the product of half-bytes by x
        A0 = _mm_shuffle_epi8(H.Lookup0, A0);
        A1 = _mm_shuffle_epi8(H.Lookup1, A1);
        // add them up
        A0 = _mm_xor_si128(A0, A1);
        // add them to the destination
        pmSrc[i] = _mm_xor_si128(pmCorrection[i], A0);
      };
    }
  };
};

/** sum the entries of pSrc1 and pSrc2, multiply the product by x and store it in pDest
pDest[i]=(pSrc1[i]+pSrc2[i])*x

The implementation is based on the identity (y0+\alpha^4*y1)*x=y0*x+(\alpha^4*y1)*x,
where y0 and y1 are the half-bytes of y=y0+\alpha^4*y1
The products y0*x, \alpha^4*y1*x are stored in helper tables

*/

void MultiplySum(
    int x, /// scale factor
    const GFValue* pSrc1, /// source data block 1
    const GFValue* pSrc2, /// source data block 1
    GFValue* pDest, /// destination array. Must be aligned
    unsigned Size /// block size
) {
  if (x < 0) {
    memset(pDest, 0, sizeof(GFValue) * Size);
    // nothing to do
    return;
  };
  const MultiplyHelper& H = pHelper[x];
  unsigned PackedCount = Size * sizeof(GFValue) / sizeof(__m128i);
  __m128i* pmSrc1 = (__m128i*) pSrc1;
  __m128i* pmSrc2 = (__m128i*) pSrc2;
  __m128i* pmDest = (__m128i*) pDest;

  LOCKEDADD(opGFMulAdd, Size);

  for (unsigned i = 0; i < PackedCount; i++) {
    // load the source  value
    __m128i A = _mm_lddqu_si128(pmSrc1 + i);
    __m128i B = _mm_lddqu_si128(pmSrc2 + i);
    A = _mm_xor_si128(A, B);
    // extract the first half byte
    __m128i A0 = _mm_and_si128(A, Mask0F);
    // extract the second half byte
    __m128i A1 = _mm_srli_epi16(A, 4);
    A1 = _mm_and_si128(A1, Mask0F);
    // compute the product of half-bytes by x
    A0 = _mm_shuffle_epi8(H.Lookup0, A0);
    A1 = _mm_shuffle_epi8(H.Lookup1, A1);
    // add them up
    A0 = _mm_xor_si128(A0, A1);
    // write it back
    _mm_store_si128(pmDest + i, A0);
  };
};
