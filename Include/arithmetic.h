/*********************************************************
 * arithmetic.h  - header file for arithmetical functions
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/
#ifndef ARITHMETIC_H
#define ARITHMETIC_H

#include "config.h"

#include <emmintrin.h>
#include <immintrin.h>

#ifdef AVX
// use AVX code
#define ARITHMETIC_ALIGNMENT 32
typedef __m256 MMType;
#define MMXOR(A, B) _mm256_xor_ps(A, B)
#else
// use SSE3 code
#define ARITHMETIC_ALIGNMENT 16
typedef __m128i MMType;
#define MMXOR(A, B) _mm_xor_si128(A, B)
#endif

#define MMSTORE(DestAddr, Src) *(MMType*) (DestAddr) = Src
#define MMLOAD(addr) *(const MMType*) (addr)

/// memory allocation as needed by the arithmetic functions
unsigned char* AlignedMalloc(size_t Size);

/// memory deallocation for aligned pointers
void AlignedFree(void* ptr);

/// XOR arrays A and B, storing the result in A
/// the size of the arrays must be a multiple of 32
/// the arrays must be aligned to ARITHMETIC_ALIGNMENT boundary
void XOR(unsigned char* pA, const unsigned char* pB, unsigned Size);

/// OR arrays A, B and C, storing the result in D
/// Everything must be aligned
void XOR(const unsigned char* pA, const unsigned char* pB, const unsigned char* pC, unsigned char* pD, unsigned Size);

/// XOR arrays A and B, storing the result in C
void XOR(const unsigned char* pA, const unsigned char* pB, unsigned char* pC, unsigned Size);

/// XOR arrays A and B, and XOR the result to  C
void XORXOR(const unsigned char* pA, const unsigned char* pB, unsigned char* pC, unsigned Size);

/// this is GF(2^m) arithmetic
/// it supports m<=8
typedef unsigned char GFValue;

GFValue multBy2(GFValue v, unsigned w);

GFValue singleMult(GFValue v1, GFValue v2, unsigned w);

/// construct GF(2^m) tables. The primitive polynomial is taken from a table
///@return false if m is too large
void InitGF(unsigned m);

/// multiply each value in pSrc by \alpha^x and store the result in pDest
void Multiply(
    int x, /// scale factor
    const GFValue* pSrc, /// source data block
    GFValue* pDest, /// destination block
    unsigned UnitSize /// size of the blocks
);

/// lookup tables needed to implement multiplication by x
struct MultiplyHelper {
  union {
    __m128i Lookup0;
    /// y*x for all possible y=\sum_{i=0}^3 y_i\alpha^i
    GFValue LookupLow[16];
  };

  union {
    __m128i Lookup1;
    /// y*x for all possible y=\sum_{i=4}^7 y_i\alpha^i
    GFValue LookupHigh[16];
  };
};

extern MultiplyHelper* pHelper;
extern const __m128i* pMask0F;

/// elementary multiplication
inline void Multiply(
    int x, /// scale factor
    const __m128i& A,
    __m128i& Dest
) {
  const MultiplyHelper& H = pHelper[x];
  // extract the first half byte
  __m128i A0 = _mm_and_si128(A, *pMask0F);
  // extract the second half byte
  __m128i A1 = _mm_srli_epi16(A, 4);
  A1 = _mm_and_si128(A1, *pMask0F);
  // compute the product of half-bytes by x
  A0 = _mm_shuffle_epi8(H.Lookup0, A0);
  A1 = _mm_shuffle_epi8(H.Lookup1, A1);
  // add them up
  Dest = _mm_xor_si128(A0, A1);
};

/// multiply each value in pSrc by \alpha^x and add the result to pDest
// pDest[i]=pDest[i]^pSrc[i]*x
void MultiplyAdd(
    int x, /// scale factor
    const GFValue* pSrc, /// source data block
    GFValue* pDest, /// destination block to be updated
    unsigned UnitSize /// size of the blocks
);

/// multiply each value in pSrc by x and add to the results values from pCorrection
// pSrc[i]=(pSrc[i]*x)^pCorrection[i]
void AddMultiply(
    int x, /// scale factor
    GFValue* pSrc, /// data block to be updated
    const GFValue* pCorrection, /// correction
    unsigned Size /// block size
);

/// sum the entries of pSrc1 and pSrc2, multiply the product by x and store it in pDest
// pDest[i]=(pSrc1[i]+pSrc2[i])*x
void MultiplySum(
    int x, /// scale factor
    const GFValue* pSrc1, /// source data block 1
    const GFValue* pSrc2, /// source data block 1
    GFValue* pDest, /// destination array. Must be aligned
    unsigned Size /// block size
);

extern unsigned* GF;
extern int* LogTable;
extern int FieldSize_1;
extern unsigned Extension;

#endif
