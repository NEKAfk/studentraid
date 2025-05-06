#include "RDP.h"

#include "arithmetic.h"
#include "misc.h"
#include <algorithm>
#include<array>
#include <cstring>
#include <iostream>
#include <numeric>
#include <vector>

#include <ranges>

static long long mod(long long n, long long m) {
    return (n % m + m) % m;
}

CRDPProcessor::CRDPProcessor(RDPParams *pParams)
    : CRAIDProcessor(pParams->CodeDimension + MAX_ERASURES,
                     pParams->PrimeNumber - 1, pParams, sizeof(*pParams)),
      m_PrimeNumber(pParams->PrimeNumber), m_pRDPXORBuffer(nullptr),
      BLOCK_SIZE(m_StripeUnitSize * m_PrimeNumber), BUFFER_SIZE(BLOCK_SIZE * BLOCKS_PER_BUFFER) {
    if (m_StripeUnitSize % ARITHMETIC_ALIGNMENT) {
        throw Exception("Stripe size must be a multiple of #ARITHMETIC_ALIGNMENT");
    }
}

void CRDPProcessor::iterative_restore(long long diag_step, long long diag_symbol, long long horizontal_symbol,
                                      std::array<unsigned char *, 2> RecoverResults,
                                      std::array<unsigned char *, MAX_ERASURES> pXORBuffer,
                                      long long diag, long long step, unsigned diag_restore,
                                      unsigned horizontal_restore, bool restore_missing_diag = false) const {
    long long row = mod(diag - diag_step * diag_symbol, m_PrimeNumber);
    while (diag != m_PrimeNumber - 1) {
        if (RecoverResults[diag_restore]) {
            memcpy(RecoverResults[diag_restore] + row * m_StripeUnitSize,
                   pXORBuffer[diag_step] + diag * m_StripeUnitSize, m_StripeUnitSize);
        }
        if (restore_missing_diag && (row + (MAX_ERASURES - diag_step) * diag_symbol) % m_PrimeNumber !=
             m_PrimeNumber - 1) {
             XOR(pXORBuffer[MAX_ERASURES - diag_step]
                 + ((row + (MAX_ERASURES - diag_step) * diag_symbol) % m_PrimeNumber) *
                 m_StripeUnitSize,
                 RecoverResults[diag_restore] + row * m_StripeUnitSize,
                 m_StripeUnitSize);
         }
        XOR(pXORBuffer[0] + row * m_StripeUnitSize, pXORBuffer[diag_step] + diag * m_StripeUnitSize,
            m_StripeUnitSize);
        if (RecoverResults[horizontal_restore]) {
            memcpy(RecoverResults[horizontal_restore] + row * m_StripeUnitSize,
                   pXORBuffer[0] + row * m_StripeUnitSize, m_StripeUnitSize);
        }
        if (restore_missing_diag && (row + (MAX_ERASURES - diag_step) * horizontal_symbol) % m_PrimeNumber !=
             m_PrimeNumber - 1) {
             XOR(pXORBuffer[MAX_ERASURES - diag_step]
                 + ((row + (MAX_ERASURES - diag_step) * horizontal_symbol) % m_PrimeNumber)
                 * m_StripeUnitSize,
                 RecoverResults[horizontal_restore] + row * m_StripeUnitSize,
                 m_StripeUnitSize);
         }
        diag = mod(step + diag, m_PrimeNumber);
        if (diag != m_PrimeNumber - 1) {
            XOR(pXORBuffer[diag_step] + diag * m_StripeUnitSize, pXORBuffer[0] + row * m_StripeUnitSize,
                m_StripeUnitSize);
        }
        row = mod(step + row, m_PrimeNumber);
    }
}

CRDPProcessor::~CRDPProcessor() {
    AlignedFree(m_pRDPXORBuffer);
}

bool CRDPProcessor::Attach(CDiskArray *pArray, unsigned ConcurrentThreads) {
    m_pRDPXORBuffer = AlignedMalloc(
        ConcurrentThreads * BUFFER_SIZE);
    return CRAIDProcessor::Attach(pArray, ConcurrentThreads);
}

void CRDPProcessor::ResetErasures() {
    CRAIDProcessor::ResetErasures();
}

bool CRDPProcessor::CheckCodeword(unsigned long long StripeID, unsigned ErasureSetID, size_t ThreadID) {
    if (GetNumOfErasures(ErasureSetID)) {
        return true;
    }
    std::array<unsigned char *, MAX_ERASURES> pXORBuffer{};
    for (unsigned i = 0; i < MAX_ERASURES; i++) {
        pXORBuffer[i] = m_pRDPXORBuffer + (i + ThreadID * BLOCKS_PER_BUFFER) * BLOCK_SIZE;
    }
    unsigned char *pReadBuffer = m_pRDPXORBuffer + (MAX_ERASURES + ThreadID * BLOCKS_PER_BUFFER) * BLOCK_SIZE;
    bool Result = ReadStripeUnit(StripeID, ErasureSetID, 0, 0, m_StripeUnitsPerSymbol, pXORBuffer[0]);
    for (unsigned i = 1; i < MAX_ERASURES; i++) {
        memcpy(pXORBuffer[i], pXORBuffer[0], m_StripeUnitsPerSymbol * m_StripeUnitSize);
    }
    for (long long i = 1; i < m_PrimeNumber; i++) {
        Result &= ReadStripeUnit(StripeID, ErasureSetID, i, 0, m_StripeUnitsPerSymbol, pReadBuffer);
        for (long long k = 0; k < m_StripeUnitsPerSymbol; k++) {
            for (long long j = 0; j < MAX_ERASURES; j++) {
                long long Diag = (k + j * i) % m_PrimeNumber;
                if (Diag == m_PrimeNumber - 1) continue;
                XOR(pXORBuffer[j] + Diag * m_StripeUnitSize,
                    pReadBuffer + k * m_StripeUnitSize,
                    m_StripeUnitSize);
            }
        }
    }
    for (unsigned j = 1; j < MAX_ERASURES; j++) {
        Result &= ReadStripeUnit(StripeID, ErasureSetID, m_PrimeNumber - 1 + j, 0, m_StripeUnitsPerSymbol,
                                 pReadBuffer);
        XOR(pXORBuffer[j], pReadBuffer, m_StripeUnitSize * m_StripeUnitsPerSymbol);
    }
    if (!Result) {
        return false;
    } else {
        return std::ranges::all_of(pXORBuffer,
                                   [this](auto res) {
                                       return std::accumulate(res, res + m_StripeUnitSize * m_StripeUnitsPerSymbol,
                                                              0,
                                                              [](unsigned char lhs, unsigned char rhs) {
                                                                  return lhs | rhs;
                                                              }) == 0;
                                   });
    }
}

bool CRDPProcessor::EncodeStripe(unsigned long long StripeID, unsigned ErasureSetID, const unsigned char *pData,
                                 size_t ThreadID) {
    std::array<unsigned char *, MAX_ERASURES> pXORBuffer{};
    pXORBuffer[0] = std::ranges::all_of(std::views::iota(0u, MAX_ERASURES),
                                        [&](unsigned i) {
                                            return IsErased(ErasureSetID, m_PrimeNumber - 1 + i);
                                        })
                        ? nullptr
                        : m_pRDPXORBuffer + ThreadID * BLOCKS_PER_BUFFER * BLOCK_SIZE;
    for (int i = 1; i < MAX_ERASURES; i++) {
        pXORBuffer[i] = IsErased(ErasureSetID, m_PrimeNumber - 1 + i)
                            ? nullptr
                            : m_pRDPXORBuffer + (i + ThreadID * BLOCKS_PER_BUFFER) * BLOCK_SIZE;
    }

    bool Result = true;
    if (!IsErased(ErasureSetID, 0)) {
        Result &= WriteStripeUnit(StripeID, ErasureSetID, 0, 0, m_StripeUnitsPerSymbol, pData);
    }
    for (unsigned i = 0; i < MAX_ERASURES; i++) {
        if (pXORBuffer[i]) memcpy(pXORBuffer[i], pData, m_StripeUnitsPerSymbol * m_StripeUnitSize);
    }
    pData += m_StripeUnitSize * m_StripeUnitsPerSymbol;
    for (unsigned i = 1; i < m_Dimension; i++, pData += m_StripeUnitSize * m_StripeUnitsPerSymbol) {
        if (!IsErased(ErasureSetID, i)) {
            Result &= WriteStripeUnit(StripeID, ErasureSetID, i, 0, m_StripeUnitsPerSymbol, pData);
        }
        for (long long j = 0; j < MAX_ERASURES; j++) {
            if (!pXORBuffer[j])
                continue;
            for (long long k = 0; k < m_StripeUnitsPerSymbol; k++) {
                long long Diag = (k + j * i) % m_PrimeNumber;
                if (Diag == m_PrimeNumber - 1) continue;
                XOR(pXORBuffer[j] + Diag * m_StripeUnitSize,
                    pData + k * m_StripeUnitSize,
                    m_StripeUnitSize);
            }
        }
    }
    if (!IsErased(ErasureSetID, m_PrimeNumber - 1)) {
        Result &= WriteStripeUnit(StripeID, ErasureSetID, m_PrimeNumber - 1, 0, m_StripeUnitsPerSymbol,
                                  pXORBuffer[0]);
    }
    for (long long j = 1; j < MAX_ERASURES; j++) {
        if (IsErased(ErasureSetID, m_PrimeNumber - 1 + j))
            continue;
        for (long long k = 0; k < m_StripeUnitsPerSymbol; k++) {
            long long Diag = (k + j * (m_PrimeNumber - 1)) % m_PrimeNumber;
            if (Diag == m_PrimeNumber - 1) continue;
            XOR(pXORBuffer[j] + Diag * m_StripeUnitSize,
                pXORBuffer[0] + k * m_StripeUnitSize,
                m_StripeUnitSize);
        }
        Result &= WriteStripeUnit(StripeID, ErasureSetID, m_PrimeNumber - 1 + j, 0, m_StripeUnitsPerSymbol,
                                  pXORBuffer[j]);
    }
    return Result;
}

bool CRDPProcessor::DecodeDataSubsymbols(unsigned long long StripeID, unsigned ErasureSetID, unsigned SymbolID,
                                         unsigned SubsymbolID, unsigned Subsymbols2Decode, unsigned char *pDest,
                                         size_t ThreadID) {
    const unsigned NumberOfErasures = GetNumOfErasures(ErasureSetID);
    if (!IsErased(ErasureSetID, SymbolID)) {
        return ReadStripeUnit(StripeID, ErasureSetID, SymbolID, SubsymbolID, Subsymbols2Decode, pDest);
    }

    bool Result = true;
    unsigned char *pReadBuffer = m_pRDPXORBuffer + (MAX_ERASURES + ThreadID * BLOCKS_PER_BUFFER) * BLOCK_SIZE;
    memset(pDest, 0, m_StripeUnitSize * Subsymbols2Decode);
    long long diag_step = SymbolID < m_PrimeNumber ? 0 : SymbolID - m_PrimeNumber + 1;

    if (NumberOfErasures == 1 ||
        std::ranges::all_of(std::views::iota(0u, NumberOfErasures),
                            [&](unsigned i) {
                                return GetErasedPosition(ErasureSetID, i) >= m_PrimeNumber;
                            }) ||
        (SymbolID < m_PrimeNumber &&
         std::ranges::count_if(std::views::iota(0u, NumberOfErasures),
                               [&](unsigned i) {
                                   return GetErasedPosition(ErasureSetID, i) >= m_PrimeNumber;
                               }) == NumberOfErasures - 1)) {
        for (unsigned i = 0; i < m_PrimeNumber; ++i) {
            if (IsErased(ErasureSetID, i))
                continue;

            for (long long j = SubsymbolID; j < SubsymbolID + Subsymbols2Decode; j++) {
                long long Shift = mod(j - diag_step * i, m_PrimeNumber);
                if (Shift == m_PrimeNumber - 1)
                    continue;
                Result &= ReadStripeUnit(StripeID, ErasureSetID, i, Shift, 1, pReadBuffer);
                XOR(pDest + (j - SubsymbolID) * m_StripeUnitSize, pReadBuffer, m_StripeUnitSize);
            }
        }
        return Result;
    } else if (std::ranges::count_if(std::views::iota(0u, NumberOfErasures),
                                     [&](unsigned i) {
                                         return GetErasedPosition(ErasureSetID, i) >= m_PrimeNumber;
                                     }) == NumberOfErasures - 1) {
        long long k = GetErasedPosition(ErasureSetID, std::ranges::min(std::views::iota(0u, NumberOfErasures),
                                                                       {}, [&](unsigned i) {
                                                                           return GetErasedPosition(ErasureSetID, i);
                                                                       }));
        for (long long i = 0; i < m_PrimeNumber; ++i) {
            if (IsErased(ErasureSetID, i))
                continue;

            for (long long j = SubsymbolID; j < SubsymbolID + Subsymbols2Decode; j++) {
                long long Shift = mod(j - diag_step * i, m_PrimeNumber);
                if (Shift == m_PrimeNumber - 1)
                    continue;
                Result &= ReadStripeUnit(StripeID, ErasureSetID, i, Shift, 1, pReadBuffer);
                XOR(pDest + (j - SubsymbolID) * m_StripeUnitSize, pReadBuffer, m_StripeUnitSize);

                Shift = mod(j - diag_step * k, m_PrimeNumber);
                if (Shift == m_PrimeNumber - 1)
                    continue;

                Result &= ReadStripeUnit(StripeID, ErasureSetID, i, Shift, 1, pReadBuffer);
                XOR(pDest + (j - SubsymbolID) * m_StripeUnitSize, pReadBuffer, m_StripeUnitSize);
            }
        }
        return Result;
    } else {
        unsigned char *TmpResult = m_pRDPXORBuffer + (MAX_ERASURES + 1 + ThreadID * BLOCKS_PER_BUFFER) * BLOCK_SIZE;
        Result = DecodeDataSymbols(StripeID, ErasureSetID, SymbolID, 1, TmpResult, ThreadID);
        memcpy(pDest, TmpResult + SubsymbolID * m_StripeUnitSize, Subsymbols2Decode * m_StripeUnitSize);
        return Result;
    }
}

bool CRDPProcessor::ReadAndCalcSyndroms(unsigned long long StripeID, unsigned ErasureSetID, unsigned SymbolID,
                                        unsigned Symbols2Decode, unsigned char *pDest, size_t ThreadID,
                                        std::array<unsigned char *, MAX_ERASURES> pXORBuffer,
                                        bool CalcMissingDiag = false) {
    bool Result = true;
    unsigned char *pReadBuffer = m_pRDPXORBuffer + (ThreadID * BLOCKS_PER_BUFFER + MAX_ERASURES) * BLOCK_SIZE;

    for (long long i = 0; i < m_Length; ++i) {
        if (IsErased(ErasureSetID, i)) {
            continue;
        }
        unsigned char *pBuffer = nullptr;
        if (i >= SymbolID && i < SymbolID + Symbols2Decode) {
            pBuffer = pDest + (i - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol;
        } else {
            if (i >= m_PrimeNumber && !pXORBuffer[i - m_PrimeNumber + 1])
                continue;
            pBuffer = pReadBuffer;
        }
        Result &= ReadStripeUnit(StripeID, ErasureSetID, i, 0, m_StripeUnitsPerSymbol, pBuffer);

        if (i >= m_PrimeNumber) {
            XOR(pXORBuffer[i - m_PrimeNumber + 1], pBuffer, m_StripeUnitSize * m_StripeUnitsPerSymbol);
            if (CalcMissingDiag) {
                for (long long j = 0; j < m_StripeUnitsPerSymbol; ++j) {
                    XOR(pXORBuffer[i - m_PrimeNumber + 1] + (m_PrimeNumber - 1) * m_StripeUnitSize,
                        pBuffer + j * m_StripeUnitSize, m_StripeUnitSize);
                }
            }
        } else {
            for (long long j = 0; j < MAX_ERASURES; j++) {
                if (!pXORBuffer[j])
                    continue;
                for (long long k = 0; k < m_StripeUnitsPerSymbol; k++) {
                    long long Diag = (k + i * j) % m_PrimeNumber;
                    if (!CalcMissingDiag && Diag == m_PrimeNumber - 1)
                        continue;
                    XOR(
                        pXORBuffer[j] + Diag * m_StripeUnitSize,
                        pBuffer + k * m_StripeUnitSize,
                        m_StripeUnitSize);
                }
            }
        }
    }
    return Result;
}

bool CRDPProcessor::DecodeOneErasure(unsigned long long StripeID, unsigned ErasureSetID, unsigned SymbolID,
                                     unsigned Symbols2Decode, unsigned char *pDest, size_t ThreadID,
                                     long long ErasurePosition = -1) {
    if (ErasurePosition == -1) {
        ErasurePosition = GetErasedPosition(ErasureSetID, 0);
    }
    long long Shift = ErasurePosition < m_PrimeNumber ? 0 : ErasurePosition - m_PrimeNumber + 1;
    std::array<unsigned char *, MAX_ERASURES> pXORBuffer{};
    pXORBuffer[Shift] = pDest + (ErasurePosition - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol;
    memset(pXORBuffer[Shift], 0, m_StripeUnitSize * m_StripeUnitsPerSymbol);
    return ReadAndCalcSyndroms(StripeID, ErasureSetID, SymbolID, Symbols2Decode, pDest, ThreadID, pXORBuffer);
}

bool CRDPProcessor::DecodeTwoErasures(unsigned long long StripeID, unsigned ErasureSetID, unsigned SymbolID,
                                      unsigned Symbols2Decode, unsigned char *pDest, size_t ThreadID,
                                      long long diag_step = 1, std::array<int, 2> ErasurePosition = {-1, -1}) {
    bool Result = true;
    if (ErasurePosition[0] == -1) {
        ErasurePosition = {GetErasedPosition(ErasureSetID, 0), GetErasedPosition(ErasureSetID, 1)};
    }
    std::ranges::sort(ErasurePosition);
    std::array<unsigned char *, MAX_ERASURES> pXORBuffer{};

    if (ErasurePosition[0] >= m_PrimeNumber ||
        (ErasurePosition[1] >= m_PrimeNumber && (
             ErasurePosition[1] < SymbolID || ErasurePosition[1] >= SymbolID + Symbols2Decode))) {
        for (auto pos: ErasurePosition) {
            if ((pos >= SymbolID) && (pos < SymbolID + Symbols2Decode)) {
                long long idx = pos < m_PrimeNumber ? 0 : pos - m_PrimeNumber + 1;
                pXORBuffer[idx] = pDest + (pos - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol;
                memset(pXORBuffer[idx], 0, m_StripeUnitSize * m_StripeUnitsPerSymbol);
            }
        }
        return ReadAndCalcSyndroms(StripeID, ErasureSetID, SymbolID, Symbols2Decode, pDest, ThreadID, pXORBuffer);
    } else if (ErasurePosition[1] >= m_PrimeNumber) {
        pXORBuffer[0] = (ErasurePosition[0] >= SymbolID) && (ErasurePosition[0] < SymbolID + Symbols2Decode)
                            ? pDest + (ErasurePosition[0] - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol
                            : m_pRDPXORBuffer + ThreadID * BLOCKS_PER_BUFFER * BLOCK_SIZE;
        memset(pXORBuffer[0], 0, m_StripeUnitSize * m_StripeUnitsPerSymbol);
        pXORBuffer[ErasurePosition[1] - m_PrimeNumber + 1] =
                pDest + (ErasurePosition[1] - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol;
        memset(pXORBuffer[ErasurePosition[1] - m_PrimeNumber + 1], 0, m_StripeUnitSize * m_StripeUnitsPerSymbol);
        Result = ReadAndCalcSyndroms(StripeID, ErasureSetID, SymbolID, Symbols2Decode, pDest, ThreadID, pXORBuffer);
        for (long long j = 0; j < m_StripeUnitsPerSymbol; j++) {
            long long Diag = (j + (ErasurePosition[1] - m_PrimeNumber + 1) * ErasurePosition[0]) % m_PrimeNumber;
            if (Diag == m_PrimeNumber - 1)
                continue;
            XOR(
                pXORBuffer[ErasurePosition[1] - m_PrimeNumber + 1] + Diag * m_StripeUnitSize,
                pXORBuffer[0] + j * m_StripeUnitSize,
                m_StripeUnitSize);
        }
        return Result;
    }

    pXORBuffer[0] = m_pRDPXORBuffer + ThreadID * BLOCKS_PER_BUFFER * BLOCK_SIZE;
    memset(pXORBuffer[0], 0, m_StripeUnitSize * m_StripeUnitsPerSymbol);
    pXORBuffer[diag_step] = m_pRDPXORBuffer + (ThreadID * BLOCKS_PER_BUFFER + diag_step) * BLOCK_SIZE;
    memset(pXORBuffer[diag_step], 0, m_StripeUnitSize * m_StripeUnitsPerSymbol);
    Result = ReadAndCalcSyndroms(StripeID, ErasureSetID, SymbolID, Symbols2Decode, pDest, ThreadID, pXORBuffer);
    const std::array<unsigned char *, 2> RecoverResults = {
        (ErasurePosition[0] >= SymbolID) && (ErasurePosition[0] < SymbolID + Symbols2Decode)
            ? pDest + (ErasurePosition[0] - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol
            : nullptr,
        (ErasurePosition[1] >= SymbolID) && (ErasurePosition[1] < SymbolID + Symbols2Decode)
            ? pDest + (ErasurePosition[1] - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol
            : nullptr
    };

    long long step = diag_step * (ErasurePosition[1] - ErasurePosition[0]);

    iterative_restore(diag_step, ErasurePosition[1], ErasurePosition[0], RecoverResults,
                      pXORBuffer, mod(diag_step * ErasurePosition[0] - 1, m_PrimeNumber), -step, 1, 0);
    iterative_restore(diag_step, ErasurePosition[0], ErasurePosition[1], RecoverResults,
                      pXORBuffer, mod(diag_step * ErasurePosition[1] - 1, m_PrimeNumber), step, 0, 1);

    return Result;
}

bool CRDPProcessor::DecodeThreeErasures(unsigned long long StripeID, unsigned ErasureSetID, unsigned SymbolID,
                                        unsigned Symbols2Decode, unsigned char *pDest, size_t ThreadID) {
    bool Result = true;
    std::array<int, 3> ErasurePosition = {
        GetErasedPosition(ErasureSetID, 0), GetErasedPosition(ErasureSetID, 1), GetErasedPosition(ErasureSetID, 2)
    };
    std::ranges::sort(ErasurePosition);

    std::array<unsigned char *, MAX_ERASURES> pXORBuffer{};

    if (ErasurePosition[1] >= m_PrimeNumber) {
        if (std::none_of(ErasurePosition.begin() + 1, ErasurePosition.end(), [&](int i) {
            return i >= SymbolID && i < SymbolID + Symbols2Decode;
        })) {
            return DecodeOneErasure(StripeID, ErasureSetID, SymbolID, Symbols2Decode, pDest, ThreadID,
                                    ErasurePosition[0]);
        } else {
            pXORBuffer[0] = (ErasurePosition[0] >= SymbolID) && (ErasurePosition[0] < SymbolID + Symbols2Decode)
                                ? pDest + (ErasurePosition[0] - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol
                                : m_pRDPXORBuffer + ThreadID * BLOCKS_PER_BUFFER * BLOCK_SIZE;
            memset(pXORBuffer[0], 0, m_StripeUnitSize * m_StripeUnitsPerSymbol);
            for (unsigned i = 1; i < MAX_ERASURES; i++) {
                if (ErasurePosition[i] >= SymbolID && ErasurePosition[i] < SymbolID + Symbols2Decode) {
                    pXORBuffer[ErasurePosition[i] - m_PrimeNumber + 1] =
                            pDest + (ErasurePosition[i] - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol;
                    memset(pXORBuffer[ErasurePosition[i] - m_PrimeNumber + 1], 0,
                           m_StripeUnitSize * m_StripeUnitsPerSymbol);
                }
            }
            Result = ReadAndCalcSyndroms(StripeID, ErasureSetID, SymbolID, Symbols2Decode, pDest, ThreadID, pXORBuffer);
            for (long long i = 1; i < MAX_ERASURES; i++) {
                if (!pXORBuffer[i])
                    continue;
                for (long long j = 0; j < m_StripeUnitsPerSymbol; j++) {
                    long long Diag = (j + i * ErasurePosition[0]) % m_PrimeNumber;
                    if (Diag == m_PrimeNumber - 1)
                        continue;
                    XOR(
                        pXORBuffer[i] + Diag * m_StripeUnitSize,
                        pXORBuffer[0] + j * m_StripeUnitSize,
                        m_StripeUnitSize);
                }
            }
            return Result;
        }
    } else if (ErasurePosition[2] >= m_PrimeNumber) {
        long long diag_step = ErasurePosition[2] == m_PrimeNumber ? 2 : 1;
        if (ErasurePosition[2] < SymbolID || ErasurePosition[2] >= SymbolID + Symbols2Decode) {
            return DecodeTwoErasures(StripeID, ErasureSetID, SymbolID, Symbols2Decode, pDest, ThreadID, diag_step,
                                     {ErasurePosition[0], ErasurePosition[1]});
        } else {
            pXORBuffer[0] = m_pRDPXORBuffer + ThreadID * BLOCKS_PER_BUFFER * BLOCK_SIZE;
            memset(pXORBuffer[0], 0, m_StripeUnitSize * m_StripeUnitsPerSymbol);
            pXORBuffer[diag_step] = m_pRDPXORBuffer + (ThreadID * BLOCKS_PER_BUFFER + diag_step) * BLOCK_SIZE;
            memset(pXORBuffer[diag_step], 0, m_StripeUnitSize * m_StripeUnitsPerSymbol);
            pXORBuffer[MAX_ERASURES - diag_step] =
                    pDest + (ErasurePosition[2] - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol;
            memset(pXORBuffer[MAX_ERASURES - diag_step], 0, m_StripeUnitSize * m_StripeUnitsPerSymbol);

            Result = ReadAndCalcSyndroms(StripeID, ErasureSetID, SymbolID, Symbols2Decode, pDest, ThreadID, pXORBuffer);
            const std::array<unsigned char *, 2> RecoverResults = {
                (ErasurePosition[0] >= SymbolID) && (ErasurePosition[0] < SymbolID + Symbols2Decode)
                    ? pDest + (ErasurePosition[0] - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol
                    : nullptr,
                (ErasurePosition[1] >= SymbolID) && (ErasurePosition[1] < SymbolID + Symbols2Decode)
                    ? pDest + (ErasurePosition[1] - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol
                    : nullptr
            };

            long long step = diag_step * (ErasurePosition[1] - ErasurePosition[0]);
            iterative_restore(diag_step, ErasurePosition[1], ErasurePosition[0], RecoverResults, pXORBuffer,
                              mod(diag_step * ErasurePosition[0] - 1, m_PrimeNumber), -step, 1, 0, true);
            iterative_restore(diag_step, ErasurePosition[0], ErasurePosition[1], RecoverResults, pXORBuffer,
                              mod(diag_step * ErasurePosition[1] - 1, m_PrimeNumber), step, 0, 1, true);
            return Result;
        }
    }

    for (unsigned i = 0; i < MAX_ERASURES; i++) {
        pXORBuffer[i] = m_pRDPXORBuffer + (i + ThreadID * BLOCKS_PER_BUFFER) * BLOCK_SIZE;
        memset(pXORBuffer[i], 0, BLOCK_SIZE);
    }

    Result = ReadAndCalcSyndroms(StripeID, ErasureSetID, SymbolID, Symbols2Decode, pDest, ThreadID, pXORBuffer, true);

    long long l{}, m{}, r{};
    for (unsigned i = 0; i < MAX_ERASURES; i++) {
        if (ErasurePosition[i] >= SymbolID && ErasurePosition[i] < SymbolID + Symbols2Decode) {
            m = ErasurePosition[i];
            l = i == 0 ? ErasurePosition[1] : ErasurePosition[0];
            r = i == MAX_ERASURES - 1 ? ErasurePosition[1] : ErasurePosition[MAX_ERASURES - 1];
            break;
        }
    }
    std::vector<unsigned> fst(m_PrimeNumber, 0), snd(m_PrimeNumber, 0), prod(m_PrimeNumber, 0);
    long long idx1 = m_PrimeNumber - 1, idx2 = m_PrimeNumber - 1;
    long long i1 = mod(std::min(m, l) - 1, m_PrimeNumber), i2 = mod(std::min(r, m) - 1, m_PrimeNumber);

    for (unsigned s = 1; s <= m_PrimeNumber - 1; ++s) {
        fst[mod(idx1 - 2 * abs(m - l), m_PrimeNumber)] =
                fst[idx1] ^ (i1 == 0 ? 1u : 0u) ^ (mod(i1 - abs(m - l), m_PrimeNumber) == 0 ? 1u : 0u);
        snd[mod(idx2 - 2 * abs(r - m), m_PrimeNumber)] =
                snd[idx2] ^ (i2 == 0 ? 1u : 0u) ^ (mod(i2 - abs(r - m), m_PrimeNumber) == 0 ? 1u : 0u);
        idx1 = mod(idx1 - 2 * abs(m - l), m_PrimeNumber);
        idx2 = mod(idx2 - 2 * abs(r - m), m_PrimeNumber);
        i1 = mod(i1 - 2 * abs(m - l), m_PrimeNumber);
        i2 = mod(i2 - 2 * abs(r - m), m_PrimeNumber);
    }
    for (long long i = 0; i < m_PrimeNumber; i++) {
        for (long long j = 0; j < m_PrimeNumber; j++) {
            prod[mod(i + j, m_PrimeNumber)] = prod[mod(i + j, m_PrimeNumber)] ^ (fst[i] & snd[j]);
        }
    }
    if (prod[m_PrimeNumber - 1]) {
        std::ranges::for_each(prod, [](auto &&e) { e = 1 - e; });
    }

    unsigned char *Recover = pDest + (m - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol;
    memset(Recover, 0, m_StripeUnitSize * m_StripeUnitsPerSymbol);

    unsigned Missing = prod[mod(-1 - (l + r), m_PrimeNumber)];
    for (long long i = 0; i < m_PrimeNumber - 1; ++i) {
        if (prod[mod(i - (l + r), m_PrimeNumber)] ^ Missing) {
            for (long long j = 0; j < m_StripeUnitsPerSymbol; ++j) {
                XOR(Recover + j * m_StripeUnitSize,
                    pXORBuffer[0] + mod(j - i, m_PrimeNumber) * m_StripeUnitSize,
                    m_StripeUnitSize);
                XOR(Recover + j * m_StripeUnitSize,
                    pXORBuffer[0] + (m_PrimeNumber - 1 - i) * m_StripeUnitSize,
                    m_StripeUnitSize);
            }
        }
    }

    Missing = prod[mod(-1 - l, m_PrimeNumber)] ^ prod[mod(-1 - r, m_PrimeNumber)];
    for (long long i = 0; i < m_PrimeNumber - 1; ++i) {
        if ((prod[mod(i - l, m_PrimeNumber)] ^ prod[mod(i - r, m_PrimeNumber)]) ^ Missing) {
            for (long long j = 0; j < m_StripeUnitsPerSymbol; ++j) {
                XOR(Recover + j * m_StripeUnitSize,
                    pXORBuffer[1] + mod(j - i, m_PrimeNumber) * m_StripeUnitSize,
                    m_StripeUnitSize);
                XOR(Recover + j * m_StripeUnitSize,
                    pXORBuffer[1] + (m_PrimeNumber - 1 - i) * m_StripeUnitSize,
                    m_StripeUnitSize);
            }
        }
    }

    for (long long i = 0; i < m_PrimeNumber - 1; ++i) {
        if (prod[i]) {
            for (long long j = 0; j < m_StripeUnitsPerSymbol; ++j) {
                XOR(Recover + j * m_StripeUnitSize,
                    pXORBuffer[2] + mod(j - i, m_PrimeNumber) * m_StripeUnitSize,
                    m_StripeUnitSize);
                XOR(Recover + j * m_StripeUnitSize,
                    pXORBuffer[2] + (m_PrimeNumber - 1 - i) * m_StripeUnitSize,
                    m_StripeUnitSize);
            }
        }
    }

    if ((l < SymbolID || l >= SymbolID + Symbols2Decode) && (r < SymbolID || r >= SymbolID + Symbols2Decode)) {
        return Result;
    }

    XOR(pXORBuffer[0], Recover, m_StripeUnitSize * m_StripeUnitsPerSymbol);
    for (long long i = 0; i < m_StripeUnitsPerSymbol; ++i) {
        long long Diag = (i + m) % m_PrimeNumber;
        if (Diag == m_PrimeNumber - 1)
            continue;
        XOR(pXORBuffer[1] + Diag * m_StripeUnitSize, Recover + i * m_StripeUnitSize, m_StripeUnitSize);
    }

    long long diag_step = 1;

    const std::array<unsigned char *, 2> RecoverResults = {
        (l >= SymbolID) && (l < SymbolID + Symbols2Decode)
            ? pDest + (l - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol
            : nullptr,
        (r >= SymbolID) && (r < SymbolID + Symbols2Decode)
            ? pDest + (r - SymbolID) * m_StripeUnitSize * m_StripeUnitsPerSymbol
            : nullptr
    };

    long long step = diag_step * (r - l);
    iterative_restore(diag_step, r, l, RecoverResults, pXORBuffer, mod(diag_step * l - 1, m_PrimeNumber), -step, 1, 0);
    iterative_restore(diag_step, l, r, RecoverResults, pXORBuffer, mod(diag_step * r - 1, m_PrimeNumber), step, 0, 1);
    return Result;
}

bool CRDPProcessor::DecodeDataSymbols(unsigned long long StripeID, unsigned ErasureSetID, unsigned SymbolID,
                                      unsigned Symbols2Decode, unsigned char *pDest, size_t ThreadID) {
    const unsigned NumberOfErasures = GetNumOfErasures(ErasureSetID);
    if (std::ranges::none_of(std::views::iota(SymbolID, SymbolID + Symbols2Decode),
                             [&](unsigned i) {
                                 return IsErased(ErasureSetID, i);
                             })) {
        bool Result = true;
        for (unsigned S = SymbolID; S < SymbolID + Symbols2Decode;
             S++, pDest += m_StripeUnitSize * m_StripeUnitsPerSymbol) {
            Result &= ReadStripeUnit(StripeID, ErasureSetID, S, 0, m_StripeUnitsPerSymbol, pDest);
        }
        return Result;
    }

    switch (NumberOfErasures) {
        case 1:
            return DecodeOneErasure(StripeID, ErasureSetID, SymbolID, Symbols2Decode, pDest, ThreadID);
        case 2:
            return DecodeTwoErasures(StripeID, ErasureSetID, SymbolID, Symbols2Decode, pDest, ThreadID);
        case 3:
            return DecodeThreeErasures(StripeID, ErasureSetID, SymbolID, Symbols2Decode, pDest, ThreadID);
        default:
            return false;
    }
}


bool CRDPProcessor::GetEncodingStrategy(unsigned ErasureSetID, unsigned StripeUnitID, unsigned Subsymbols2Encode) {
    return true; //CRAIDProcessor::GetEncodingStrategy(ErasureSetID, StripeUnitID, Subsymbols2Encode);
}

bool CRDPProcessor::UpdateInformationSymbols(unsigned long long StripeID, unsigned ErasureSetID, unsigned StripeUnitID,
                                             unsigned Units2Update, const unsigned char *pData, size_t ThreadID) {
    bool Result = true;
    if (std::ranges::all_of(std::views::iota(0u, MAX_ERASURES), [&](unsigned i) {
        return IsErased(ErasureSetID, m_PrimeNumber - 1 + i);
    })) {
        unsigned WritedUnits = std::min(m_StripeUnitsPerSymbol - StripeUnitID % m_StripeUnitsPerSymbol, Units2Update);
        unsigned CurrentSymbol = StripeUnitID / m_StripeUnitsPerSymbol;
        if (WritedUnits != 0 && !IsErased(ErasureSetID, CurrentSymbol)) {
            Result &= WriteStripeUnit(StripeID, ErasureSetID, CurrentSymbol,
                                      StripeUnitID % m_StripeUnitsPerSymbol, WritedUnits, pData);
        }
        pData += WritedUnits * m_StripeUnitSize;
        CurrentSymbol++;
        Units2Update -= WritedUnits;
        for (unsigned i = 0; i < Units2Update / m_StripeUnitsPerSymbol; ++i,
                                                                        pData += m_StripeUnitSize *
                                                                        m_StripeUnitsPerSymbol, CurrentSymbol++) {
            if (!IsErased(ErasureSetID, CurrentSymbol)) {
                Result &= WriteStripeUnit(StripeID, ErasureSetID, CurrentSymbol,
                                          0, m_StripeUnitsPerSymbol, pData);
            }
        }
        Units2Update %= m_StripeUnitsPerSymbol;
        if (Units2Update != 0) {
            if (!IsErased(ErasureSetID, CurrentSymbol)) {
                Result &= WriteStripeUnit(StripeID, ErasureSetID, CurrentSymbol,
                                          0, Units2Update, pData);
            }
        }
    }
    // if none inf erased - just xor
    return Result;
}
