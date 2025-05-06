/*********************************************************
 * RAID5.h  - implementation of a RAID-5 processor
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/
#include <iostream>
#include <string.h>
#include "misc.h"
#include "RAID5.h"
#include "arithmetic.h"


using namespace std;

///initialize coding-related parameters
CRAID5Processor::CRAID5Processor(RAID5Params *P ///the configuration file
): CRAIDProcessor(P->CodeDimension + 1, 1, P, sizeof(*P)), m_pXORBuffer(0) {
    if (m_StripeUnitSize % ARITHMETIC_ALIGNMENT)
        throw Exception("Stripe size must be a multiple of #ARITHMETIC_ALIGNMENT");
};

CRAID5Processor::~CRAID5Processor() {
    AlignedFree(m_pXORBuffer);;
};


///attach to the disk array
///Prepare for multi-threaded processing
///@return true on success
bool CRAID5Processor::Attach(CDiskArray *pArray, ///the disk array
                             unsigned ConcurrentThreads
                             ///the number of concurrent processing threads that will make calls to the processor
) {
    m_pXORBuffer = AlignedMalloc(ConcurrentThreads * m_StripeUnitSize * 2);
    return CRAIDProcessor::Attach(pArray, ConcurrentThreads);
};


/**
 * If the symbol is not erased, read it from the disk. Otherwise, obtain it as a XOR of data on other disks
 */
bool CRAID5Processor::DecodeDataSymbols(unsigned long long StripeID, ///the stripe to be processed
                                        unsigned ErasureSetID, ///identifies the load balancing offset
                                        unsigned SymbolID, ///the first symbol to be processed
                                        unsigned Symbols2Decode,
                                        ///the number of symbols within this symbol to be decoded
                                        unsigned char *pDest,
                                        ///destination array. Must have size at least Subsymbols2Decode*m_StripeUnitSize
                                        size_t ThreadID ///the ID of the calling thread
) {
    bool Result = true;
    if ((GetNumOfErasures(ErasureSetID) == 0) ||
        (GetErasedPosition(ErasureSetID, 0) < (int) SymbolID) || (
            GetErasedPosition(ErasureSetID, 0) >= int(SymbolID + Symbols2Decode))) {
        //read the data as is
        for (unsigned S = SymbolID; S < SymbolID + Symbols2Decode; S++, pDest += m_StripeUnitSize) {
            Result &= ReadStripeUnit(StripeID, ErasureSetID, S, 0, 1, pDest);
        };
        return Result;
    } else {
        //read all symbols and XOR them to obtain the erased one
        unsigned S = GetErasedPosition(ErasureSetID, 0);
        unsigned i = (S == 0) ? 1 : 0;
        unsigned char *pXORBuffer = pDest + (S - SymbolID) * m_StripeUnitSize;
        unsigned char *pReadBuffer = m_pXORBuffer + (ThreadID * 2 + 1) * m_StripeUnitSize;
        //the first symbol can be loaded directly to the buffer
        Result &= ReadStripeUnit(StripeID, ErasureSetID, i, 0, 1, pXORBuffer);
        if ((i >= SymbolID) && (i < SymbolID + Symbols2Decode))
            memcpy(pDest + (i - SymbolID) * m_StripeUnitSize, pXORBuffer, m_StripeUnitSize);
        i++;
        for (; i < m_Length; i++) {
            if (i == S) continue;
            if ((i >= SymbolID) && (i < SymbolID + Symbols2Decode)) {
                //this is a payload symbol which has to be read
                unsigned char *pCurDest = pDest + (i - SymbolID) * m_StripeUnitSize;
                Result &= ReadStripeUnit(StripeID, ErasureSetID, i, 0, 1, pCurDest);
                XOR(pXORBuffer, pCurDest, m_StripeUnitSize);
            } else {
                //this symbol is needed only to compute the checksum
                Result &= ReadStripeUnit(StripeID, ErasureSetID, i, 0, 1, pReadBuffer);
                XOR(pXORBuffer, pReadBuffer, m_StripeUnitSize);
            };
        };
        //the erased symbol is now recovered
        return Result;
    };
};


/** Compute the parity symbol for the whole stripe, and write it down together with the payload data
*/
bool CRAID5Processor::EncodeStripe(unsigned long long StripeID, ///the stripe to be encoded
                                   unsigned ErasureSetID, ///identifies the load balancing offset
                                   const unsigned char *pData, ///the data to be envoced
                                   size_t ThreadID ///the ID of the calling thread
) {
    unsigned char *pXORBuffer = m_pXORBuffer + ThreadID * 2 * m_StripeUnitSize;

    bool Result = true;
    if (!IsErased(ErasureSetID, 0)) {
        Result &= WriteStripeUnit(StripeID, ErasureSetID, 0, 0, 1, pData);
    };
    memcpy(pXORBuffer, pData, m_StripeUnitSize);
    pData += m_StripeUnitSize;
    for (unsigned i = 1; i < m_Dimension; i++) {
        if (!IsErased(ErasureSetID, i)) {
            Result &= WriteStripeUnit(StripeID, ErasureSetID, i, 0, 1, pData);
        };
        XOR(pXORBuffer, pData, m_StripeUnitSize);
        pData += m_StripeUnitSize;
    };
    //write the parity symbol
    if (!IsErased(ErasureSetID, m_Dimension)) {
        Result &= WriteStripeUnit(StripeID, ErasureSetID, m_Dimension, 0, 1, pXORBuffer);
    };
    return Result;
};


/**Modify some information symbols and recompute the check sum.
 * The old values of the symbols will be read, XORed with the new ones and together to obtain delta,
 * then the old parity symbol will be read, XORed with the delta and written back
 */
bool CRAID5Processor::UpdateInformationSymbols(unsigned long long StripeID, ///the stripe to be updated,
                                               unsigned ErasureSetID, ///identifies the load balancing offset
                                               unsigned StripeUnitID, ///the first stripe unit to be updated
                                               unsigned Units2Update, ///the number of units to be updated
                                               const unsigned char *pData, ///new payload data symbols
                                               size_t ThreadID ///the ID of the calling thread
) {
    bool Result = true;
    //if the check symbol is erased, we do not need to update it
    if (IsErased(ErasureSetID, m_Dimension)) {
        //write the data as is
        for (unsigned i = 0; i < Units2Update; i++)
            Result &= WriteStripeUnit(StripeID, ErasureSetID, i + StripeUnitID, 0, 1, pData + i * m_StripeUnitSize);
    } else {
        //the parity check symbol has to be updated
        unsigned char *pXORBuffer = m_pXORBuffer + ThreadID * 2 * m_StripeUnitSize;
        unsigned char *pReadBuffer = m_pXORBuffer + (ThreadID * 2 + 1) * m_StripeUnitSize;
        unsigned S = GetErasedPosition(ErasureSetID, 0);
        if ((S >= StripeUnitID) && (S < StripeUnitID + Units2Update)) {
            //there is an erasure, and we have to update the erased symbol.
            //the updated parity check value is given by \sum_{i\not \in U} A_i +\sum_{i\in U} A_i'
            //U is the set of symbols to be updated, A_i are the old symbol values,
            //A_i' are the new symbol values
            memset(pXORBuffer, 0, m_StripeUnitSize);
            //process the symbols not to be updated
            for (unsigned i = 0; i < StripeUnitID; i++) {
                Result &= ReadStripeUnit(StripeID, ErasureSetID, i, 0, 1, pReadBuffer);
                XOR(pXORBuffer, pReadBuffer, m_StripeUnitSize);
            };
            for (unsigned i = StripeUnitID + Units2Update; i < m_Dimension; i++) {
                Result &= ReadStripeUnit(StripeID, ErasureSetID, i, 0, 1, pReadBuffer);
                XOR(pXORBuffer, pReadBuffer, m_StripeUnitSize);
            };
            //process the symbols to be updated
            for (unsigned i = 0; i < Units2Update; i++) {
                XOR(pXORBuffer, pData + i * m_StripeUnitSize, m_StripeUnitSize);
                if (S == StripeUnitID + i)
                    continue; //we cannot write to the failed disk
                else
                    Result &= WriteStripeUnit(StripeID, ErasureSetID, i + StripeUnitID, 0, 1,
                                              pData + i * m_StripeUnitSize);
            };
        } else {
            //the updated parity check value is given by S'=S +\sum_{i\in U} A_i'
            //load the old parity check symbol
            Result &= ReadStripeUnit(StripeID, ErasureSetID, m_Dimension, 0, 1, pXORBuffer);
            for (unsigned i = 0; i < Units2Update; i++) {
                XOR(pXORBuffer, pData + i * m_StripeUnitSize, m_StripeUnitSize);
                Result &= ReadStripeUnit(StripeID, ErasureSetID, i + StripeUnitID, 0, 1, pReadBuffer);
                XOR(pXORBuffer, pReadBuffer, m_StripeUnitSize);
                Result &= WriteStripeUnit(StripeID, ErasureSetID, i + StripeUnitID, 0, 1, pData + i * m_StripeUnitSize);
            };
        };
        Result &= WriteStripeUnit(StripeID, ErasureSetID, m_Dimension, 0, 1, pXORBuffer);
    };
    return Result;
};

/** Check if the sum of all codeword symbols is equal zero
* @return true on success
*/
bool CRAID5Processor::CheckCodeword(unsigned long long StripeID, ///identifies the codeword to be validated
                                    unsigned ErasureSetID, ///identifies the load balancing offset
                                    size_t ThreadID ///calling thread ID
) {
    if (GetNumOfErasures(ErasureSetID))
        //there is no way to check it for consistency
        return true;
    unsigned char *pXORBuffer = m_pXORBuffer + ThreadID * 2 * m_StripeUnitSize;
    unsigned char *pReadBuffer = m_pXORBuffer + (1 + ThreadID * 2) * m_StripeUnitSize;
    bool Result = ReadStripeUnit(StripeID, ErasureSetID, 0, 0, 1, pXORBuffer);
    for (unsigned i = 1; i < m_Length; i++) {
        Result &= ReadStripeUnit(StripeID, ErasureSetID, i, 0, 1, pReadBuffer);
        XOR(pXORBuffer, pReadBuffer, m_StripeUnitSize);
    };
    if (!Result)
        return false;
    else {
        unsigned char S = 0;
        for (unsigned i = 0; i < m_StripeUnitSize; i++)
            S |= pXORBuffer[i];
        return S == 0;
    };
};
