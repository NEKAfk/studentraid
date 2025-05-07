/*********************************************************
 * RAIDProcessor.cpp  - implementation of the generic wrapper for various RAID processing algorithms
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/
#include "RAIDProcessor.h"

#include "RAIDconfig.h"
#include "array.h"
#include "misc.h"

#include <stdlib.h>
#include <string.h>

#include <iostream>

using namespace std;

/// initialize coding-related parameters
CRAIDProcessor::CRAIDProcessor(
    unsigned Length, /// the length of the array code
    unsigned StripeUnitsPerSymbol, /// number of subsymbols per codeword symbol
    RAIDParams* pParams, /// Configuration of the RAID array
    unsigned ConfigSize /// size of the configuration entry
)
    : m_pParams(pParams)
    , m_ConfigSize(ConfigSize)
    , m_Length(Length)
    , m_Dimension(pParams->CodeDimension)
    , m_StripeUnitSize(pParams->StripeUnitSize)
    , m_StripeUnitsPerSymbol(StripeUnitsPerSymbol)
    , m_pArray(0)
    , m_pNumOfOfflineDisks(0)
    , m_ppOfflineDisks(0)
    , m_pUpdateBuffer(0)
    , m_InterleavingOrder(pParams->InterleavingOrder) {
  if (!m_Dimension || !m_StripeUnitSize || !m_StripeUnitsPerSymbol || !m_InterleavingOrder) {
    throw Exception(
        "Invalid initialization for RAID processor:\n"
        "Dimension=%d, StripeUnitSize=%d, StripeUnitsPersymbol=%d, InterleavingOrder=%d",
        m_Dimension,
        m_StripeUnitSize,
        m_StripeUnitsPerSymbol,
        m_InterleavingOrder
    );
  }
  m_pNumOfOfflineDisks = new unsigned[m_InterleavingOrder];
  m_ppOfflineDisks = new unsigned*[m_InterleavingOrder];
  memset(m_ppOfflineDisks, 0, sizeof(unsigned*) * m_InterleavingOrder);
};

CRAIDProcessor::~CRAIDProcessor() {
  for (unsigned i = 0; i < m_InterleavingOrder; i++) {
    delete[] m_ppOfflineDisks[i];
  }
  delete[] m_ppOfflineDisks;
  delete[] m_pNumOfOfflineDisks;
  delete[] m_pUpdateBuffer;
  delete m_pParams;
};

/**Attach to the disk array,
 * allocate memory for data encoding,
 * inspect the disks and find those not being online
 */
bool CRAIDProcessor::Attach(
    CDiskArray* pArray, /// the disk array
    unsigned ConcurrentThreads
    /// the number of concurrent processing threads that will make calls to the processor
) {
  m_pArray = pArray;
  m_pUpdateBuffer = new unsigned char[ConcurrentThreads * m_Dimension * m_StripeUnitsPerSymbol * m_StripeUnitSize];

  ResetErasures();
  return true;
};

/** Mark the non-online disks as erased
 */
void CRAIDProcessor::ResetErasures() {
  // get the offline disks for each subarray
  for (unsigned j = 0; j < m_InterleavingOrder; j++) {
    m_pNumOfOfflineDisks[j] = 0;
    for (unsigned i = 0; i < m_Length; i++) {
      if (m_pArray->m_pDisks[j * m_Length + i].GetDiskState() != dsOnline) {
        m_pNumOfOfflineDisks[j]++;
      }
    };
    // enumerate the offline disks
    if (m_ppOfflineDisks[j]) {
      delete[] m_ppOfflineDisks[j];
    }
    if (m_pNumOfOfflineDisks[j]) {
      m_ppOfflineDisks[j] = new unsigned[m_pNumOfOfflineDisks[j]];
      unsigned S = 0;
      for (unsigned i = 0; i < m_Length; i++) {
        if (m_pArray->m_pDisks[j * m_Length + i].GetDiskState() != dsOnline) {
          m_ppOfflineDisks[j][S++] = i;
        }
      };
    } else {
      m_ppOfflineDisks[j] = 0;
    }
  };
};

/** Check if the corresponding disk is not online.
 * */
bool CRAIDProcessor::IsErased(
    unsigned ErasureSetID, /// the erasure combination (identifies the load balancing offset)
    unsigned i
) const {
  unsigned SubarrayID = ErasureSetID / m_Length;
  i = (i + ErasureSetID) % m_Length;
  return m_pArray->m_pDisks[i + SubarrayID * m_Length].GetDiskState() != dsOnline;
}

/** Check if all possible erasure patterns are correctable
 *
 *
 * */
bool CRAIDProcessor::IsMountable() {
  // make sure that all cyclic shifts of the erasure pattern are correctable
  bool Result = true;
  for (unsigned i = 0; i < m_pArray->m_NumOfDisks; i++) {
    // prepare to correct all erasure patterns
    Result &= IsCorrectable(i);
  }
  return Result;
};

/**Read a number of stripe units from the disk. Implements cyclic mapping of codeword symbols onto the disks.
 * The offset is given by ErasureSetID
 *
 * */
bool CRAIDProcessor::ReadStripeUnit(
    unsigned long long StripeID, /// identifies the codeword (stripe)
    unsigned ErasureSetID, /// identifies the load balancing offset
    unsigned SymbolID, /// identifies the disk to be accessed
    unsigned StripeUnitID, /// identifies the first subsymbol to be read
    unsigned Units2Read, /// number of stripe units to be loaded
    void* pDest /// the destination buffer. Must have size  Units2Read*m_StripeUnitSize
) {
  unsigned SubarrayID = ErasureSetID / m_Length;
  SymbolID = (SymbolID + ErasureSetID) % m_Length;
  return m_pArray->m_pDisks[SymbolID + SubarrayID * m_Length]
      .ReadData(StripeID * m_StripeUnitsPerSymbol + StripeUnitID, Units2Read, pDest);
};

/**Write a number of stripe units to the disk. Implements cyclic mapping of codeword symbols onto the disks.
 * The offset is given by ErasureSetID
 *
 * */
bool CRAIDProcessor::WriteStripeUnit(
    unsigned long long StripeID, /// identifies the codeword (stripe)
    unsigned ErasureSetID, /// identifies the load balancing offset
    unsigned SymbolID, /// identifies the disk to be accessed
    unsigned StripeUnitID, /// identifies the first subsymbol to be read
    unsigned Units2Write, /// number of stripe units to be loaded
    const void* pSrc /// the data to be written (Units2Read*m_StripeUnitSize bytes)
) {
  unsigned SubarrayID = ErasureSetID / m_Length;
  SymbolID = (SymbolID + ErasureSetID) % m_Length;
  return m_pArray->m_pDisks[SymbolID + SubarrayID * m_Length]
      .WriteData(StripeID * m_StripeUnitsPerSymbol + StripeUnitID, Units2Write, pSrc);
};

/** Translates the read request into a number of decoder calls.
 * This method essentially splits the Read call into a number of Decode calls
 *
 *
 */
bool CRAIDProcessor::ReadData(
    unsigned long long StripeID, /// the stripe to be read
    unsigned StripeUnitID, /// the first payload  stripe unit to read
    unsigned SubarrayID, /// identifies the subarray to be used
    unsigned NumOfUnits, /// the number of units to read
    unsigned char* pDest,
    /// destination buffer. Must have size at least NumOfUnits*m_StripeUnitSize
    size_t ThreadID /// calling thread ID
) {
  unsigned FirstSymbolID = StripeUnitID / m_StripeUnitsPerSymbol;
  unsigned FirstSymbolOffset = StripeUnitID % m_StripeUnitsPerSymbol;
  unsigned ErasureSetID = (StripeID % m_Length) + SubarrayID * m_Length;
  bool Result = true;
  if (FirstSymbolOffset) {
    // incomplete symbol read request
    unsigned Units2Read = min(m_StripeUnitsPerSymbol - FirstSymbolOffset, NumOfUnits);
    if (!Units2Read) {
      return Result;
    }
    Result &=
        DecodeDataSubsymbols(StripeID, ErasureSetID, FirstSymbolID, FirstSymbolOffset, Units2Read, pDest, ThreadID);
    pDest += Units2Read * m_StripeUnitSize;
    if (NumOfUnits < Units2Read) {
      return false;
    }
    NumOfUnits -= Units2Read;
    FirstSymbolID++;
  };
  // read the entire symbols
  unsigned Symbols2Decode = NumOfUnits / m_StripeUnitsPerSymbol;
  if (Symbols2Decode) {
    Result &= DecodeDataSymbols(StripeID, ErasureSetID, FirstSymbolID, Symbols2Decode, pDest, ThreadID);
  }
  if (NumOfUnits < Symbols2Decode * m_StripeUnitsPerSymbol) {
    return false;
  }
  NumOfUnits -= Symbols2Decode * m_StripeUnitsPerSymbol;
  pDest += Symbols2Decode * m_StripeUnitsPerSymbol * m_StripeUnitSize;
  FirstSymbolID += Symbols2Decode;
  // check if we still have something to decode
  if (NumOfUnits > 0) {
    // incomplete symbol read request
    /*if (NumOfUnits>m_Length) {
        return false;
    }*/
    if (!NumOfUnits) {
      return Result;
    }
    Result &= DecodeDataSubsymbols(StripeID, ErasureSetID, FirstSymbolID, 0, NumOfUnits, pDest, ThreadID);
  };
  return Result;
};

/**Translate write call into a number of Encode calls
 * The encoding strategy is determined by the GetEncodingStrategy function. If needed,
 * this method will get all non-affected the data from the disk and re-encode it
 * */
bool CRAIDProcessor::WriteData(
    unsigned long long StripeID, /// the stripe to be written
    unsigned StripeUnitID, /// the first payload stripe unit to write
    unsigned SubarrayID, /// identifies the subarray to be used
    unsigned NumOfUnits, /// the number of units to write
    const unsigned char* pSrc,
    /// source data . Must have size at least NumOfUnits*m_StripeUnitSize
    size_t ThreadID /// calling thread ID
) {
  unsigned ErasureSetID = (StripeID % m_Length) + SubarrayID * m_Length;
  bool Result = true;
  if (GetEncodingStrategy(ErasureSetID, StripeUnitID, NumOfUnits)) {
    if (NumOfUnits == m_Dimension * m_StripeUnitsPerSymbol) {
      Result &= EncodeStripe(StripeID, ErasureSetID, pSrc, ThreadID);
    } else {
      unsigned char* pBuffer = m_pUpdateBuffer + ThreadID * m_Dimension * m_StripeUnitsPerSymbol * m_StripeUnitSize;
      if (StripeUnitID) {
        // fetch the data residing before the new data
        Result &= ReadData(StripeID, 0, SubarrayID, StripeUnitID, pBuffer, ThreadID);
      };
      memcpy(pBuffer + StripeUnitID * m_StripeUnitSize, pSrc, NumOfUnits * m_StripeUnitSize);
      int TrailingUnits = m_Dimension * m_StripeUnitsPerSymbol - (StripeUnitID + NumOfUnits);
      if (TrailingUnits > 0) {
        // fetch the data residing after the new data
        Result &= ReadData(
            StripeID,
            StripeUnitID + NumOfUnits,
            SubarrayID,
            TrailingUnits,
            pBuffer + (StripeUnitID + NumOfUnits) * m_StripeUnitSize,
            ThreadID
        );
      };
      Result &= EncodeStripe(StripeID, ErasureSetID, pBuffer, ThreadID);
    };
    return Result;
  } else {
    // update selected symbols
    bool Res = UpdateInformationSymbols(StripeID, ErasureSetID, StripeUnitID, NumOfUnits, pSrc, ThreadID);
    return Res;
  };
}
