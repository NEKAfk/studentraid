/*********************************************************
 * RAID5.h  - header file for a RAID-5 processor
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/
#ifndef RAID5_H
#define RAID5_H

#include "RAIDProcessor.h"

/*///there is nothing special to configure for RAID5
struct RAID5Params:public RAIDParams
{
  RAID5Params(unsigned Dimension,unsigned stripeUnitSize):RAIDParams(rtRAID5,Dimension,stripeUnitSize)
  {

  };
};*/

class CRAID5Processor : public CRAIDProcessor {
  /// the buffer used for parity computation
  unsigned char* m_pXORBuffer;

protected:
  /// Check if it is possible to correct a given combination of erasures
  /// If yes, the method should initialize the internal data structures
  /// and be ready to do the actual erasure correction. This combination of erasures
  ///  is uniquely identified by ErasureID
  ///@return true if the specified combination of erasures is correctable
  virtual bool IsCorrectable(
      unsigned ErasureSetID /// identifies the erasure combination. This will not exceed m_Length-1
  ) {
    return GetNumOfErasures(ErasureSetID) <= 1;
  };

  /// reset the erasure correction engine
  ///  this will be called if the set of failed disks changes
  virtual void ResetErasures() {
    CRAIDProcessor::ResetErasures();
  };

  /// This is a stub which should never be called
  virtual bool DecodeDataSubsymbols(
      unsigned long long StripeID, /// the stripe to be processed
      unsigned ErasureSetID, /// identifies the load balancing offset
      unsigned SymbolID, /// the symbol to be processed
      unsigned SubsymbolID, /// the first subsymbol to be processed
      unsigned Subsymbols2Decode,
      /// the number of subsymbols within this symbol to be decoded
      unsigned char* pDest,
      /// destination array. Must have size at least Subsymbols2Decode*m_StripeUnitSize
      size_t ThreadID /// the ID of the calling thread
  ) {
    return false;
  };

  /// decode a number of payload subsymbols from a given symbol
  ///@return true on success
  virtual bool DecodeDataSymbols(
      unsigned long long StripeID, /// the stripe to be processed
      unsigned ErasureSetID, /// identifies the load balancing offset
      unsigned SymbolID, /// the first symbol to be processed
      unsigned Symbols2Decode, /// the number of subsymbols within this symbol to be decoded
      unsigned char* pDest,
      /// destination array. Must have size at least Subsymbols2Decode*m_StripeUnitSize
      size_t ThreadID /// the ID of the calling thread
  );

  /// encode and write the whole stripe
  ///@return true on success
  virtual bool EncodeStripe(
      unsigned long long StripeID, /// the stripe to be encoded
      unsigned ErasureSetID, /// identifies the load balancing offset
      const unsigned char* pData, /// the data to be envoced
      size_t ThreadID /// the ID of the calling thread
  );

  /// update some information symbols and the corresponding check symbols
  ///@return true on success
  virtual bool UpdateInformationSymbols(
      unsigned long long StripeID, /// the stripe to be updated,
      unsigned ErasureSetID, /// identifies the load balancing offset
      unsigned StripeUnitID, /// the first stripe unit to be updated
      unsigned Units2Update, /// the number of units to be updated
      const unsigned char* pData, /// new payload data symbols
      size_t ThreadID /// the ID of the calling thread
  );

  /// make sure that the codeword is a legal one
  ///@return true on success
  bool CheckCodeword(
      unsigned long long StripeID, /// identifies the codeword to be validated
      unsigned ErasureSetID, /// identifies the load balancing offset
      size_t ThreadID /// calling thread ID
  );

public:
  /// initialize coding-related parameters
  CRAID5Processor(RAID5Params* P /// the configuration file
  );

  ~CRAID5Processor();

  /// attach to the disk array
  /// Prepare for multi-threaded processing
  ///@return true on success
  virtual bool Attach(
      CDiskArray* pArray, /// the disk array
      unsigned ConcurrentThreads
      /// the number of concurrent processing threads that will make calls to the processor
  );
};
#endif
