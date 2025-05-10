/*********************************************************
 * RDP.h  - implementation of a RAID-5 processor
 * Author: N. Zavialov 408627@niuitmo.ru
 * ********************************************************/
#ifndef RDP_H
#define RDP_H
#include "RAIDProcessor.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <ostream>
#include <ranges>
#include <unordered_map>
#include <vector>

class CRDPProcessor final : public CRAIDProcessor {
private:
  long long m_PrimeNumber;
  unsigned char* m_pRDPXORBuffer;
  long long m_Shortering;

  static constexpr size_t MAX_ERASURES = 3;
  static constexpr size_t BLOCKS_PER_BUFFER = MAX_ERASURES + 2;
  size_t BLOCK_SIZE;
  size_t BUFFER_SIZE;

  std::vector<unsigned> fst, snd;
  std::vector<std::vector<unsigned>> prod;

      /**
       * Helper function to read symbols that are not erased and calculate syndroms
       * @param StripeID the stripe to be processed
       * @param ErasureSetID identifies the load balancing offset
       * @param SymbolID the symbol to be processed
       * @param Symbols2Decode the number of subsymbols within this symbol to be decoded
       * @param pDest destination array. Must have size at least Subsymbols2Decode*m_StripeUnitSize
       * @param ThreadID the ID of the calling thread
       * @param pXORBuffer pointers to buffer where syndromes will be stored(nullptr for unused syndromes)
       * @param CalcMissingDiag optional argument for calculating diagonal syndrome of the (m_PrimeNumber-1) diagonal
       * @return true on a success, false otherwise
       */
      bool
      ReadAndCalcSyndroms(
          unsigned long long StripeID,
          unsigned ErasureSetID,
          unsigned SymbolID,
          unsigned Symbols2Decode,
          unsigned char* pDest,
          size_t ThreadID,
          std::array<unsigned char*, 3> pXORBuffer,
          bool CalcMissingDiag
      );

  /**
   * Helper function to restore three erasures
   * @param StripeID the stripe to be processed
   * @param ErasureSetID identifies the load balancing offset
   * @param SymbolID the symbol to be processed
   * @param Symbols2Decode the number of subsymbols within this symbol to be decoded
   * @param pDest destination array. Must have size at least Subsymbols2Decode*m_StripeUnitSize
   * @param ThreadID the ID of the calling thread
   * @param ErasurePosition optional argument if we already know what symbol we want to restore
   * @return true on a success, false otherwise
   */
  bool DecodeOneErasure(
      unsigned long long StripeID,
      unsigned ErasureSetID,
      unsigned SymbolID,
      unsigned Symbols2Decode,
      unsigned char* pDest,
      size_t ThreadID,
      long long ErasurePosition
  );

  /**
   * Helper function to restore three erasures
   * @param StripeID the stripe to be processed
   * @param ErasureSetID identifies the load balancing offset
   * @param SymbolID the symbol to be processed
   * @param Symbols2Decode the number of subsymbols within this symbol to be decoded
   * @param pDest destination array. Must have size at least Subsymbols2Decode*m_StripeUnitSize
   * @param ThreadID the ID of the calling thread
   * @param diag_step diagonal that will be used to restore erasures
   * @param ErasurePosition optional argument if we already know what symbols we want to restore
   * @return true on a success, false otherwise
   */
  bool DecodeTwoErasures(
      unsigned long long StripeID,
      unsigned ErasureSetID,
      unsigned SymbolID,
      unsigned Symbols2Decode,
      unsigned char* pDest,
      size_t ThreadID,
      long long diag_step,
      std::array<int, 2> ErasurePosition
  );

  /**
   * Helper function to restore three erasures
   * @param StripeID the stripe to be processed
   * @param ErasureSetID identifies the load balancing offset
   * @param SymbolID the symbol to be processed
   * @param Symbols2Decode the number of subsymbols within this symbol to be decoded
   * @param pDest destination array. Must have size at least Subsymbols2Decode*m_StripeUnitSize
   * @param ThreadID the ID of the calling thread
   * @return true on a success, false otherwise
   */
  bool DecodeThreeErasures(
      unsigned long long StripeID,
      unsigned ErasureSetID,
      unsigned SymbolID,
      unsigned Symbols2Decode,
      unsigned char* pDest,
      size_t ThreadID
  );

  /**
   * Restore two information symbols
   * @param diag_step describes which diagonal to use(offset of the diagonal)
   * @param diag_symbol ID of symbols recovered by diagonal syndrome
   * @param horizontal_symbol ID of symbols recovered by horizontal syndrome
   * @param RecoverResults data to store recovered symbols
   * @param pXORBuffer buffer with calculated syndromes
   * @param diag current diagonal to started restoring
   * @param step difference between two symbols IDs
   * @param diag_restore which of the RecoverResults will be recovered by diagonal syndrome
   * @param horizontal_restore which of the RecoverResults will be recovered by horizontal syndrome
   * @param restore_missing_diag should we recalc syndrome of unused diagonal after restoring symbols
   */
  void iterative_restore(
      long long diag_step,
      long long diag_symbol,
      long long horizontal_symbol,
      std::array<unsigned char*, 2> RecoverResults,
      std ::array<unsigned char*, 3> pXORBuffer,
      long long diag,
      long long step,
      unsigned diag_restore,
      unsigned horizontal_restore,
      bool restore_missing_diag
  ) const;

protected:
  bool IsCorrectable(unsigned ErasureSetID) override;

  bool DecodeDataSubsymbols(
      unsigned long long StripeID,
      unsigned ErasureSetID,
      unsigned SymbolID,
      unsigned SubsymbolID,
      unsigned Subsymbols2Decode,
      unsigned char* pDest,
      size_t ThreadID
  ) override;

  bool DecodeDataSymbols(
      unsigned long long StripeID,
      unsigned ErasureSetID,
      unsigned SymbolID,
      unsigned Symbols2Decode,
      unsigned char* pDest,
      size_t ThreadID
  ) override;

  bool EncodeStripe(unsigned long long StripeID, unsigned ErasureSetID, const unsigned char* pData, size_t ThreadID)
      override;

  bool UpdateInformationSymbols(
      unsigned long long StripeID,
      unsigned ErasureSetID,
      unsigned StripeUnitID,
      unsigned Units2Update,
      const unsigned char* pData,
      size_t ThreadID
  ) override;

  bool CheckCodeword(unsigned long long StripeID, unsigned ErasureSetID, size_t ThreadID) override;

  bool GetEncodingStrategy(unsigned ErasureSetID, unsigned StripeUnitID, unsigned Subsymbols2Encode) override;

public:
  /**
   * Constructor for RDP Processor
   * @param pParams parameters of the code
   */
  CRDPProcessor(RDPParams* pParams);

  ~CRDPProcessor() override;

  bool Attach(CDiskArray* pArray, unsigned ConcurrentThreads) override;

  void ResetErasures() override;
};

#endif // RDP_H
