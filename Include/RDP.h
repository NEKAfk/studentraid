#ifndef RDP_H
#define RDP_H
#include "RAIDProcessor.h"

#include <array>

class CRDPProcessor final : public CRAIDProcessor {
private:
  long long m_PrimeNumber;
  unsigned char* m_pRDPXORBuffer;

  static constexpr size_t MAX_ERASURES = 3;
  static constexpr size_t BLOCKS_PER_BUFFER = MAX_ERASURES + 2;
  size_t BLOCK_SIZE;
  size_t BUFFER_SIZE;

  bool DecodeOneErasure(
      unsigned long long StripeID,
      unsigned ErasureSetID,
      unsigned SymbolID,
      unsigned Symbols2Decode,
      unsigned char* pDest,
      size_t ThreadID,
      long long ErasurePosition
  );

protected:
  bool IsCorrectable(unsigned ErasureSetID) override {
    return GetNumOfErasures(ErasureSetID) <= MAX_ERASURES;
  }

  bool DecodeDataSubsymbols(
      unsigned long long StripeID,
      unsigned ErasureSetID,
      unsigned SymbolID,
      unsigned SubsymbolID,
      unsigned Subsymbols2Decode,
      unsigned char* pDest,
      size_t ThreadID
  ) override;

  bool ReadAndCalcSyndroms(
      unsigned long long StripeID,
      unsigned ErasureSetID,
      unsigned SymbolID,
      unsigned Symbols2Decode,
      unsigned char* pDest,
      size_t ThreadID,
      std::array<unsigned char*, 3> pXORBuffer,
      bool CalcMissingDiag
  );

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

  bool DecodeThreeErasures(
      unsigned long long StripeID,
      unsigned ErasureSetID,
      unsigned SymbolID,
      unsigned Symbols2Decode,
      unsigned char* pDest,
      size_t ThreadID
  );

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
  CRDPProcessor(RDPParams* pParams);

  ~CRDPProcessor() override;

  bool Attach(CDiskArray* pArray, unsigned ConcurrentThreads) override;

  void ResetErasures() override;
};

#endif // RDP_H
