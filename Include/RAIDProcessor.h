/*********************************************************
 * RAIDProcessor.h  - header file for a generic wrapper for various RAID processing algorithms
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/
#ifndef RAIDPROCESSOR_H
#define RAIDPROCESSOR_H


#include <stdlib.h>


class CDiskArray;


///the basic configuration set for a RAID code
///avoid padding of struct members
#pragma pack(push)
#pragma pack(1)
struct RAIDParams {
    ///type of the code
    int Type;
    ///dimension of the code
    unsigned CodeDimension;
    ///size of one stripe unit
    unsigned StripeUnitSize;
    ///number of independent arrays
    unsigned InterleavingOrder;

    RAIDParams(int type, unsigned Dimension, unsigned interleavingOrder, unsigned stripeUnitSize): Type(type),
        CodeDimension(Dimension),
        StripeUnitSize(stripeUnitSize), InterleavingOrder(interleavingOrder) {
    };
};
#pragma pack(pop)


///this is a base class for all RAID data processing algorithms
///It implements also cyclic load balancing across the drives
///each derived class must be able to support a given number of parallel calls
class CRAIDProcessor {
    ///the full configuration record
    RAIDParams *m_pParams;
    ///size of the configuration record
    unsigned m_ConfigSize;
    ///provides interface for reading and writing data on disks
    CDiskArray *m_pArray;
    ///number of offline disks in each of the subarrays
    unsigned *m_pNumOfOfflineDisks;
    ///IDs of offline disks for each of the subarrays
    unsigned **m_ppOfflineDisks;
    ///the temporary buffer for data update
    unsigned char *m_pUpdateBuffer;

protected:
    ///length of the array code
    unsigned m_Length;
    ///the number of stripe units constituting a single codeword symbol
    unsigned m_StripeUnitsPerSymbol;
    ///the number of information (payload) symbols per codeword (stripe)
    unsigned m_Dimension;
    ///the number of arrays operating jointly
    unsigned m_InterleavingOrder;
    ///the number of bytes in each stripe units
    unsigned m_StripeUnitSize;

    ///Read from disks a contiguous set of stripe units corresponding to the same symbol
    ///In other words, read a number of subsymbols corresponding to some symbol
    ///This function will not check if the requested range crosses the stripe boundary
    ///This function provides also load balancing by cyclically mapping symbols onto disks (offset is given by ErasureSetID)
    ///@return true on success
    bool ReadStripeUnit(unsigned long long StripeID, ///identifies the codeword (stripe)
                        unsigned ErasureSetID, ///identifies the load balancing offset
                        unsigned SymbolID, ///identifies the disk to be accessed
                        unsigned StripeUnitID, ///identifies the first subsymbol to be read
                        unsigned Units2Read, ///number of stripe units to be loaded
                        void *pDest ///the destination buffer. Must have size  Units2Read*m_StripeUnitSize
    );

    ///Write to disks a contiguous set of stripe units corresponding to the same symbol
    ///In other words, read a number of subsymbols corresponding to some symbol
    ///This function will not check if the requested range crosses the stripe boundary
    ///This function provides also load balancing by cyclically mapping symbols onto disks (offset is given by ErasureSetID)
    ///@return true on success
    bool WriteStripeUnit(unsigned long long StripeID, ///identifies the codeword (stripe)
                         unsigned ErasureSetID, ///identifies the load balancing offset
                         unsigned SymbolID, ///identifies the disk to be accessed
                         unsigned StripeUnitID, ///identifies the first subsymbol to be read
                         unsigned Units2Write, ///number of stripe units to be loaded
                         const void *pSrc ///the data to be written (Units2Read*m_StripeUnitSize bytes)
    );

    ///Check if it is possible to correct a given combination of erasures
    ///If yes, the method should initialize the internal data structures
    ///and be ready to do the actual erasure correction. This combination of erasures
    /// is uniquely identified by ErasureID
    ///@return true if the specified combination of erasures is correctable
    virtual bool IsCorrectable(
        unsigned ErasureSetID ///identifies the erasure combination. This will not exceed m_Length-1
    ) =0;

    ///get the total number of erasures
    unsigned GetNumOfErasures(unsigned ErasureSetID) const {
        return m_pNumOfOfflineDisks[ErasureSetID / m_Length];
    };
    ///@return the i-th erased symbol, or -1 if it does not exist
    int GetErasedPosition(unsigned ErasureSetID, ///the erasure combination
                          unsigned i ///ID of the erased symbol
    ) const {
        unsigned SubarrayID = ErasureSetID / m_Length;
        ErasureSetID = ErasureSetID % m_Length;
        if (i >= m_pNumOfOfflineDisks[SubarrayID])
            return -1;
        //cyclic load balancing for disks
        int Res = m_ppOfflineDisks[SubarrayID][i] - ErasureSetID;
        if (Res < 0)
            Res += m_Length;
        return Res;
    };

    ///@return true if the i-th symbol is erased
    bool IsErased(unsigned ErasureSetID, ///the erasure combination (identifies the load balancing offset)
                  unsigned i) const;

    ///decode a number of payload subsymbols from a given symbol
    ///@return true on success
    virtual bool DecodeDataSubsymbols(unsigned long long StripeID, ///the stripe to be processed
                                      unsigned ErasureSetID, ///identifies the load balancing offset
                                      unsigned SymbolID, ///the symbol to be processed
                                      unsigned SubsymbolID, ///the first subsymbol to be processed
                                      unsigned Subsymbols2Decode,
                                      ///the number of subsymbols within this symbol to be decoded
                                      unsigned char *pDest,
                                      ///destination array. Must have size at least Subsymbols2Decode*m_StripeUnitSize
                                      size_t ThreadID ///the ID of the calling thread
    ) =0;

    ///decode a number of payload subsymbols from a given symbol
    ///@return true on success
    virtual bool DecodeDataSymbols(unsigned long long StripeID, ///the stripe to be processed
                                   unsigned ErasureSetID, ///identifies the load balancing offset
                                   unsigned SymbolID, ///the first symbol to be processed
                                   unsigned Symbols2Decode, ///the number of subsymbols within this symbol to be decoded
                                   unsigned char *pDest,
                                   ///destination array. Must have size at least Subsymbols2Decode*m_StripeUnitSize
                                   size_t ThreadID ///the ID of the calling thread
    ) =0;

    ///make a decision about the optimal data encoding method
    ///@return true if one should fetch all the data from disk, update it and write back. Otherwise, only check symbols will be updated
    virtual bool GetEncodingStrategy(unsigned ErasureSetID, ///identifies the load balancing offset
                                     unsigned StripeUnitID, ///the first stripe unit to be updated
                                     unsigned Subsymbols2Encode ///the number of subsymbols to be encoded
    ) {
        return (Subsymbols2Encode > 2 * m_Dimension * m_StripeUnitsPerSymbol / 3);
    };

    ///encode and write the whole CRAIDProcessorstripe
    ///@return true on success
    virtual bool EncodeStripe(unsigned long long StripeID, ///the stripe to be encoded
                              unsigned ErasureSetID, ///identifies the load balancing offset
                              const unsigned char *pData, ///the data to be envoced
                              size_t ThreadID ///the ID of the calling thread
    ) =0;

    ///update some information symbols and the corresponding check symbols
    ///@return true on success
    virtual bool UpdateInformationSymbols(unsigned long long StripeID, ///the stripe to be updated,
                                          unsigned ErasureSetID, ///identifies the load balancing offset
                                          unsigned StripeUnitID, ///the first stripe unit to be updated
                                          unsigned Units2Update, ///the number of units to be updated
                                          const unsigned char *pData, ///new payload data symbols
                                          size_t ThreadID ///the ID of the calling thread
    ) =0;

    ///check if the codeword is consistent
    virtual bool CheckCodeword(unsigned long long StripeID, ///the stripe to be checked
                               unsigned ErasureSetID, ///identifies the load balancing offset
                               size_t ThreadID ///identifies the calling thread
    ) =0;

public:
    ///initialize coding-related parameters
    CRAIDProcessor(unsigned Length, ///the length of the array code
                   unsigned StripeUnitsPerSymbol, ///number of subsymbols per codeword symbol
                   RAIDParams *pParams, ///Configuration of the RAID array
                   unsigned ConfigSize ///size of the configuration entry
    );

    virtual ~CRAIDProcessor();

    ///@return the length of the code
    unsigned GetCodeLength() const {
        return m_Length;
    };
    ///@return stripe unit size
    unsigned GetStripeUnitSize() const {
        return m_StripeUnitSize;
    };
    ///@return the number of stripe units per symbol
    unsigned GetStripeUnitsPerSymbol() const {
        return m_StripeUnitsPerSymbol;
    };
    ///@return the number of symbols per codeword
    unsigned GetDimension() const {
        return m_Dimension;
    };
    ///@return the number of independent arrays operating jointly
    unsigned GetInterleavingOrder() const {
        return m_InterleavingOrder;
    };
    ///get the full configuration record of the code
    ///@return record size
    unsigned GetConfiguration(const void *&pData) {
        pData = m_pParams;
        return m_ConfigSize;
    };
    ///@return RAID type
    int GetType() const {
        return m_pParams->Type;
    };


    ///attach to the disk array
    ///Prepare for multi-threaded processing
    ///The derived classes must override this method, and the overridden one
    ///must make a call to the one in the parent class AFTER all general initialization has been done.
    ///This method will make a call to ResetErasures() method
    ///@return true on success
    virtual bool Attach(CDiskArray *pArray, ///the disk array
                        unsigned ConcurrentThreads
                        ///the number of concurrent processing threads that will make calls to the processor
    );

    ///reset the erasure correction engine
    /// this will be called if the set of failed disks changes
    /// The derived class must first call the method in the parent one
    virtual void ResetErasures();

    ///check if we have sufficient amount of online disks in the array
    ///so that the data can be recovered
    ///@return true if read and write access to the data is possible
    bool IsMountable();

    ///get some payload stripe units from a given stripe
    ///@return true on success
    bool ReadData(unsigned long long StripeID, ///the stripe to be read
                  unsigned StripeUnitID, ///the first payload  stripe unit to read
                  unsigned SubarrayID, ///identifies the subarray to be used
                  unsigned NumOfUnits, ///the number of units to read
                  unsigned char *pDest, ///destination buffer. Must have size at least NumOfUnits*m_StripeUnitSize
                  size_t ThreadID ///calling thread ID
    );

    ///write some payload stripe units to a given stripe
    ///@return true on success
    bool WriteData(unsigned long long StripeID, ///the stripe to be written
                   unsigned StripeUnitID, ///the first payload stripe unit to write
                   unsigned SubarrayID, ///identifies the subarray to be used
                   unsigned NumOfUnits, ///the number of units to write
                   const unsigned char *pSrc, ///source data . Must have size at least NumOfUnits*m_StripeUnitSize
                   size_t ThreadID ///calling thread ID
    );

    ///make sure that the codeword is a legal one
    ///@return true on success
    bool VerifyStripe(unsigned long long StripeID, ///identifies the codeword to be validated
                      size_t ThreadID ///calling thread ID
    ) {
        return CheckCodeword(StripeID, StripeID % m_Length, ThreadID);
    };
};

#include "RAIDconfig.h"


#endif
