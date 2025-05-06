/*********************************************************
 * RSRAID.h  - header file for a Reed-Solomon based RAID
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/
#ifndef RSRAID_H
#define RSRAID_H

#include "config.h"
#include "RAIDProcessor.h"
#include "arithmetic.h"

#ifdef CYCLOTOMIC_FFT
#include "CyclotomicFFT.h"
#endif

///implements a generic Reed-Solomon based RAID 
///some symbols of the codeword are designated as 
///information (m_pInfSymbols) and check (m_pCheckSymbols) ones
///Those not listed in both arrays are information symbols, 
///which are always set to zero and not stored on any disk
///In all cases, the underlying RS code has length 255
class CRSProcessor : public CRAIDProcessor {
    ///the number of check symbols
	///this field is redundant is needed to
	///improve performance
    unsigned m_Redundancy;
    ///the indices of information symbols sorted in ascending order
    int *m_pInfSymbols;
    ///the indices of check symbols sorted in ascending order
    int *m_pCheckSymbols;

    //true if cyclotomic processing is used
    bool m_CyclotomicProcessing;
    ///temporary array for cyclotomic processing
    GFValue *m_pCyclotomicTemp;
    //true if the check symbol locators are optimized ones
    bool m_OptimizedCheckLocators;

    ///the check symbol locator polynomial
    GFValue *m_pCheckLocator;
    ///values of \alpha^{1-b}/\Lambda'(1/X_i) (needed by Forney algorithm)
    ///where X_i are locators of check symbols
    int *m_pCheckLocatorsPrime;
    ///syndromes for each stripe unit
    GFValue *m_pSyndromes;
    ///the erasure evaluator polynomial
    GFValue *m_pErasureEvaluator;
    ///the erasure locator polynomial for each erasure configuration
    GFValue *m_pErasureLocators;
    ///values of \alpha^{1-b}/\Lambda'(1/X_i) (needed by Forney algorithm)
    ///where X_i are locators of erased symbols
    int *m_pErasureLocatorsPrime;
    ///buffer for fetching the codeword symbols
    GFValue *m_pSymbols;
    ///pointers to the fetched symbols
    const GFValue **m_ppSymbols;

protected:
    ///attach to the disk array
    ///Prepare for multi-threaded processing
    ///@return true on success
    virtual bool Attach(CDiskArray *pArray, ///the disk array
                        unsigned ConcurrentThreads
                        ///the number of concurrent processing threads that will make calls to the processor
    );

    ///reset the erasure correction engine
    /// this will be called if the set of failed disks changes
    virtual void ResetErasures();

    ///Check if it is possible to correct a given combination of erasures
    ///If yes, the method should initialize the internal data structures
    ///and be ready to do the actual erasure correction. This combination of erasures
    /// is uniquely identified by ErasureID
    ///@return true if the specified combination of erasures is correctable
    virtual bool IsCorrectable(
        unsigned ErasureSetID ///identifies the erasure combination. This will not exceed m_Length-1
    );

    ///This is a stub which should never be called
    virtual bool DecodeDataSubsymbols(unsigned long long StripeID, ///the stripe to be processed
                                      unsigned ErasureSetID, ///identifies the load balancing offset
                                      unsigned SymbolID, ///the symbol to be processed
                                      unsigned SubsymbolID, ///the first subsymbol to be processed
                                      unsigned Subsymbols2Decode,
                                      ///the number of subsymbols within this symbol to be decoded
                                      unsigned char *pDest,
                                      ///destination array. Must have size at least Subsymbols2Decode*m_StripeUnitSize
                                      size_t ThreadID ///the ID of the calling thread
    ) {
        return false;
    };

    ///decode a number of payload subsymbols from a given symbol
    ///@return true on success
    virtual bool DecodeDataSymbols(unsigned long long StripeID, ///the stripe to be processed
                                   unsigned ErasureSetID, ///identifies the load balancing offset
                                   unsigned SymbolID, ///the first symbol to be processed
                                   unsigned Symbols2Decode, ///the number of subsymbols within this symbol to be decoded
                                   unsigned char *pDest,
                                   ///destination array. Must have size at least Subsymbols2Decode*m_StripeUnitSize
                                   size_t ThreadID ///the ID of the calling thread
    );

    ///encode and write the whole stripe
    ///@return true on success
    virtual bool EncodeStripe(unsigned long long StripeID, ///the stripe to be encoded
                              unsigned ErasureSetID, ///identifies the load balancing offset
                              const unsigned char *pData, ///the data to be envoced
                              size_t ThreadID ///the ID of the calling thread
    );

    ///update some information symbols and the corresponding check symbols
    ///@return true on success
    virtual bool UpdateInformationSymbols(unsigned long long StripeID, ///the stripe to be updated,
                                          unsigned ErasureSetID, ///identifies the load balancing offset
                                          unsigned StripeUnitID, ///the first stripe unit to be updated
                                          unsigned Units2Update, ///the number of units to be updated
                                          const unsigned char *pData, ///new payload data symbols
                                          size_t ThreadID ///the ID of the calling thread
    );

    ///check if the codeword is consistent
    virtual bool CheckCodeword(unsigned long long StripeID, ///the stripe to be checked
                               unsigned ErasureSetID, ///identifies the load balancing offset
                               size_t ThreadID ///identifies the calling thread
    );

    ///make a decision about the optimal data encoding method
    ///@return true if one should fetch all the data from disk, update it and write back. Otherwise, only check symbols will be updated
    virtual bool GetEncodingStrategy(unsigned ErasureSetID, ///identifies the load balancing offset
                                     unsigned StripeUnitID, ///the first stripe unit to be updated
                                     unsigned Subsymbols2Encode ///the number of subsymbols to be encoded
    );

public:
    CRSProcessor(RSParams *pParams);

    ~CRSProcessor();
};

///maximal length of Reed-Solomon code
extern const unsigned RSLength;

#endif
