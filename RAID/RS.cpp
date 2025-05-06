/*********************************************************
 * RSRAID.cpp  - implementation for a Reed-Solomon based RAID
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/
#include <algorithm>
#include <string.h>
#include "misc.h"
#include "RS.h"

#ifndef STUDENTBUILD
#include "MultipointEvaluation.h"
#include "FFTProcessors.h"
#endif

using namespace std;
//RS code will be defined over GF(2^8)
#define Extension 8
//length of the parent code
const unsigned RSLength = (1u << Extension) - 1;


//compute X_i^{1-b}/\Lambda'(1/X_i) (needed by Forney algorithm)
int GetForneyMultiple(unsigned LambdaDegree, //the degree of the erasure locator polynomial
                      const GFValue *pLambda, //the polynomial itself
                      int b, //first root
                      int X //the erasure locator (log)
) {
    //X_i^(1-b)
    int Y = (X * (1 - b)) % FieldSize_1;
    //1/X_i^2
    X = (X) ? FieldSize_1 - X : 0;
    X += X;
    if (X >= FieldSize_1) X -= FieldSize_1;

    //evaluate \Lambda' using Horner rule
    //in GF(2^m) one can skip even-numbered components
    GFValue Res = 0;
    for (int j = (LambdaDegree - 1) & ~1; j >= 0; j -= 2) {
        //multiply by 1/X_i^2
        if (Res) {
            int L = LogTable[Res] + X;
            Res = GF[1 + L];
        };
        Res ^= pLambda[j + 1];
    };
    //compute X_i^{1-b}/\Lambda'(1/X_i)
    int R = Y - LogTable[Res];
    if (R < 0)R += FieldSize_1;
    return R;
};


/**Initialize the set of information and check symbols,
setup check symbol locator polynomial*/
CRSProcessor::CRSProcessor(RSParams *pParams): CRAIDProcessor(pParams->CodeDimension + pParams->Redundancy,
                                                              1, pParams, sizeof(RSParams)),
                                               m_Redundancy(pParams->Redundancy),
                                               m_pErasureLocatorsPrime(0),
                                               m_pSyndromes(0), m_pErasureEvaluator(0), m_pErasureLocators(0),
                                               m_pSymbols(0), m_ppSymbols(0), m_pCyclotomicTemp(0) {
    if (m_Dimension >= m_Length)
        throw Exception("Dimension exceeds Reed-Solomon code length");
    if (m_Length > RSLength)
        throw Exception("Reed-Solomon code length exceeds %d", RSLength);
    if (!m_Redundancy)
        throw Exception("Invalid redundancy %d for Reed-Solomon code", m_Redundancy);

#ifndef STUDENTBUILD
    if (pParams->CyclotomicProcessing&&!((m_Redundancy>=MINCYCLOTOMICREDUNDANCY)&&(m_Redundancy<=MAXCYCLOTOMICREDUNDANCY)))
        throw Exception("The cyclotomic processor does not support redundancy %d",m_Redundancy);
    if (pParams->OptimizedCheckLocators&&(m_Redundancy>MAXOPTIMIZEDREDUNDANCY))
        throw Exception("Optimized check locators are not available for redundancy %d",m_Redundancy);
    m_CyclotomicProcessing=pParams->CyclotomicProcessing;
    m_OptimizedCheckLocators=pParams->OptimizedCheckLocators;
#endif


    InitGF(Extension);

    m_pInfSymbols = new int[m_Dimension];
    m_pCheckSymbols = new int[m_Redundancy];
#ifndef STUDENTBUILD
    if (m_OptimizedCheckLocators)
    {
        //optimized locator selection
        memcpy(m_pCheckSymbols,CheckLocators[m_Redundancy-1].pLocators,sizeof(int)*m_Redundancy);
    }else
#endif
    {
        for (unsigned i = 0; i < m_Redundancy; i++)
            m_pCheckSymbols[i] = RSLength - m_Redundancy + i;
    };

    //the information symbols can be any except those reserved for check symbols
    unsigned j = 0;
#ifndef STUDENTBUILD
     if (m_CyclotomicProcessing)
     {
            const int* pOrdering=CyclotomicOrderings[m_Redundancy-1];
            for(unsigned i=0;i<RSLength;i++)
            {
                //make sure that this symbol is not used as a check one
                bool Valid=true;
                for(unsigned s=0;s<m_Redundancy;s++)
                    if (m_pCheckSymbols[s]==pOrdering[i])
                    {
                        Valid=false;
                        break;
                    };
                if (Valid)
                {
                    m_pInfSymbols[j++]=pOrdering[i];
                    if (j==m_Dimension)
                        break;
                };
            };

     }else
#endif
    {
        for (unsigned i = 0; i < RSLength; i++) {
            //make sure that this symbol is not used as a check one
            bool Valid = true;
            for (unsigned s = 0; s < m_Redundancy; s++)
                if (m_pCheckSymbols[s] == i) {
                    Valid = false;
                    break;
                };
            if (Valid) {
                m_pInfSymbols[j++] = i;
                if (j == m_Dimension)
                    break;
            };
        };
    };

    //construct the erasure locator polynomial for efficient encoding
    m_pCheckLocator = new GFValue[m_Redundancy + 1];
    memset(m_pCheckLocator + 1, 0, sizeof(GFValue) * m_Redundancy);
    m_pCheckLocator[0] = 1;
    ///\prod_i (1-xX_i)
    for (unsigned i = 0; i < m_Redundancy; i++) {
        for (unsigned j = i + 1; j > 0; j--) {
            if (m_pCheckLocator[j - 1]) {
                int L = LogTable[m_pCheckLocator[j - 1]] + m_pCheckSymbols[i];
                if (L >= FieldSize_1)
                    L -= FieldSize_1;
                m_pCheckLocator[j] ^= GF[1 + L];
            };
        };
    };
    m_pCheckLocatorsPrime = new int[m_Redundancy];
    for (unsigned i = 0; i < m_Redundancy; i++)
        m_pCheckLocatorsPrime[i] = GetForneyMultiple(m_Redundancy, m_pCheckLocator, 0, m_pCheckSymbols[i]);

    m_pErasureLocators = new GFValue[(m_Redundancy + 1) * m_Length * m_InterleavingOrder];
    m_pErasureLocatorsPrime = new int[m_Redundancy * m_Length * m_InterleavingOrder];
};


CRSProcessor::~CRSProcessor() {
    delete[]m_pInfSymbols;
    delete[]m_pCheckSymbols;
    delete[]m_pCheckLocator;
    delete[]m_pErasureLocators;
    delete[]m_ppSymbols;
    delete[]m_pErasureLocatorsPrime;
    delete[]m_pCheckLocatorsPrime;
    AlignedFree(m_pErasureEvaluator);
    AlignedFree(m_pSyndromes);
    AlignedFree(m_pSymbols);
    AlignedFree(m_pCyclotomicTemp);
};


/**
Allocate memory for syndrome and erasure evaluator polynomials
return true on success
*/
bool CRSProcessor::Attach(CDiskArray *pArray, ///the disk array
                          unsigned ConcurrentThreads
                          ///the number of concurrent processing threads that will make calls to the processor
) {
    m_pSyndromes = AlignedMalloc(m_Redundancy * m_StripeUnitSize * ConcurrentThreads);
    m_pErasureEvaluator = AlignedMalloc(m_Redundancy * m_StripeUnitSize * ConcurrentThreads);
    m_pSymbols = AlignedMalloc(m_Length * m_StripeUnitSize * ConcurrentThreads);
    if (m_CyclotomicProcessing) {
#ifndef STUDENTBUILD
        m_pCyclotomicTemp=AlignedMalloc(sizeof(GFValue)*m_StripeUnitSize*ConcurrentThreads*CYCLOTOMIC_TEMP_SIZE);
#endif
    };


    m_ppSymbols = new const GFValue *[RSLength * ConcurrentThreads];
    memset(m_ppSymbols, 0, RSLength * ConcurrentThreads * sizeof(GFValue *));
    return CRAIDProcessor::Attach(pArray, ConcurrentThreads);
};
///reset the erasure correction engine
/// this will be called if the set of failed disks changes
void CRSProcessor::ResetErasures() {
    CRAIDProcessor::ResetErasures();
};


/**Construct the erasure locator polynomial for a given combination
of erasures
@return true if the specified combination of erasures is correctable
*/
bool CRSProcessor::IsCorrectable(
    unsigned ErasureSetID ///identifies the erasure combination. This will not exceed m_Length-1
) {
    if (GetNumOfErasures(ErasureSetID) == 0) return true;
    if (GetNumOfErasures(ErasureSetID) > m_Redundancy) return false;
    //setup the erasure locator polynomial for decoding
    GFValue *pLambda = m_pErasureLocators + ErasureSetID * (m_Redundancy + 1);
    pLambda[0] = 1;
    unsigned t = GetNumOfErasures(ErasureSetID);
    memset(pLambda + 1, 0, sizeof(GFValue) * t);
    //Lambda(x)=\prod_{i=0}^{t-1}(1-xX_i)
    for (unsigned i = 0; i < t; i++) {
        unsigned DiskID = GetErasedPosition(ErasureSetID, i);
        //translate it into the locator value
        int Locator;
        if (DiskID < m_Dimension)
            Locator = m_pInfSymbols[DiskID];
        else
            Locator = m_pCheckSymbols[DiskID - m_Dimension];
        //multiply by (1-xX_i)
        for (unsigned j = i + 1; j > 0; j--) {
            if (pLambda[j - 1]) {
                int L = LogTable[pLambda[j - 1]] + Locator;
                if (L >= FieldSize_1)
                    L -= FieldSize_1;
                pLambda[j] ^= GF[1 + L];
            };
        };
    };
    //compute \alpha^{1-b}/\Lambda'(1/X_i) (needed by Forney algorithm)
    int *pLambdaPrime = m_pErasureLocatorsPrime + ErasureSetID * m_Redundancy;
    for (unsigned i = 0; i < t; i++) {
        unsigned DiskID = GetErasedPosition(ErasureSetID, i);
        int Locator;
        if (DiskID < m_Dimension)
            Locator = m_pInfSymbols[DiskID];
        else
            Locator = m_pCheckSymbols[DiskID - m_Dimension];
        pLambdaPrime[i] = GetForneyMultiple(t, pLambda, 0, Locator);
    };

    return true;
};

/** Compute
S_i=\sum_{j=0}^{n-1} y_j \alpha^{ij}, Low<=i<High.
The operation will be performed for each of UnitSize words
*/
void ComputeSyndrome(const GFValue *const*ppData,
                     ///pointers to y_j. Each y_j is an array of size UnitSize. If ppData[i]=0, the symbols are assumed to be 0
                     GFValue *pSyndromes, /// S_i (packed). Each S_i is an array of size UnitSize
                     unsigned Low, ///low range of syndrome indices
                     unsigned High, ///upper limit on syndrome indices
                     unsigned UnitSize ///size of one data unit
) {
    if (ppData[0]) {
        for (unsigned i = 0; i < High - Low; i++)
            memcpy(pSyndromes + i * UnitSize, ppData[0], UnitSize * sizeof(GFValue));
    } else memset(pSyndromes, 0, UnitSize * (High - Low) * sizeof(GFValue));

    for (int i = 1; i < RSLength; i++) {
        if (ppData[i]) {
            int L1 = (Low * i) % FieldSize_1;
            for (unsigned j = 0; j < High - Low; j++, L1 += i) {
                if (L1 >= FieldSize_1)
                    L1 -= FieldSize_1;
                MultiplyAdd(L1, ppData[i], pSyndromes + j * UnitSize, UnitSize);
            };
        };
    };
};


/**Compute the erasure evaluator polynomial as
\Gamma(x)=\Lambda(x)*S(x)\bmod x^t

This is performed simultaneously for all symbols in the data block*/
void GetErasureEvaluator(const GFValue *pSyndromes, ///S(x) stored blockwise
                         const GFValue *pLambda, ///\Lambda(x)
                         GFValue *pGamma, ///the erasure evaluator polynomial to be computed
                         unsigned MaxErrors, ///the degree of \Lambda(x)
                         unsigned UnitSize ///size of one data block
) {
    for (unsigned i = 0; i < MaxErrors; i++) {
        GFValue *pDest = pGamma + i * UnitSize;
        //\Lambda(0)=1, so we can just copy
        memcpy(pDest, pSyndromes + i * UnitSize, sizeof(GFValue) * UnitSize);
        //\Gamma_i=\sum_{j=0}^i S_{i-j}\Lambda_j
        for (unsigned j = 1; j <= i; j++) {
            if (pLambda[j])
                MultiplyAdd(LogTable[pLambda[j]], pSyndromes + (i - j) * UnitSize, pDest, UnitSize);
        };
    };
};

/** evaluate a block of polynomials at a given point
The degree must be >=0
*/
void Evaluate(const GFValue *pPolynomials, ///polynomials stored blockwise
              unsigned Degree, ///the degree of the polynomials
              int x, ///the polynomials will be evaluated here
              GFValue *pValue, ///output will be placed here
              unsigned UnitSize ///the number of polynomials to be evaluated (i.e. block size)
) {
    memcpy(pValue, pPolynomials, sizeof(GFValue) * UnitSize);
    int L = x;
    for (unsigned i = 1; i <= Degree; i++, L += x) {
        if (L >= FieldSize_1) L -= FieldSize_1;
        MultiplyAdd(L, pPolynomials + i * UnitSize, pValue, UnitSize);
    };
}


/**
Fetch the non-erased symbols as is. If there are erased symbols,
compute the syndrome S(x)=\sum_{i=0}^{t-1} S_i,
where t is the number of erasures, S_i=\sum_{j=0}^{n-1} y_j X_j^i, X_j is the locator of the j-th symbol, y_j is assumed
to be zero if it is erased,
construct the erasure evaluator polynomial \Gamma(x)=S(x)\Lambda(x)\bmod x^{t},
and recover the erasures via Forney algorithm
*/
bool CRSProcessor::DecodeDataSymbols(unsigned long long StripeID, ///the stripe to be processed
                                     unsigned ErasureSetID, ///identifies the load balancing offset
                                     unsigned SymbolID, ///the first symbol to be processed
                                     unsigned Symbols2Decode,
                                     ///the number of subsymbols within this symbol to be decoded
                                     unsigned char *pDest,
                                     ///destination array. Must have size at least Subsymbols2Decode*m_StripeUnitSize
                                     size_t ThreadID ///the ID of the calling thread
) {
    bool NeedsDecoding = false;
    //pointers to fetched data
    const GFValue **ppData = m_ppSymbols + RSLength * ThreadID;
    for (unsigned i = 0; i < Symbols2Decode; i++) {
        unsigned S = SymbolID + i;
        if (IsErased(ErasureSetID, S)) {
            NeedsDecoding = true;
            ppData[m_pInfSymbols[S]] = 0;
        } else {
            //fetch it
            if (!ReadStripeUnit(StripeID, ErasureSetID, SymbolID + i, 0, 1, pDest + i * m_StripeUnitSize))
                return false;
            //save the pointer if we need it for decoding
            ppData[m_pInfSymbols[S]] = pDest + i * m_StripeUnitSize;
        };
    };
    if (NeedsDecoding) {
        GFValue *pFetchBuffer = m_pSymbols + ThreadID * m_Length * m_StripeUnitSize;
        //fetch all surviving information symbols
        for (unsigned i = 0; i < SymbolID; i++) {
            if (IsErased(ErasureSetID, i)) {
                ppData[m_pInfSymbols[i]] = 0;
            } else {
                //fetch it
                ppData[m_pInfSymbols[i]] = pFetchBuffer + i * m_StripeUnitSize;
                if (!ReadStripeUnit(StripeID, ErasureSetID, i, 0, 1, pFetchBuffer + i * m_StripeUnitSize))
                    return false;
            };
        };
        for (unsigned i = SymbolID + Symbols2Decode; i < m_Dimension; i++) {
            if (IsErased(ErasureSetID, i)) {
                ppData[m_pInfSymbols[i]] = 0;
            } else {
                //fetch it
                ppData[m_pInfSymbols[i]] = pFetchBuffer + i * m_StripeUnitSize;
                if (!ReadStripeUnit(StripeID, ErasureSetID, i, 0, 1, pFetchBuffer + i * m_StripeUnitSize))
                    return false;
            };
        };
        //fetch all surviving check symbols
        for (unsigned i = 0; i < m_Redundancy; i++) {
            if (IsErased(ErasureSetID, m_Dimension + i)) {
                ppData[m_pCheckSymbols[i]] = 0;
            } else {
                //fetch it
                ppData[m_pCheckSymbols[i]] = pFetchBuffer + (m_Dimension + i) * m_StripeUnitSize;
                if (!ReadStripeUnit(StripeID, ErasureSetID, m_Dimension + i, 0, 1,
                                    pFetchBuffer + (m_Dimension + i) * m_StripeUnitSize))
                    return false;
            };
        };
        GFValue *pSyndrome = m_pSyndromes + ThreadID * m_Redundancy * m_StripeUnitSize;
        GFValue *pErasureEvaluator = m_pErasureEvaluator + ThreadID * m_Redundancy * m_StripeUnitSize;
#ifndef STUDENTBUILD
		if (m_CyclotomicProcessing)
		{
		//   ComputeSyndrome(ppData,pSyndrome,m_FirstRoot,m_FirstRoot+m_Redundancy,m_StripeUnitSize);
			ComputeSyndromeCyclotomic(ppData,pSyndrome,m_Redundancy,m_pCyclotomicTemp+ThreadID*CYCLOTOMIC_TEMP_SIZE*m_StripeUnitSize,pErasureEvaluator,m_StripeUnitSize);
		}else
#endif
        ComputeSyndrome(ppData, pSyndrome, 0, m_Redundancy, m_StripeUnitSize);

        GetErasureEvaluator(pSyndrome, m_pErasureLocators + ErasureSetID * (m_Redundancy + 1), pErasureEvaluator,
                            GetNumOfErasures(ErasureSetID), m_StripeUnitSize);
        //recover the erasures
        const int *pErasureLocatorsPrime = m_pErasureLocatorsPrime + ErasureSetID * m_Redundancy;
        for (unsigned i = 0; i < GetNumOfErasures(ErasureSetID); i++) {
            unsigned S = GetErasedPosition(ErasureSetID, i);
            if (S < SymbolID) continue;
            if (S >= SymbolID + Symbols2Decode)continue;
            int X = (m_pInfSymbols[S]) ? FieldSize_1 - m_pInfSymbols[S] : 0;
            GFValue *pCurDest = pDest + (S - SymbolID) * m_StripeUnitSize;
            //\Gamma(1/X_i)
            Evaluate(pErasureEvaluator, GetNumOfErasures(ErasureSetID) - 1, X, pCurDest, m_StripeUnitSize);
            //\alpha^{1-b}\Gamma(1/X_i)/\Lambda'(1/X_i)
            Multiply(pErasureLocatorsPrime[i], pCurDest, pCurDest, m_StripeUnitSize);
        };
    };
    return true;
};
///encode and write the whole stripe
///@return true on success
bool CRSProcessor::EncodeStripe(unsigned long long StripeID, ///the stripe to be encoded
                                unsigned ErasureSetID, ///identifies the load balancing offset
                                const unsigned char *pData, ///the data to be envoced
                                size_t ThreadID ///the ID of the calling thread
) {
    //initialize the information symbol positions and compute check ones via erasure decoding
    const GFValue **ppData = m_ppSymbols + RSLength * ThreadID;
    for (unsigned i = 0; i < m_Dimension; i++) {
        ppData[m_pInfSymbols[i]] = pData + i * m_StripeUnitSize;
        //send the data to disk
        WriteStripeUnit(StripeID, ErasureSetID, i, 0, 1, pData + i * m_StripeUnitSize);
    };
    for (unsigned i = 0; i < m_Redundancy; i++)
        ppData[m_pCheckSymbols[i]] = 0;

    GFValue *pSyndrome = m_pSyndromes + ThreadID * m_Redundancy * m_StripeUnitSize;
    GFValue *pErasureEvaluator = m_pErasureEvaluator + ThreadID * m_Redundancy * m_StripeUnitSize;
#ifndef STUDENTBUILD
    if (m_CyclotomicProcessing)
    {
//        ComputeSyndrome(ppData,pSyndrome,m_FirstRoot,m_FirstRoot+m_Redundancy,m_StripeUnitSize);
		ComputeSyndromeCyclotomic(ppData,pSyndrome,m_Redundancy,m_pCyclotomicTemp+ThreadID*CYCLOTOMIC_TEMP_SIZE*m_StripeUnitSize,pErasureEvaluator,m_StripeUnitSize);
    }else
#endif
    ComputeSyndrome(ppData, pSyndrome, 0, m_Redundancy, m_StripeUnitSize);


    GetErasureEvaluator(pSyndrome, m_pCheckLocator, pErasureEvaluator, m_Redundancy, m_StripeUnitSize);
#ifndef STUDENTBUILD
   if (m_OptimizedCheckLocators)
    {
        //use pSyndrome as a temporary storage
        //\Gamma(1/X_i)
        CheckLocators[m_Redundancy-1].Evaluator(pErasureEvaluator,pSyndrome,m_StripeUnitSize);
        //recover the erased check symbols
        for(unsigned i=0;i<m_Redundancy;i++)
        {
            int X=(m_pCheckSymbols[i])?FieldSize_1-m_pCheckSymbols[i]:0;
            //Evaluate(pErasureEvaluator,m_Redundancy-1,X,pSyndrome+i*m_StripeUnitSize,m_StripeUnitSize);

            //X_i^{1-b}\Gamma(1/X_i)/\Lambda'(1/X_i)
            Multiply(m_pCheckLocatorsPrime[i],pSyndrome+i*m_StripeUnitSize,pSyndrome+i*m_StripeUnitSize,m_StripeUnitSize);
            //send check symbols to disk
            WriteStripeUnit(StripeID,ErasureSetID,m_Dimension+i,0,1,pSyndrome+i*m_StripeUnitSize);
        };

    }else
#endif
    {
        //recover the erased check symbols
        for (unsigned i = 0; i < m_Redundancy; i++) {
            int X = (m_pCheckSymbols[i]) ? FieldSize_1 - m_pCheckSymbols[i] : 0;
            //use pSyndrome as a temporary storage
            //\Gamma(1/X_i)
            Evaluate(pErasureEvaluator, m_Redundancy - 1, X, pSyndrome, m_StripeUnitSize);
            //X_i^{1-b}\Gamma(1/X_i)/\Lambda'(1/X_i)
            Multiply(m_pCheckLocatorsPrime[i], pSyndrome, pSyndrome, m_StripeUnitSize);
            //send check symbols to disk
            WriteStripeUnit(StripeID, ErasureSetID, m_Dimension + i, 0, 1, pSyndrome);
        };
    };

    return true;
};

/** If the range of symbols to be updated does not include any erasures, use the standard decision rule.
Otherwise, do full fetch-combine-encode cycle
*/
bool CRSProcessor::GetEncodingStrategy(unsigned ErasureSetID, ///identifies the load balancing offset
                                       unsigned StripeUnitID, ///the first stripe unit to be updated
                                       unsigned Subsymbols2Encode ///the number of subsymbols to be encoded
) {
    for (unsigned i = 0; i < GetNumOfErasures(ErasureSetID); i++) {
        if (GetErasedPosition(ErasureSetID, i) < (int) StripeUnitID)
            continue; //it does not belong into the update set
        if (GetErasedPosition(ErasureSetID, i) >= int(StripeUnitID + Subsymbols2Encode))
            continue; //it does not belong into the update set

        //erased symbol has to be updated, do full encoding
        return true;
    };
    return CRAIDProcessor::GetEncodingStrategy(ErasureSetID, StripeUnitID, Subsymbols2Encode);
};
/**update some information symbols and the corresponding check symbols
This will fetch old values of the symbols to be updated, compute the corresponding syndrome and
   @return true on success

   */
bool CRSProcessor::UpdateInformationSymbols(unsigned long long StripeID, ///the stripe to be updated,
                                            unsigned ErasureSetID, ///identifies the load balancing offset
                                            unsigned StripeUnitID, ///the first stripe unit to be updated
                                            unsigned Units2Update, ///the number of units to be updated
                                            const unsigned char *pData, ///new payload data symbols
                                            size_t ThreadID ///the ID of the calling thread
) {
    //assume here that there are no erasures
    GFValue *pFetchBuffer = m_pSymbols + ThreadID * m_Length * m_StripeUnitSize;
    const GFValue **ppData = m_ppSymbols + RSLength * ThreadID;
    memset(ppData, 0, RSLength * sizeof(ppData[0]));
    bool Result = true;
    for (unsigned i = 0; i < Units2Update; i++) {
        //find the difference between new and old values
        GFValue *pCurSymbol = pFetchBuffer + i * m_StripeUnitSize;
        Result &= ReadStripeUnit(StripeID, ErasureSetID, StripeUnitID + i, 0, 1, pCurSymbol);
        XOR(pCurSymbol, pData + i * m_StripeUnitSize, m_StripeUnitSize);
        ppData[m_pInfSymbols[StripeUnitID + i]] = pCurSymbol;
        //save the new value
        Result &= WriteStripeUnit(StripeID, ErasureSetID, StripeUnitID + i, 0, 1, pData + i * m_StripeUnitSize);
    };
    GFValue *pSyndrome = m_pSyndromes + ThreadID * m_Redundancy * m_StripeUnitSize;
    GFValue *pErasureEvaluator = m_pErasureEvaluator + ThreadID * m_Redundancy * m_StripeUnitSize;
#ifndef STUDENTBUILD
    if (m_CyclotomicProcessing)
    {
     //   ComputeSyndrome(ppData,pSyndrome,m_FirstRoot,m_FirstRoot+m_Redundancy,m_StripeUnitSize);
		ComputeSyndromeCyclotomic(ppData,pSyndrome,m_Redundancy,m_pCyclotomicTemp+ThreadID*CYCLOTOMIC_TEMP_SIZE*m_StripeUnitSize,pErasureEvaluator,m_StripeUnitSize);

    }
    elseUpdate
#endif
    ComputeSyndrome(ppData, pSyndrome, 0, m_Redundancy, m_StripeUnitSize);

    GetErasureEvaluator(pSyndrome, m_pCheckLocator, pErasureEvaluator, m_Redundancy, m_StripeUnitSize);
#ifndef STUDENTBUILD
    if (m_OptimizedCheckLocators)
    {
        //use pSyndrome as a temporary storage
        //\Gamma(1/X_i)
        CheckLocators[m_Redundancy-1].Evaluator(pErasureEvaluator,pSyndrome,m_StripeUnitSize);
        //recover the erased check symbols
        for(unsigned i=0;i<m_Redundancy;i++)
        {
            int X=(m_pCheckSymbols[i])?FieldSize_1-m_pCheckSymbols[i]:0;
            //fetch old value
            Result&=ReadStripeUnit(StripeID,ErasureSetID,m_Dimension+i,0,1,pFetchBuffer);
            //X_i^{1-b}\Gamma(1/X_i)/\Lambda'(1/X_i)
            MultiplyAdd(m_pCheckLocatorsPrime[i],pSyndrome+i*m_StripeUnitSize,pFetchBuffer,m_StripeUnitSize);
            //send check symbols to disk
            WriteStripeUnit(StripeID,ErasureSetID,m_Dimension+i,0,1,pFetchBuffer);
        };

    }else
#endif
    {
        //recover the erased check symbols
        for (unsigned i = 0; i < m_Redundancy; i++) {
            if (IsErased(ErasureSetID, m_Dimension + i))
                //no need to update this symbol
                continue;
            int X = (m_pCheckSymbols[i]) ? FieldSize_1 - m_pCheckSymbols[i] : 0;
            //use pSyndrome as a temporary storage
            //\Gamma(1/X_i)
            Evaluate(pErasureEvaluator, m_Redundancy - 1, X, pSyndrome, m_StripeUnitSize);
            //fetch old value
            Result &= ReadStripeUnit(StripeID, ErasureSetID, m_Dimension + i, 0, 1, pFetchBuffer);
            //add to it X_i^{1-b}\Gamma(1/X_i)/\Lambda'(1/X_i)
            MultiplyAdd(m_pCheckLocatorsPrime[i], pSyndrome, pFetchBuffer, m_StripeUnitSize);
            //write it back
            Result &= WriteStripeUnit(StripeID, ErasureSetID, m_Dimension + i, 0, 1, pFetchBuffer);
        };
    };

    return true;
};
/**
   Fetch all codeword symbols, compute the syndrome and check if it is zero
*/
bool CRSProcessor::CheckCodeword(unsigned long long StripeID, ///the stripe to be checked
                                 unsigned ErasureSetID, ///identifies the load balancing offset
                                 size_t ThreadID ///identifies the calling thread
) {
    if (GetNumOfErasures(ErasureSetID))
        return true;
    GFValue *pFetchBuffer = m_pSymbols + ThreadID * m_Length * m_StripeUnitSize;
    const GFValue **ppData = m_ppSymbols + RSLength * ThreadID;
    //fetch information symbols
    for (unsigned i = 0; i < m_Dimension; i++) {
        GFValue *pCurSymbol = pFetchBuffer + i * m_StripeUnitSize;
        if (!ReadStripeUnit(StripeID, ErasureSetID, i, 0, 1, pCurSymbol)) return false;
        ppData[m_pInfSymbols[i]] = pCurSymbol;
    };
    //fetch check symbols
    for (unsigned i = 0; i < m_Redundancy; i++) {
        GFValue *pCurSymbol = pFetchBuffer + (m_Dimension + i) * m_StripeUnitSize;
        if (!ReadStripeUnit(StripeID, ErasureSetID, m_Dimension + i, 0, 1, pCurSymbol)) return false;
        ppData[m_pCheckSymbols[i]] = pCurSymbol;
    };
    GFValue *pSyndrome = m_pSyndromes + ThreadID * m_Redundancy * m_StripeUnitSize;
    ComputeSyndrome(ppData, pSyndrome, 0, m_Redundancy, m_StripeUnitSize);
    GFValue X = 0;
    for (unsigned i = 0; i < m_Redundancy * m_StripeUnitSize; i++)
        X |= pSyndrome[i];

    return X == 0;
}
