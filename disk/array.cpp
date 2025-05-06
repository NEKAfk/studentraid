/*********************************************************
 * array.cpp  - implementation of the  RAID emulator
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/

#include <iostream>
#include <time.h>
#include <string.h>
#include "misc.h"
#include "array.h"
#include "arithmetic.h"

using namespace std;

///initialize the array. The array parameters
///will be extracted from the processor object

CDiskArray::CDiskArray(unsigned NumberOfDisks,
                       ///the number of disks to be emulated. This must be not less
                                             ///than the length of the array code implemented by Processor
                                             ///All extra disks will be ignored
                       DiskConf const *pDiskFiles, ///configuration of the emulated disks
                       size_t DiskCapacity, ///the capacity of a single disk
                       CRAIDProcessor &Processor, ///provides encoding and decoding functionality
                       unsigned NumOfThreads ///the number of concurrent processing threads
) : m_NumOfThreads(NumOfThreads), m_Engine(Processor),
    m_MountState(msUnmounted), m_NumOfDisks(NumberOfDisks),
    m_StripeUnitSize(Processor.GetStripeUnitSize()),
    m_UnitsPerStripePrim(Processor.GetStripeUnitsPerSymbol() * Processor.GetDimension()),
    m_UnitsPerStripe(m_UnitsPerStripePrim * Processor.GetInterleavingOrder()),
    m_NumOfStripes(DiskCapacity / (Processor.GetStripeUnitSize() *
                                   Processor.GetStripeUnitsPerSymbol())),
    m_StripeSize(m_UnitsPerStripe * m_StripeUnitSize), m_Locker(NumOfThreads) {
    if (Processor.GetCodeLength() * Processor.GetInterleavingOrder() > m_NumOfDisks)
        throw Exception("Not enough disks for a given code (minimum %d is required)",
                        Processor.GetCodeLength() * Processor.GetInterleavingOrder());
    else m_NumOfDisks = Processor.GetCodeLength() * Processor.GetInterleavingOrder();

    void const *pCodeConfig;
    unsigned CodeConfigSize = Processor.GetConfiguration(pCodeConfig);
    //attach the disks
    m_pDisks = new CDisk[m_NumOfDisks];
    time_t LastArrayMount = 0;
    for (unsigned i = 0; i < m_NumOfDisks; i++) {
        if (m_pDisks[i].Initialize(pDiskFiles[i].pFileName, i, m_StripeUnitSize,
                                   m_NumOfStripes * Processor.GetStripeUnitsPerSymbol(),
                                   CodeConfigSize)) {
            //check if the array configuration stored on disk is the same as the one of the processor
            void const *pCodeConfig2;
            unsigned CodeConfigSize2 = m_pDisks[i].GetArrayData(pCodeConfig2);
            if ((CodeConfigSize2 != CodeConfigSize) || memcmp(pCodeConfig, pCodeConfig2, CodeConfigSize)) {
                //cerr << "Array configuration mismatch for disk " << i << endl;
                m_pDisks[i].SetDiskState(dsInvalid);
            };


            if (m_pDisks[i].GetDiskState() == dsOffline) {
                //identify the latest mounted disks
                time_t LastMount = m_pDisks[i].GetLastUnmountTime();
                if (LastMount > LastArrayMount)
                    LastArrayMount = LastMount;
            };
        };
    };
    unsigned NumOfInitializedDisks = 0;
    unsigned NumOfOnlineDisks = 0;
    //take online the latest mounted disks
    for (unsigned i = 0; i < m_NumOfDisks; i++) {
        if ((m_pDisks[i].GetDiskState() == dsOffline) && pDiskFiles[i].Online) {
            NumOfInitializedDisks++;
            if (m_pDisks[i].GetLastUnmountTime() == LastArrayMount) {
                m_pDisks[i].SetDiskState(dsOnline);
                NumOfOnlineDisks++;
            } else
            //there were data modifications since the last disk mount
                m_pDisks[i].SetDiskState(dsInvalid);
        };
    };
    //make final initialization of the coding engine
    m_Engine.Attach(this, NumOfThreads);
    if (NumOfInitializedDisks == 0)
        m_ArrayState = asUninitialized;
    else {
        if (NumOfOnlineDisks == m_NumOfDisks)
            m_ArrayState = asNormal;
        else {
            ///check if the array is mountable
            if (Processor.IsMountable())
                m_ArrayState = asDegraded;
            else
                m_ArrayState = asFailed;
        };
    };
    m_pPartialRWBuffer = AlignedMalloc(m_StripeUnitSize * m_NumOfThreads);
};

CDiskArray::~CDiskArray() {
    Unmount();
    delete[]m_pDisks;
    AlignedFree(m_pPartialRWBuffer);
};

///enable data access
///@return true on success, false if not mountable
bool CDiskArray::Mount(bool Write ///true if write access is needed
) {
    if ((m_ArrayState == asFailed) || (m_ArrayState == asUninitialized))
        return false;
    if (m_MountState != msUnmounted)
        return false;
    bool Result = true;
    //mount the underlying disks
    for (unsigned i = 0; i < m_NumOfDisks; i++)
        if (m_pDisks[i].GetDiskState() == dsOnline)
            Result &= m_pDisks[i].Mount(Write);

    if (Result)
        m_MountState = (Write) ? msReadWrite : msRead;
    else
        //this should not happen
        throw Exception("Unexpected mount failure");
    return Result;
};


///disable data access
///@return true on success
bool CDiskArray::Unmount() {
    if (m_MountState == msUnmounted)
        return false;
    m_MountState = msUnmounted;
    //unmount all the disks and put the timestamp if necessary
    bool Result = true;
    time_t Timestamp = time(NULL);
    for (unsigned i = 0; i < m_NumOfDisks; i++)
        Result &= m_pDisks[i].Unmount(Timestamp);

    return Result;
};

///initialize the array. It must be unmounted
///@return true on success
bool CDiskArray::Init() {
    if (m_MountState != msUnmounted)
        //illegal array access
        return false;
    m_ArrayState = asUninitialized;
    bool Result = true;
    const void *pArrayData;
    unsigned DataSize = m_Engine.GetConfiguration(pArrayData);
    for (unsigned i = 0; i < m_NumOfDisks; i++) {
        if (m_pDisks[i].GetDiskState() == dsOnline)
            m_pDisks[i].SetDiskState(dsOffline);
        m_pDisks[i].SetArrayData(pArrayData, DataSize);
        Result &= m_pDisks[i].ResetDisk();
    };
    if (Result) {
        //reset the erasure configuration
        m_Engine.ResetErasures();
        //check if the array is mountable
        if (m_Engine.IsMountable())
            m_ArrayState = asNormal;
        else {
            m_ArrayState = asFailed;
            throw Exception("Cannot mount the initialized array");
        };
    };
    return Result;
};


///read a number of stripe units. The array must be mounted
///@return true on success
bool CDiskArray::Read(unsigned long long StripeUnitID, ///the first stripe unit
                      unsigned long long Units2Read, ///the number of stripe units to be read
                      unsigned char *pDest, ///destination buffer. Must have size for Units2Read*m_StripeUnitSize bytes
                      size_t ThreadID ///the ID of a calling thread obtained from m_Locker
) {
    if (m_MountState == msUnmounted)
        return false;
    unsigned long long StripeID = StripeUnitID / m_UnitsPerStripe;
    unsigned UnitID = StripeUnitID % m_UnitsPerStripe;
    bool Result = true;
    /*    if (UnitID)
        {
            //partial read
            unsigned CurUnits2Read=(unsigned)min((unsigned long long)(m_UnitsPerStripe-UnitID),Units2Read);
            Result&=m_Engine.ReadData(StripeID,UnitID,CurUnits2Read,pDest,ThreadID);
            pDest+=CurUnits2Read*m_StripeUnitSize;
            StripeID++;
            Units2Read-=CurUnits2Read;
            UnitID=0;
        };
        while(Result&&(Units2Read>=m_UnitsPerStripe))
        {
            Result&=m_Engine.ReadData(StripeID++,0,m_UnitsPerStripe,pDest,ThreadID);
            pDest+=m_StripeSize;
            Units2Read-=m_UnitsPerStripe;
        };
        if(Result&&(Units2Read>0))
        {
            Result&=m_Engine.ReadData(StripeID++,0,(unsigned)Units2Read,pDest,ThreadID);
        };*/
    unsigned InterleavedID = UnitID / m_UnitsPerStripePrim;
    unsigned CurUnit = UnitID % m_UnitsPerStripePrim;
    while (Result && Units2Read) {
        unsigned CurUnits2Read = (unsigned) min((unsigned long long) (m_UnitsPerStripePrim - CurUnit), Units2Read);
        Result &= m_Engine.ReadData(StripeID, CurUnit, InterleavedID, CurUnits2Read, pDest, ThreadID);
        pDest += CurUnits2Read * m_StripeUnitSize;
        Units2Read -= CurUnits2Read;
        CurUnit = 0;
        InterleavedID++;
        if (InterleavedID == m_Engine.GetInterleavingOrder()) {
            InterleavedID = 0;
            StripeID++;
        };
    };
    return Result;
};

///write a number of stripe units. The array must be write-mounted
///@return true on success
bool CDiskArray::Write(unsigned long long StripeUnitID, ///the first stripe unit
                       unsigned long long Units2Write, ///the number of stripe units to be written
                       const unsigned char *pSrc,
                       ///source buffer. Must have size for Units2Write*m_StripeUnitSize bytes
                       size_t ThreadID ///the ID of a calling thread obtained from m_Locker
) {
    if (m_MountState != msReadWrite)
        return false;
    unsigned long long StripeID = StripeUnitID / m_UnitsPerStripe;
    unsigned UnitID = StripeUnitID % m_UnitsPerStripe;
    bool Result = true;
    /* if (UnitID)
     {
         //partial write
         unsigned CurUnits2Write=(unsigned)min((unsigned long long)(m_UnitsPerStripe-UnitID),Units2Write);
         Result&=m_Engine.WriteData(StripeID,UnitID,CurUnits2Write,pSrc,ThreadID);
         pSrc+=CurUnits2Write*m_StripeUnitSize;
         StripeID++;
         Units2Write-=CurUnits2Write;
         UnitID=0;
     };
     while(Result&&(Units2Write>=m_UnitsPerStripe))
     {
         Result&=m_Engine.WriteData(StripeID++,0,m_UnitsPerStripe,pSrc,ThreadID);
         pSrc+=m_StripeSize;
         Units2Write-=m_UnitsPerStripe;
     };
     if(Result&&(Units2Write>0))
     {
         Result&=m_Engine.WriteData(StripeID++,0,(unsigned)Units2Write,pSrc,ThreadID);
     };*/
    unsigned InterleavedID = UnitID / m_UnitsPerStripePrim;
    unsigned CurUnit = UnitID % m_UnitsPerStripePrim;
    while (Result && Units2Write) {
        unsigned CurUnits2Write = (unsigned) min((unsigned long long) (m_UnitsPerStripePrim - CurUnit), Units2Write);
        Result &= m_Engine.WriteData(StripeID, CurUnit, InterleavedID, CurUnits2Write, pSrc, ThreadID);
        pSrc += CurUnits2Write * m_StripeUnitSize;
        Units2Write -= CurUnits2Write;
        CurUnit = 0;
        InterleavedID++;
        if (InterleavedID == m_Engine.GetInterleavingOrder()) {
            InterleavedID = 0;
            StripeID++;
        };
    };

    return Result;
};

///check if the array is consistend
///@return true on success
bool CDiskArray::Check() {
    eMountState OldState = m_MountState;
    size_t LockID = m_Locker.Lock(0, m_NumOfStripes);
    Unmount();
    //mount disks read-only
    for (unsigned i = 0; i < m_NumOfDisks; i++)
        m_pDisks[i].Mount(false);

    bool Result = true;
    for (unsigned long long S = 0; S < m_NumOfStripes; S++) {
        bool R = m_Engine.VerifyStripe(S, 0);
        if (!R) {
            cerr << "Invalid stripe " << S << endl;
            Result = false;
        };
    };
    if (OldState != msUnmounted)
        Mount(OldState == msReadWrite);
    m_Locker.Unlock(LockID);
    return Result;
};


/** If the requested range does not fit into an integer number of stripe units,
 * read the incomplete ones and extract the required information from them.
 * The remaining data is read via a huge Read call
 * @return the actual number of bytes read, or -1 in case of error
 */
long long CDiskArray::read(tHandle &fd, ///file description, i.e. current position
                           long long Bytes2Read, ///the number of bytes to be read
                           unsigned char *pDest ///destination address
) {
    long long NewPos = fd + Bytes2Read;
    if ((unsigned long long) NewPos > GetCapacity())
        NewPos = GetCapacity();
    Bytes2Read = NewPos - fd;
    if (Bytes2Read < 0)
        //this should never happen
        return -1;
    unsigned long long S = fd / m_StripeUnitSize;
    unsigned Offset = fd % m_StripeUnitSize;
    size_t ThreadID = m_Locker.Lock(fd / m_StripeSize, NewPos / m_StripeSize + ((NewPos % m_StripeSize) ? 1 : 0));
    if (Offset) {
        //partial stripe unit read is necessary
        unsigned char *pTemp = m_pPartialRWBuffer + ThreadID * m_StripeUnitSize;
        if (!Read(S, 1, pTemp, ThreadID)) {
            m_Locker.Unlock(ThreadID);
            return -1;
        };
        unsigned L = m_StripeUnitSize - Offset;
        if (L > Bytes2Read)
            L = (unsigned) Bytes2Read;
        memcpy(pDest, pTemp + Offset, L);
        fd += L;
        pDest += L;
        S++;
    };
    unsigned long long Stripes2Read = (NewPos - fd) / m_StripeUnitSize;
    if (!Read(S, Stripes2Read, pDest, ThreadID)) {
        m_Locker.Unlock(ThreadID);
        return -1;
    };
    S += Stripes2Read;
    pDest += Stripes2Read * m_StripeUnitSize;
    fd += Stripes2Read * m_StripeUnitSize;
    if (fd < NewPos) {
        //partial stripe read is necessary
        unsigned char *pTemp = m_pPartialRWBuffer + ThreadID * m_StripeUnitSize;
        if (!Read(S, 1, pTemp, ThreadID)) {
            m_Locker.Unlock(ThreadID);
            return -1;
        };
        memcpy(pDest, pTemp, (NewPos - fd));
        fd = NewPos;
    };
    m_Locker.Unlock(ThreadID);
    return Bytes2Read;
};


/** Write a number of bytes. If the requested write range does not fit into an
 * integer number of stripe units, the incomplete ones will be read, partially
 * updated and written back
 @return the actual number of bytes read, or -1 in case of error
 */
long long CDiskArray::write(tHandle &fd, ///file description, i.e. current position
                            long long Bytes2Write, ///the number of bytes to be written
                            const unsigned char *pSrc ///source address, must be aligned
) {
    long long NewPos = fd + Bytes2Write;
    if ((unsigned long long) NewPos > GetCapacity())
        NewPos = GetCapacity();
    Bytes2Write = NewPos - fd;
    if (Bytes2Write < 0)
        //this should never happen
        return -1;
    unsigned long long S = fd / m_StripeUnitSize;
    unsigned Offset = fd % m_StripeUnitSize;
    size_t ThreadID = m_Locker.Lock(fd / m_StripeSize, NewPos / m_StripeSize + ((NewPos % m_StripeSize) ? 1 : 0));
    if (Offset) {
        //partial stripe write is necessary
        unsigned char *pTemp = m_pPartialRWBuffer + ThreadID * m_StripeUnitSize;
        if (!Read(S, 1, pTemp, ThreadID)) {
            m_Locker.Unlock(ThreadID);
            return -1;
        };
        unsigned L = m_StripeUnitSize - Offset;
        if (L > Bytes2Write)
            L = (unsigned) Bytes2Write;
        memcpy(pTemp + Offset, pSrc, L);
        if (!Write(S, 1, pTemp, ThreadID)) {
            m_Locker.Unlock(ThreadID);
            return -1;
        };
        fd += L;
        pSrc += L;
        S++;
    };
    unsigned long long Stripes2Write = (NewPos - fd) / m_StripeUnitSize;
    if (!Write(S, Stripes2Write, pSrc, ThreadID)) {
        m_Locker.Unlock(ThreadID);
        return -1;
    };
    S += Stripes2Write;
    pSrc += Stripes2Write * m_StripeUnitSize;
    fd += Stripes2Write * m_StripeUnitSize;
    if (fd < NewPos) {
        //partial stripe write is necessary
        unsigned char *pTemp = m_pPartialRWBuffer + ThreadID * m_StripeUnitSize;
        if (!Read(S, 1, pTemp, ThreadID)) {
            m_Locker.Unlock(ThreadID);
            return -1;
        };
        memcpy(pTemp, pSrc, NewPos - fd);
        if (!Write(S, 1, pTemp, ThreadID)) {
            m_Locker.Unlock(ThreadID);
            return -1;
        };
        fd = NewPos;
    };
    m_Locker.Unlock(ThreadID);
    return Bytes2Write;
};
