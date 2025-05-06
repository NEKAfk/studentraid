/*********************************************************
 * array.h  - header file for RAID emulator
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/
#ifndef ARRAY_H
#define ARRAY_H

#include <string>
#include "disk.h"
#include "RAIDProcessor.h"
#include "locker.h"


///possible states of a disk array

enum eArrayState {
    asUninitialized, ///the array has not yet been initialized
    asFailed, ///the data were lost
    asDegraded, ///some disks have failed
    asNormal ///all disks work fine
};


///disk configuration record

struct DiskConf {
    ///name of the configuration file
    const char *pFileName;
    ///true of the disk is online
    bool Online;
};


///A redundant array of independent disks
///This is a wrapper for various FEC algorithms
///This class provides multi-threaded access to the data

class CDiskArray {
    ///the total number of disks in the array
    unsigned m_NumOfDisks;
    ///the smallest data unit the array can process atomically. This is the same as disk block size
    unsigned m_StripeUnitSize;
    ///the number of payload stripe units within each stripe of a single subarray
    unsigned m_UnitsPerStripePrim;
    ///the number of payload stripe units within each stripe
    unsigned m_UnitsPerStripe;
    ///the total number of stripes within the array
    unsigned long long m_NumOfStripes;
    ///size of each payload stripe in bytes
    unsigned m_StripeSize;
    ///the number of working threads
    unsigned m_NumOfThreads;
    ///current mount state
    eMountState m_MountState;
    ///current array state
    eArrayState m_ArrayState;

    ///underlying disks
    CDisk *m_pDisks;
    ///the computational engine. It must be able to support m_NumOfThreads concurrent calls
    CRAIDProcessor &m_Engine;
    ///temporary buffer for partial stripe unit read/write operations
    unsigned char *m_pPartialRWBuffer;
    ///provides stripe range locking
    CRangeLocker m_Locker;
    ///CRAIDProcessor will directly access m_pDisks
    friend class CRAIDProcessor;

    ///read a number of stripe units. The array must be mounted
    ///@return true on success
    bool Read(unsigned long long StripeUnitID, ///the first stripe unit
              unsigned long long Units2Read, ///the number of stripe units to be read
              unsigned char *pDest, ///destination buffer. Must have size for Units2Read*m_StripeUnitSize bytes
              size_t ThreadID ///the ID of a calling thread obtained from m_Locker
    );

    ///write a number of stripe units. The array must be write-mounted
    ///@return true on success
    bool Write(unsigned long long StripeUnitID, ///the first stripe unit
               unsigned long long Units2Write, ///the number of stripe units to be written
               const unsigned char *pSrc, ///source buffer. Must have size for Units2Write*m_StripeUnitSize bytes
               size_t ThreadID ///the ID of a calling thread obtained from m_Locker
    );

public:
    ///initialize the array. The array parameters
    ///will be extracted from the processor object
    CDiskArray(unsigned NumberOfDisks,
               ///the number of disks to be emulated. This must be not less
                          ///than the length of the array code implemented by Processor
                          ///All extra disks will be ignored
               DiskConf const *pDiskFiles, ///configuration of the emulated disks
               size_t DiskCapacity, ///the capacity of a single disk
               CRAIDProcessor &Processor, ///provides encoding and decoding functionality
               unsigned NumOfThreads ///the number of concurrent processing threads
    );

    virtual ~CDiskArray();

    ///initialize the array. It must be unmounted
    ///@return true on success
    bool Init();

    ///@return current array state

    eArrayState GetState() const {
        return m_ArrayState;
    };
    ///@return RAID type
    eRAIDTypes GetType() const {
        return (eRAIDTypes) m_Engine.GetType();
    };
    ///get the total number of disks
    unsigned GetNumOfDisks() const {
        return m_NumOfDisks;
    };
    ///@return true if the i-th disk is online
    bool IsDiskOnline(unsigned i) const {
        return m_pDisks[i].GetDiskState() == dsOnline;
    };
    ///@return the number of subarrays
    unsigned GetNumOfSubarrays() const {
        return m_Engine.GetInterleavingOrder();
    };

    ///enable data access
    ///@return true on success, false if not mountable
    bool Mount(bool Write ///true if write access is needed
    );

    ///disable data access
    ///@return true on success
    bool Unmount();

    ///check if the array is consistend
    ///@return true on success
    bool Check();

    ///get the payload array capacity

    unsigned long long GetCapacity() const {
        return m_NumOfStripes * m_UnitsPerStripe * m_StripeUnitSize;
    };
    ///@return stripe unit size

    unsigned GetStripeUnitSize() const {
        return m_StripeUnitSize;
    };
    ///the virtual file handle
    typedef long long tHandle;
    ///open a "file" for read and write

    tHandle open() const {
        return 0;
    };
    ///seek to a given position
    ///@return the offset location from the beginning of the file, or -1 in case of error

    long long seek(tHandle &fd, ///current position
                   unsigned long long offset, ///seek offset
                   int whence) {
        switch (whence) {
            case SEEK_SET:
                fd = offset;
                return fd;
            case SEEK_CUR: {
                if (fd + offset > GetCapacity())
                    return -1; //seek beyond the end of drive
                fd += offset;
                return fd;
            };
            case SEEK_END: {
                fd = GetCapacity();
                return fd;
            };
            default:
                return -1;
        };
    };

    ///read the data at a given position, updating it
    ///@return the actual number of bytes read, or -1 in case of error
    long long read(tHandle &fd, ///file description, i.e. current position
                   long long Bytes2Read, ///the number of bytes to be read
                   unsigned char *pDest ///destination address, must be aligned
    );

    ///write the data at a given position, updating it
    ///@return the actual number of bytes written, or -1 in case of error
    long long write(tHandle &fd, ///file description, i.e. current position
                    long long Bytes2Write, ///the number of bytes to be read
                    const unsigned char *pSrc ///source address, must be aligned
    );
};


#endif
