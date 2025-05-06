/*********************************************************
 * locker.h  - header file for a range locker 
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/

#ifndef LOCKER_H
#define LOCKER_H

#include <algorithm>
#include "sync.h"


///this class provides thread locking for critical sections given by an 
///integer interval. 

class CRangeLocker {
    ///maximal number of threads which can be executed concurrently
    unsigned m_MaxThreads;

    ///possible lock states
    enum eLockStates {
        lsInvalid, ///the lock is invalid
        lsLocked, ///the lock has been given and a thread is running,
        lsUnlocked
        ///the lock has been released, and the entry is waiting for all relevant threads to process this event
    };

    //locked range

    struct LockedRange {
        ///all entries x with Low<=x<High will be locked
        unsigned long long Low;
        ///all entries x with Low<=x<High will be locked
        unsigned long long High;
        ///the variable used to signal that the lock state has changed
        tCondVariable Condition;
        ///current lock state
        eLockStates State;
        ///the number of threads waiting for it
        volatile unsigned WaitCount;
        ///pointer to the next entry in the active node list
        LockedRange *pNext;
        ///pointer to the previous entry in the active node list
        LockedRange *pPrev;
    };

    ///a pool of locks
    LockedRange *m_pLockPool;
    ///a stack of unused locks
    LockedRange **m_ppFreeLocks;
    ///index of the first free element in the buffer
    ///the value is equal to the current number of active locks
    unsigned m_FreeStackTop;
    ///pointer to a double-linked list of granted locks
    LockedRange *m_pActiveLocks;

    ///used to signal about m_FreeStackTop changes
    tCondVariable m_FreePoolSig;
    ///the global mutex used to protect the internal data structures
    tCriticalSection m_GlobalMutex;

    ///wait for the locked range to become unlocked
    void Wait(LockedRange *pLock ///a lock entry in the list
    );

    ///remove an entry from the active lock list
    void Release(LockedRange *pLock);

public:
    CRangeLocker(unsigned MaxThreads ///maximal number of threads to obtain locks simultaneously
    );

    ~CRangeLocker();

    ///lock the specified range [RangeLow,RangeHigh). The function will wait if
    /// a part of this range is locked by another thread
    ///@return the unique ID of the lock (<m_MaxQueueSize)
    size_t Lock(const unsigned long long RangeLow, ///lower bound
                const unsigned long long RangeHigh ///upper bound
    );

    ///unlock the range
    void Unlock(size_t LockID ///the ID value returned by Lock
    );
};


#endif
