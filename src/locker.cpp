/*********************************************************
 * locker.cpp  - implementation for a range locker 
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/

#include "locker.h"
#include "misc.h"


/**
 * Allocate locking structures, initialize the mutexes
 */
CRangeLocker::CRangeLocker(unsigned MaxThreads ///maximal number of threads to obtain locks simultaneously
): m_MaxThreads(MaxThreads),
   m_FreeStackTop(0), m_pActiveLocks(0) {
    m_pLockPool = new LockedRange[MaxThreads];
    m_ppFreeLocks = new LockedRange *[MaxThreads];

    if (!InitCS(m_GlobalMutex))
        throw Exception("Global mutex initialization failed");
    if (!InitCond(m_FreePoolSig))
        throw Exception("Buffer condition initialization failed");
    for (unsigned i = 0; i < MaxThreads; i++) {
        if (!InitCond(m_pLockPool[i].Condition))
            throw Exception("Range condition initialization failed");
        m_ppFreeLocks[i] = m_pLockPool + i;
        m_pLockPool[i].State = lsInvalid;
    };
};

CRangeLocker::~CRangeLocker() {
    DestroyCS(m_GlobalMutex);
    DestroyCond(m_FreePoolSig);
    for (unsigned i = 0; i < m_MaxThreads; i++) {
        DestroyCond(m_pLockPool[i].Condition);
    };
    delete[]m_pLockPool;
    delete[]m_ppFreeLocks;
};


/**
 * Increment wait counter, wait until the unlock happens, decrement wait counter, and, if
 * it becomes zero, release the entry
 */
void CRangeLocker::Wait(LockedRange *pLock ///a lock entry in the list
) {
    pLock->WaitCount++;
    while (pLock->State != lsUnlocked)
        CondWait(pLock->Condition, m_GlobalMutex);
    pLock->WaitCount--;
    if ((pLock->State == lsUnlocked) && !pLock->WaitCount)
        Release(pLock);
};


/** 1. Lock the global data structures
    2. wait for some free space in the queue
    3. Insert an entry into the list of active locks
    4. Search the list for overlapping ranges
    5. If a locked overlapping range is found, wait for it and go to 4
    6. Otherwise, grant the lock
    */
size_t CRangeLocker::Lock(const unsigned long long RangeLow, ///lower bound
                          const unsigned long long RangeHigh ///upper bound
) {
    LockCS(m_GlobalMutex);
    //make sure there are not so many active threads
    while (m_FreeStackTop >= m_MaxThreads)
        CondWait(m_FreePoolSig, m_GlobalMutex);
    bool Block = true;
    //check if we have to wait for someone
    while (Block) {
        Block = false;
        LockedRange *pCurRange = m_pActiveLocks;
        while (pCurRange) {
            if (pCurRange->State == lsLocked) {
                //check if we intersect with this range
                if ((RangeHigh > pCurRange->Low) && (RangeLow < pCurRange->High)) {
                    Block = true;
                    Wait(pCurRange);
                    break; //we will need to re-inspect the list from the beginning
                };
            };
            pCurRange = pCurRange->pNext;
        };
    };
    //no conflicts with active locks, grant one
    LockedRange *pRange = m_ppFreeLocks[m_FreeStackTop++];
    pRange->Low = RangeLow;
    pRange->High = RangeHigh;
    pRange->WaitCount = 0;
    //insert it into the beginning of the list of locks
    if (m_pActiveLocks)
        m_pActiveLocks->pPrev = pRange;
    pRange->pNext = m_pActiveLocks;
    m_pActiveLocks = pRange;
    pRange->pPrev = NULL;
    pRange->State = lsLocked;
    UnlockCS(m_GlobalMutex);
    return pRange - m_pLockPool;
};

/** Notify all threads about unlock and remove the element from the list
of active locks*/
void CRangeLocker::Unlock(size_t LockID ///the ID value returned by Lock
) {
    LockCS(m_GlobalMutex);
    //cerr<<"Relese "<<ThreadID<<endl;
    LockedRange &Lock = m_pLockPool[LockID];
    Lock.State = lsUnlocked;
    //remove it from the list of active entries
    if (Lock.pNext)
        Lock.pNext->pPrev = Lock.pPrev;
    if (Lock.pPrev)
        Lock.pPrev->pNext = Lock.pNext;
    else
        //this was the first element in the list
        m_pActiveLocks = Lock.pNext;
    //notify all threads about the unlock
    CondWakeAll(m_pLockPool[LockID].Condition);
    if (!m_pLockPool[LockID].WaitCount)
        //nobody is waiting for it
        Release(m_pLockPool + LockID);
    UnlockCS(m_GlobalMutex);
}

/*
 * remove an entry from the active lock list,
 * put the entry into the list of free ones
 * and notify somebody about this
 */
void CRangeLocker::Release(LockedRange *pLock) {
    pLock->State = lsInvalid;
    m_ppFreeLocks[--m_FreeStackTop] = pLock;
    CondWake(m_FreePoolSig);
}
