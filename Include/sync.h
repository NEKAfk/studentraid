/*********************************************************
 * locker.h  - header file for synchronization and threading primitives
 * Provides a portability layer between POSIX threads and WinAPI
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/
#ifndef SYNC_H
#define SYNC_H

#ifdef WIN32

// requires at least Windows Vista or Windows Server 2008
#include <Windows.h>

typedef CRITICAL_SECTION tCriticalSection;
typedef CONDITION_VARIABLE tCondVariable;

/// initialize a critical section object
inline bool InitCS(tCriticalSection& CS) {
  InitializeCriticalSection(&CS);
  return true;
};

/// destroy a critical section
inline bool DestroyCS(tCriticalSection& CS) {
  DeleteCriticalSection(&CS);
  return true;
};

/// lock a critical section
inline bool LockCS(tCriticalSection& CS) {
  EnterCriticalSection(&CS);
  return true;
};

// release a critical section
inline bool UnlockCS(tCriticalSection& CS) {
  LeaveCriticalSection(&CS);
  return true;
};

/// initialize a condition object
inline bool InitCond(tCondVariable& Var) {
  InitializeConditionVariable(&Var);
  return true;
};

/// destroy a condition variable
inline bool DestroyCond(tCondVariable& Var) {
  // nothing to do
  return true;
};

/// atomically release the critical section, wait for the condition to be signalled
/// and take back the CS
inline bool CondWait(tCondVariable& C, tCriticalSection& M) {
  SleepConditionVariableCS(&C, &M, INFINITE);
  return true;
}

/// atomically wake everyone waiting for the condition
inline bool CondWakeAll(tCondVariable& C) {
  WakeAllConditionVariable(&C);
  return true;
};

/// wake a thread waiting for the condition
inline bool CondWake(tCondVariable& C) {
  WakeConditionVariable(&C);
  return true;
};

#else
#include <pthread.h>
typedef pthread_cond_t tCondVariable;
typedef pthread_mutex_t tCriticalSection;

/// initialize a critical section object
inline bool InitCS(tCriticalSection& CS) {
  return pthread_mutex_init(&CS, NULL) == 0;
};

/// destroy a critical section
inline bool DestroyCS(tCriticalSection& CS) {
  return pthread_mutex_destroy(&CS) == 0;
};

/// lock a critical section
inline bool LockCS(tCriticalSection& CS) {
  return pthread_mutex_lock(&CS) == 0;
};

/// initialize a condition object
inline bool InitCond(tCondVariable& Var) {
  return pthread_cond_init(&Var, NULL) == 0;
};

/// destroy a condition variable
inline bool DestroyCond(tCondVariable& Var) {
  return pthread_cond_destroy(&Var) == 0;
};

/// atomically release the critical section, wait for the condition to be signalled
/// and take back the CS
inline bool CondWait(tCondVariable& C, tCriticalSection& M) {
  return pthread_cond_wait(&C, &M) == 0;
}

// release a critical section
inline bool UnlockCS(tCriticalSection& CS) {
  return pthread_mutex_unlock(&CS) == 0;
  ;
};

/// atomically wake everyone waiting for the condition
inline bool CondWakeAll(tCondVariable& C) {
  return pthread_cond_broadcast(&C) == 0;
}

/// wake a thread waiting for the condition
inline bool CondWake(tCondVariable& C) {
  return pthread_cond_signal(&C) == 0;
};

#endif

#endif
