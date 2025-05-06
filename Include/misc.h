/*********************************************************
 * misc.h  - miscellaneous definitions for RAID emulator
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/

#ifndef MISC_H
#define MISC_H

#include <exception>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "config.h"

#ifdef WIN32

#include <io.h>
#define lseek64(F,O,S) _lseeki64(F,O,S)
#define lstat64(F,S) _stat64(F,S)
#define ftruncate64(F,S) _chsize_s(F,S)
typedef unsigned __int64 off64_t;
#define O_LARGEFILE 0
#define O_BINARY _O_BINARY
#pragma warning(disable: 4996) /* stop complaining about deprecated POSIX names*/ 
typedef struct _stat64 Stat64;
#define OPEN_FLAGS (_S_IREAD | _S_IWRITE )
#else

#include <unistd.h>
#include <errno.h>

typedef struct stat64 Stat64;

#define OPEN_FLAGS (S_IRUSR|S_IWUSR)
#define O_BINARY 0

#endif

///additional file I/O options
#define FILE_IO_OPTIONS  /*O_DIRECT|*/O_LARGEFILE|O_BINARY/*|O_SYNC*/


///and exception capable of reporting a problem
class Exception {
    char What[1024];

public:
    ///set the error message, using sprintf for message formatting
    Exception(const char *pWhat, ...) {
        va_list Args;
        va_start(Args, pWhat);
        vsprintf(What, pWhat, Args);
        va_end(Args);
    };

    const char *what() const {
        return What;
    };
};

///initialize CRC table
void InitCRC32();

///process data block and update CRC counter
void UpdateCRC32(unsigned &CRC, size_t size, const unsigned char *buf);

///get current time
void GetTimes(double &UserTime, ///user-mode process time
              double &KernelTime, ///kernel-mode process time
              double &WallClockTime ///wall-clock time
);

#ifdef OPERATION_COUNTING
//different operation types
enum eOperations { opXOR, opGFMul, opGFMulAdd, opRead, opWrite, opEnd };

//counter for each operation
extern volatile unsigned long long OPCount[];
//human-readable names for each operation
extern const char *pOpNames[];

//reset operation counters
bool ResetOpCount();


//atomic updating of a counter
#ifdef WIN32
#include <Windows.h>
#define LOCKEDADD(Type,x) InterlockedAdd64((LONG64*)&OPCount[Type],x)
#else
#define LOCKEDADD(Type,x)  __sync_fetch_and_add(&OPCount[Type],x)
#endif
#else
#define LOCKEDADD(Type,x)
#endif

#endif
