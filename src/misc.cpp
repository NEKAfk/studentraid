/*********************************************************
 * misc.cpp  - implementation of various auxiliary functions for the RAID testbed
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/
#include "misc.h"

#include <string.h>
#include <time.h>

#include <iostream>

#ifdef WIN32
#include <Windows.h>
#else
#include <sys/times.h>
#endif

/// auxiliary table for CRC computation
unsigned tb32[256];

#define POLY_32 0xEDB88320ul

/// initialize CRC table
void InitCRC32() {
  unsigned i, j;
  unsigned crc, v;

  for (crc = i = 0; i < 256; tb32[i++] = crc) {
    for (crc = i, j = 0; j < 8; j++) {
      v = -int(crc & 1), crc >>= 1, crc ^= v & POLY_32;
    }
  }
};

/// process data block and update CRC counter
void UpdateCRC32(unsigned& CRC, size_t size, const unsigned char* buf) {
  unsigned i;
  for (i = 0; i < size; i++) {
    CRC ^= *(unsigned char*) (buf + i);
    CRC = (CRC >> 8) ^ tb32[(unsigned char) CRC];
  }
};

/// get current time
void GetTimes(
    double& UserTime, /// user-mode process time
    double& KernelTime, /// kernel-mode process time
    double& WallClockTime /// wall-clock time
) {
#ifdef WIN32
  static HANDLE hProcess = NULL;
  if (!hProcess) {
    hProcess = GetCurrentProcess();
  }
  FILETIME ftDummy, ftUtime, ftKTime;
  GetProcessTimes(hProcess, &ftDummy, &ftDummy, &ftKTime, &ftUtime);
  UserTime = *(__int64*) &ftUtime * 100E-9;
  KernelTime = *(__int64*) &ftKTime * 100E-9;
#else
  tms PTime;
  times(&PTime);
  UserTime = double(PTime.tms_utime) / sysconf(_SC_CLK_TCK);
  KernelTime = double(PTime.tms_stime) / sysconf(_SC_CLK_TCK);
#endif
  time_t WTime;
  time(&WTime);
  WallClockTime = WTime;
};

// counter for each operation
volatile unsigned long long OPCount[opEnd] = {0, 0, 0, 0, 0};
// human-readable names for each operation
const char* pOpNames[opEnd] = {"XOR", "Multiply", "Multiply-XOR", "Read", "Write"};

// reset operation counters
bool ResetOpCount() {
  memset((void*) OPCount, 0, sizeof(OPCount));
  return true;
};

// force counter reset at program startup
bool Reset = ResetOpCount();
