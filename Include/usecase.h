/*********************************************************
 * usecase.h  - prototypes for  various testbed usage scenarios
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/
#ifndef USECASE_H
#define USECASE_H

#include "array.h"

/// Initialize a disk array
///@return 0 on success
int InitializeArray(CDiskArray& A);

/// write the sequence of 0,1,2,3,... to the array
/// read it back and validate
///@return 0 on success
int IntegerReadVerify(
    CDiskArray& A, /// the array to be inspected
    unsigned BlocksPerRequest /// the number of stripe units to be read/written simultaneously
);

/// store a given file in the disk array
///@return 0 on success
int StoreFile(
    CDiskArray& A, /// the array to be used
    const char* pFilename /// the name of the file to be stored
);

/// extract file stored in the disk array
///@return 0 on success
int ReadFile(
    CDiskArray& A, /// the array to be used
    const char* pFilename
);

/// check array consistency
///@return 0 on success
int Check(CDiskArray& A /// the array to be checked
);

/// run performance benchmarks
///@return 0 on success
int Benchmark(
    CDiskArray& A, /// the array to be benchmarked
    bool Random, /// true if random read/write is needed, otherwise linear
    unsigned BlockSize, /// size of the data blocks to be read/written
    bool Aligned, /// true if the read-write requests should be aligned to BlockSize multiple
    double WriteRatio, /// the fraction of write requests
    unsigned ThreadCount, /// number of threads to spawn
    unsigned MaxDuration /// maximal benchmark duration (sec)
);

#endif
