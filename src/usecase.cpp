/*********************************************************
 * usecase.h  - implementation of  various testbed usage scenarios
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/

#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <iostream>

#ifdef WIN32
#include <process.h>
#endif

#include "arithmetic.h"
#include "misc.h"
#include "sync.h"
#include "usecase.h"

using namespace std;

int InitializeArray(CDiskArray& A) {
  if (!A.Init()) {
    cout << "Array initialization failed\n";
    return 2;
  } else {
    cout << "Array initialization successful\n";
    return 0;
  };
};

/**Generate a sequence of 0,1,2,..., write it to the array,
 * verify the codewords, read back and check correctness
 * @return 0 on success
 * */
int IntegerReadVerify(
    CDiskArray& A, /// the array to be inspected
    unsigned BlocksPerRequest /// the number of stripe units to be read/written simultaneously
) {
  unsigned long long Size = A.GetCapacity();
  unsigned long long CounterSize = Size / sizeof(unsigned);
  unsigned* pData = new unsigned[Size];
  unsigned StripeUnitSize = A.GetStripeUnitSize();
  unsigned long long SizeInUnits = Size / StripeUnitSize;

  if (BlocksPerRequest) {
    CounterSize = (CounterSize / BlocksPerRequest / StripeUnitSize) * BlocksPerRequest * StripeUnitSize;
  }

  if (!A.Mount(true)) {
    cerr << "Array mount failed\n";
    return 3;
  };
  unsigned offset = time(NULL);
  for (unsigned i = 0; i < CounterSize; i++) {
    pData[i] = i + offset;
  }
  unsigned char* pcData = (unsigned char*) pData;
  CDiskArray::tHandle F = A.open();
  double StartTime, StopTime, Dummy;
  GetTimes(StartTime, Dummy, Dummy);
  if (BlocksPerRequest) {
    // CounterSize = (SizeInUnits/BlocksPerRequest)*BlocksPerRequest*StripeUnitSize;
    for (unsigned i = 0; i < (SizeInUnits / BlocksPerRequest) * BlocksPerRequest; i += BlocksPerRequest) {
      if (A.write(F, StripeUnitSize * BlocksPerRequest, pcData + i * StripeUnitSize) !=
          StripeUnitSize * BlocksPerRequest) {
        cerr << "Unit " << i << " write failed\n";
        delete[] pData;
        return 2;
      };
    }
  } else {
    if (A.write(F, Size, pcData) != Size) {
      cerr << "Write failed\n";
      delete[] pData;
      return 2;
    };
  };
  GetTimes(StopTime, Dummy, Dummy);
  cout << "Write throughput " << Size / (StopTime - StartTime) << " bytes/s" << endl;
#ifdef OPERATION_COUNTING
  cout << "Operations per byte: ";
  for (unsigned i = 0; i < opEnd; i++) {
    cout << pOpNames[i] << '(' << double(OPCount[i]) / Size << ") ";
  }
  cout << endl;
  ResetOpCount();
#endif
  StartTime = StopTime;
  if (!A.Check()) {
    cerr << "Array self-check failed\n";
  };
  GetTimes(StopTime, Dummy, Dummy);
  cout << "Check throughput " << Size / (StopTime - StartTime) << " bytes/s" << endl;
#ifdef OPERATION_COUNTING
  cout << "Operations per byte: ";
  for (unsigned i = 0; i < opEnd; i++) {
    cout << pOpNames[i] << '(' << double(OPCount[i]) / Size << ") ";
  }
  cout << endl;
  ResetOpCount();
#endif

  memset(pData, 0xff, Size);
  A.seek(F, 0, SEEK_SET);
  GetTimes(StartTime, Dummy, Dummy);
  if (BlocksPerRequest) {
    // CounterSize = (SizeInUnits/BlocksPerRequest)*BlocksPerRequest*StripeUnitSize;
    for (unsigned i = 0; i < (SizeInUnits / BlocksPerRequest) * BlocksPerRequest; i += BlocksPerRequest) {
      if (A.read(F, StripeUnitSize * BlocksPerRequest, pcData + i * StripeUnitSize) !=
          StripeUnitSize * BlocksPerRequest) {
        cerr << "Unit " << i << " read failed\n";
        delete[] pData;
        return 2;
      };
    }
  } else {
    if (A.read(F, Size, (unsigned char*) pData) != Size) {
      cerr << "Read failed\n";
      delete[] pData;
      return 2;
    };
  };
  GetTimes(StopTime, Dummy, Dummy);
  cout << "Read throughput " << Size / (StopTime - StartTime) << " bytes/s" << endl;
#ifdef OPERATION_COUNTING
  cout << "Operations per byte: ";
  for (unsigned i = 0; i < opEnd; i++) {
    cout << pOpNames[i] << '(' << double(OPCount[i]) / Size << ") ";
  }
  cout << endl;
  ResetOpCount();
#endif
  // make sure read was correct
  for (unsigned i = 0; i < CounterSize; i++) {
    if (pData[i] != i + offset) {
      cerr << "Verify failed at offset " << (i * sizeof(unsigned)) << endl;
      delete[] pData;
      return 3;
    };
  }
  delete[] pData;
  A.Unmount();
  cerr << "Verification successful\n";
  return 0;
};

/// this identifies the properties of the file stored in the array

struct FileHeader {
  /// file size
  off64_t Size;
  /// CRC32 checksum
  unsigned CRC32;
  /// header checksum
  off64_t Checksum;
};

/** store file size, CRC checksum and payload data in the array
 */
int StoreFile(
    CDiskArray& A, /// the array to be used
    const char* pFilename
) {
  if (!A.Mount(true)) {
    cerr << "Array mount failed\n";
    return 3;
  };
  int File = open(pFilename, FILE_IO_OPTIONS);
  if (File < 0) {
    cerr << "Failed to open file " << pFilename << endl;
    return 3;
  };
  off64_t FileSize = lseek64(File, 0, SEEK_END);
  lseek64(File, 0, SEEK_SET);
  unsigned char* pData = new unsigned char[FileSize];
  if (read(File, pData, (unsigned) FileSize) != FileSize) {
    cerr << "Failed to read data from " << pFilename;
    return 3;
  };
  close(File);
  FileHeader Header;
  Header.Size = FileSize;
  InitCRC32();
  Header.CRC32 = 0;
  UpdateCRC32(Header.CRC32, FileSize, pData);
  Header.Checksum = Header.Size ^ Header.CRC32;

  CDiskArray::tHandle F = A.open();
  if (A.write(F, sizeof(FileHeader), (unsigned char*) &Header) != sizeof(FileHeader)) {
    cerr << "Failed to store file header on the array\n";
    return 3;
  };
  double StartTime, StopTime, Dummy;
  GetTimes(StartTime, Dummy, Dummy);

  if (A.write(F, FileSize, pData) != FileSize) {
    cerr << "Failed to store data on the array\n";
    return 3;
  };
  GetTimes(StopTime, Dummy, Dummy);
  delete[] pData;
  cerr << "File stored successfully\n";
#ifdef OPERATION_COUNTING
  cout << "Operations per byte: ";
  for (unsigned i = 0; i < opEnd; i++) {
    cout << pOpNames[i] << '(' << double(OPCount[i]) / FileSize << ") ";
  }
  cout << endl;
  ResetOpCount();
#endif
  return 0;
};

/** read file size, CRC checksum and payload data from the array
 */
int ReadFile(
    CDiskArray& A, /// the array to be used
    const char* pFilename
) {
  if (!A.Mount(false)) {
    cerr << "Array mount failed\n";
    return 3;
  };
  FileHeader Header;
  CDiskArray::tHandle F = A.open();
  if (A.read(F, sizeof(FileHeader), (unsigned char*) &Header) != sizeof(FileHeader)) {
    cerr << "Failed to read file header from the array\n";
    return 3;
  };

  if ((Header.Size ^ Header.CRC32) != Header.Checksum) {
    cerr << "Invalid file header\n";
    return 3;
  };

  unsigned char* pData = new unsigned char[Header.Size];

  double StartTime, StopTime, Dummy;
  GetTimes(StartTime, Dummy, Dummy);
  if (A.read(F, Header.Size, (unsigned char*) pData) != Header.Size) {
    cerr << "Failed to read file data from the array\n";
    return 3;
  };
  GetTimes(StopTime, Dummy, Dummy);
  InitCRC32();
  unsigned CRC32 = 0;
  UpdateCRC32(CRC32, Header.Size, pData);
  // if (CRC32 != Header.CRC32)
  // {
  //     cerr << "File checksum mismatch\n";
  //     return 3;
  // };
  int File = open(pFilename, O_RDWR | O_CREAT | FILE_IO_OPTIONS);
  if (File < 0) {
    cerr << "Error opening file " << pFilename << endl;
    return 3;
  };
  if (write(File, pData, (unsigned) Header.Size) != Header.Size) {
    cerr << "Failed to read data from " << pFilename;
    return 3;
  };
  delete[] pData;
  cerr << "File extracted successfully\n";
#ifdef OPERATION_COUNTING
  cout << "Operations per byte: ";
  for (unsigned i = 0; i < opEnd; i++) {
    cout << pOpNames[i] << '(' << double(OPCount[i]) / Header.Size << ") ";
  }
  cout << endl;
  ResetOpCount();
#endif
  return 0;
};

/// check array consistency
///@return 0 on success

int Check(CDiskArray& A /// the array to be checked
) {
  if (A.Check()) {
    cout << "Array is consistent\n";
    return 0;
  } else {
    cout << "Array is corrupted\n";
    return 3;
  };
};

/// this structure will be used to pass the parameters to the testing thread
/// and get the results back

struct BenchmarkData {
  unsigned ThreadID;
  /// the array to be tested
  CDiskArray* pArray;
  /// true if random read/write is needed, otherwise linear
  bool Random;
  /// size of the data blocks to be read/written
  unsigned BlockSize;
  /// true if the read-write requests should be aligned to BlockSize multiple
  bool Aligned;
  /// the fraction of write requests
  double WriteRatio;
  /// the total number of bytes written
  unsigned long long BytesWritten;
  /// the total number of bytes read
  unsigned long long BytesRead;
  /// the number of IO operations
  unsigned long long IOCount;
};

// a quick and dirty random generator, which provides RAND_MAX=2^32-1

unsigned long long Rand(unsigned long long& RNGState) {
  RNGState = RNGState * 6364136223846793005ull + 1442695040888963407ull;
  return RNGState;
};

// set this to true to cause the testing threads to exit
bool BenchmarkDone = false;

/**Generates a mixture of read and write requests in linear of random order.
 * Passes back the measurement results
 */
#ifdef WIN32
unsigned __stdcall
#else
void *
#endif
    BenchThread(void* pParams /// must be a pointer to BenchmarkData
    ) {
  BenchmarkData& D = *(BenchmarkData*) pParams;
  D.BytesRead = D.BytesWritten = D.IOCount = 0;
  unsigned long long Capacity = D.pArray->GetCapacity();
  unsigned long long MaxSeek = D.pArray->GetCapacity() - D.BlockSize;
  if (D.Aligned) {
    MaxSeek /= D.BlockSize;
  }
  // make sure each thread gets a unique random number sequence
  unsigned long long RNGState = (unsigned long long) pParams;
  // the threshold used to generate random read/write process
  double RWThreshold = pow(2.0, 64) * D.WriteRatio;
  // generate the random data to be read/written
  unsigned char* pData = new unsigned char[D.BlockSize];
  for (unsigned i = 0; i < D.BlockSize; i++) {
    pData[i] = (unsigned char) Rand(RNGState);
  };
  CDiskArray::tHandle F = D.pArray->open();
  if (D.Random) {
    // random access
    while (!BenchmarkDone) {
      unsigned long long Offset = Rand(RNGState) % MaxSeek;
      if (D.Aligned) {
        D.pArray->seek(F, Offset * D.BlockSize, SEEK_SET);
      } else {
        D.pArray->seek(F, Offset, SEEK_SET);
      }
      if (Rand(RNGState) < RWThreshold) {
        // write
        D.pArray->write(F, D.BlockSize, pData);
        D.BytesWritten += D.BlockSize;
      } else {
        D.pArray->read(F, D.BlockSize, pData);
        D.BytesRead += D.BlockSize;
      };
      D.IOCount++;
    };
  } else {
    // linear access
    while (!BenchmarkDone) {
      unsigned long long Offset;
      if (D.Aligned) {
        Offset = 0;
      } else {
        Offset = Rand(RNGState) % D.BlockSize;
      }

      D.pArray->seek(F, Offset, SEEK_SET);
      while (F < (long long) Capacity) {
        if (Rand(RNGState) < RWThreshold) {
          // write
          D.pArray->write(F, D.BlockSize, pData);
          D.BytesWritten += D.BlockSize;
        } else {
          // read
          D.pArray->read(F, D.BlockSize, pData);
          D.BytesRead += D.BlockSize;
        };
        D.IOCount++;
        if (BenchmarkDone) {
          break;
        }
      };
    };
  };
  delete[] pData;
  return 0;
};

/// run performance benchmarks
int Benchmark(
    CDiskArray& A, /// the array to be benchmarked
    bool Random, /// true if random read/write is needed, otherwise linear
    unsigned BlockSize, /// size of the data blocks to be read/written
    bool Aligned, /// true if the read-write requests should be aligned to BlockSize multiple
    double WriteRatio, /// the fraction of write requests
    unsigned ThreadCount, /// number of threads to spawn
    unsigned MaxDuration /// maximal benchmark duration (sec)
) {
  if ((WriteRatio > 1) || (WriteRatio < 0)) {
    cerr << "Invalid write ratio\n";
    return 1;
  };
  if (!A.Mount(true)) {
    cerr << "Array mount failed\n";
    return 2;
  };
  cout << "Running " << ((Random) ? "random " : "linear ") << ((Aligned) ? " aligned" : "non-aligned")
       << " I/O benchmark with " << ThreadCount << " threads,  block size " << BlockSize << " and write ratio "
       << WriteRatio << endl;

  BenchmarkData* pData = new BenchmarkData[ThreadCount];
#ifdef WIN32
  HANDLE* Threads = new HANDLE[ThreadCount];
#else
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_t* Threads = new pthread_t[ThreadCount];
#endif
  double StartTimeU, StartTimeS, StartTimeW;
  GetTimes(StartTimeU, StartTimeS, StartTimeW);
  for (unsigned i = 0; i < ThreadCount; i++) {
    pData[i].pArray = &A;
    pData[i].Random = Random;
    pData[i].BlockSize = BlockSize;
    pData[i].Aligned = Aligned;
    pData[i].WriteRatio = WriteRatio;
    pData[i].ThreadID = i;
    // spawn a benchmarking thread
#ifdef WIN32
    Threads[i] = (HANDLE) _beginthreadex(NULL, 0, BenchThread, pData + i, 0, 0);
#else
    pthread_create(Threads + i, &attr, BenchThread, pData + i);
#endif
  };
#ifdef WIN32
  Sleep(MaxDuration * 1000);
#else
  pthread_attr_destroy(&attr);
  sleep(MaxDuration);
#endif
  BenchmarkDone = true;
  // wait for all threads and collect the statistics
  unsigned long long BytesWritten = 0, BytesRead = 0, IOCount = 0;
  for (unsigned i = 0; i < ThreadCount; i++) {
#ifdef WIN32
    WaitForSingleObject(Threads[i], INFINITE);
    CloseHandle(Threads[i]);
#else
    void* status;
    pthread_join(Threads[i], &status);
#endif
    BytesWritten += pData[i].BytesWritten;
    BytesRead += pData[i].BytesRead;
    IOCount += pData[i].IOCount;
  };
  double StopTimeU, StopTimeS, StopTimeW;
  GetTimes(StopTimeU, StopTimeS, StopTimeW);

  // TimeSpentU - CPU time spent in userspace
  // TimeSpentT - CPU time spent both in userspace and kernel
  // TimeSpentW - wall-clock time
  double TimeSpentU = StopTimeU - StartTimeU;
  double TimeSpentT = StopTimeS - StartTimeS + TimeSpentU;
  double TimeSpentW = StopTimeW - StartTimeW;
  cout << "\nPerformance in terms of userspace, process and wall-clock time:\n"
       << "Read throughput (bytes/s): " << BytesRead / TimeSpentU << '\t' << BytesRead / TimeSpentT << '\t'
       << BytesRead / TimeSpentW << '\n'
       << "Write throughput (bytes/s): " << BytesWritten / TimeSpentU << '\t' << BytesWritten / TimeSpentT << '\t'
       << BytesWritten / TimeSpentW << '\n'
       << "I/O operations per second: " << IOCount / TimeSpentU << '\t' << IOCount / TimeSpentT << '\t'
       << IOCount / TimeSpentW << endl;

  delete[] Threads;
  delete[] pData;
  return 0;
}
