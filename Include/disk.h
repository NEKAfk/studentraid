/*********************************************************
 * disk.h  - header file for file-based hard disk drive emulator
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/
#ifndef DISK_H
#define DISK_H

#include "config.h"
#include "sync.h"

#include <stdlib.h>
#include <time.h>

#include <string>

#ifdef USE_MMAP
// include headers for memory mapped files
#ifdef WIN32
#include <windows.h>
#else
#include <sys/mman.h>
#endif
#endif

/// Possible disk state

enum eDiskState {
  dsInvalid, /// The disk file was not properly initialized
  dsOffline, /// Disk is not available
  dsOnline /// The disk is accessible and is assumed to contain correct data
};

/// Possible mount state

enum eMountState {
  msUnmounted, /// the data cannot be accessed
  msRead, /// the data can be read
  msReadWrite /// the data can be read and written
};

// enable memory-mapped files

/** Provides block-based access interface to the hard disk emulated as an ordinary file*/
class CDisk {
  /// name of the underlying file
  const char* m_pFileName;
#ifdef USE_MMAP
  /// pointer to the file mapped to the address space
  unsigned char* m_pMap;
#ifdef WIN32
  /// handle to the file
  HANDLE m_File;
  /// handle to the mapping
  HANDLE m_Mapping;
#else
  /// the descriptor of the underlying file
  int m_File;
#endif

#else
  /// the descriptor of the underlying file
  int m_File;
#endif
  /// disk identifier
  unsigned m_DiskID;
  /// current disk status
  eDiskState m_DiskState;
  /// current mount status
  eMountState m_MountState;
  /// size of one block - the smallest accessible unit
  unsigned m_BlockSize;
  /// number of blocks available
  size_t m_NumOfBlocks;
  /// last write-unmount time
  time_t m_LastUnmount;
  /// array control data stored in the header
  void* m_pArrayData;
  /// size of the array data
  unsigned m_ArrayDataSize;
  /// start of the payload data
  unsigned m_PayloadOffset;
  /*    ///true if the data on disk were modified
      bool m_Dirty;*/
  /// read-write lock
  tCriticalSection m_Lock;

  /// enter a critical section
  void Lock();

  /// leave a critical section
  void Unlock();

  /// write an updated header,
  ///@return true on success
  bool WriteHeader();

  std::ostream Filename();

public:
  /// default constructor. Set to the invalid state
  CDisk();

  /// try to open the file. The parameters on disk will be checked
  ///  with those given to the constructor. The disk state will be set
  ///  to invalid in case of mismatch. Otherwise, it will be set to offline.
  ///  @return true on success

  bool Initialize(
      const char* pFilename, /// the name of the backend file
      unsigned DiskID, /// disk identifier within the array
      unsigned BlockSize, /// the intended block size
      size_t NumOfBlocks, /// number of blocks in the file
      unsigned ArrayDataSize /// size of the disk array configuration structure
  );

  /// this is a wrapper for Initialize()
  CDisk(
      const char* pFilename, /// the name of the backend file
      unsigned DiskID, /// disk identifier within the array
      unsigned BlockSize, /// the intended block size
      size_t NumOfBlocks, /// number of blocks in the file
      unsigned ArrayDataSize /// size of the disk array configuration structure
  );

  /// close the file and deallocate memory
  virtual ~CDisk();

  ///@return the current disk state

  eDiskState GetDiskState() const {
    return m_DiskState;
  };

  /// set the current disk state

  void SetDiskState(eDiskState State) {
    m_DiskState = State;
    // if the disk lost the online status, force unmount
    if (State != dsOnline) {
      m_MountState = msUnmounted;
    }
  };

  ///@return the current mount state

  eMountState GetMountState() const {
    return m_MountState;
  };

  ///@return the size of one payload data block

  unsigned GetBlockSize() const {
    return m_BlockSize;
  };

  ///@return the time of last disk write-unmount

  time_t GetLastUnmountTime() const {
    return m_LastUnmount;
  };

  /// set the array data block. A copy of data will be made
  void SetArrayData(const void* pData, unsigned Size);

  /// get the array data block
  ///@return data block size; 0 if not available

  unsigned GetArrayData(const void*& pData) const {
    pData = m_pArrayData;
    return m_ArrayDataSize;
  };

  /// Initialize the disk. The disk must be in dsOffline or dsInvalid state.
  /// The payload data is filled with zeroes. On success, the disk status is changed to online
  ///@return true on success
  bool ResetDisk();

  /// Try to mount the disk. The disk must be in dsOnline state and not mounted
  ///@return true on success
  bool Mount(bool Write /// true if read/write mount is needed. Otherwise, mount read only
  );

  /// Unmount the disk and set the timestamp. The disk must be in dsOnline state and mounted
  ///@return true on success
  bool Unmount(time_t Timestamp /// the unmount timestamp to be written to the disk if
  );

  /// read a number of payload data blocks. The disk must be mounted
  ///@return true on success
  bool ReadData(
      unsigned long long BlockID, /// the first block to be read
      unsigned NumOfBlocks, /// the number of data blocks to be read
      void* pDest /// destination address. Must have size for at least NumOfBlocks*GetBlockSize() bytes
  );

  /// write a number of payload data blocks, The disk must be read-write mounted
  ///@return true on success
  bool WriteData(
      unsigned long long BlockID, /// start of the destination area
      unsigned NumOfBlocks, /// the number of blocks to be written
      const void* pData /// the data to be written
  );
};

#endif
