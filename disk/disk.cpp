/*********************************************************
 * disk.cpp  - implementation of the file-based hard disk drive emulator
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/


#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <iostream>
#include "misc.h"
#include "disk.h"
#include "misc.h"


using namespace std;

///file format identifier
#define MAGICNUMBER 0x600DF00D
///disk header version number
#define DISKHEADERVERSION 1


///disk header data

struct DiskHeader {
    ///used to identify correct files; Must be equal to MAGICNUMBER
    unsigned MagicNumber;
    ///version of this structure
    unsigned HeaderVersion;
    ///disk identifier
    unsigned DiskID;
    ///size of one payload data block on disk
    unsigned BlockSize;
    ///number of blocks on disk
    size_t NumOfBlocks;
    ///last write-unmount time
    time_t LastUnmount;
    ///true if the disk was properly initialized and is assumed to contain valid data
    bool Valid;
    ///size of the array control structure
    unsigned ArrayDataSize;
};

///default constructor. Set to the invalid state

CDisk::CDisk() : m_MountState(msUnmounted), m_DiskState(dsInvalid), m_pArrayData(0),
#if defined(USE_MMAP)
                 m_pMap(0),
#ifdef WIN32
        m_File(NULL)
#else
                 m_File(-1)
#endif
#else
        m_File(-1)
#endif


{
    if (!InitCS(m_Lock))
        throw Exception("Failed to initialize disk mutex");
};
///this is a wrapper for Initialize()
CDisk::CDisk(const char *pFilename, ///the name of the backend file
             unsigned DiskID, ///disk identifier within the array
             unsigned BlockSize, ///the intended block size
             size_t NumOfBlocks, ///number of blocks in the file
             unsigned ArrayDataSize ///size of the disk array configuration structure
) {
    if (!InitCS(m_Lock))
        throw Exception("Failed to initialize disk mutex");

    Initialize(pFilename, DiskID, BlockSize, NumOfBlocks, ArrayDataSize);
};

/**try to open the file. The parameters on disk will be checked
 * with those given to the constructor. The disk state will be set
 * to invalid in case of mismatch */
bool CDisk::Initialize(const char *pFilename, ///the name of the backend file
                       unsigned DiskID, ///disk identifier within the array
                       unsigned BlockSize, ///the intended block size
                       size_t NumOfBlocks, ///number of blocks in the file
                       unsigned ArrayDataSize ///size of the disk array configuration structure
) {
    //m_Dirty=false;
    m_MountState = msUnmounted;
    m_pFileName = pFilename;
    m_ArrayDataSize = ArrayDataSize;
    m_BlockSize = BlockSize;
    m_NumOfBlocks = NumOfBlocks;
    m_DiskState = dsInvalid;
    m_DiskID = DiskID;

    m_pArrayData = malloc(ArrayDataSize);
#ifdef USE_MMAP
    m_pMap = NULL;
#ifdef WIN32
    m_File=CreateFile(pFilename,GENERIC_READ | GENERIC_WRITE,FILE_SHARE_READ,NULL,OPEN_EXISTING,0,NULL);
    if (m_File==INVALID_HANDLE_VALUE)
#else
    m_File = open(pFilename, O_RDWR | FILE_IO_OPTIONS); //if the file does not exist, it will not be created
    if (m_File < 0)
#endif
#else
    m_File = open(pFilename, O_RDWR | FILE_IO_OPTIONS); //if the file does not exist, it will not be created
    if (m_File < 0)
#endif
    {
        cerr << "Cannot open file " << pFilename << endl;
        return false;
    };
#ifdef USE_MMAP
    //do the actual memory mapping
#ifdef WIN32
    off64_t FileSize;
    GetFileSizeEx(m_File,(PLARGE_INTEGER )&FileSize);
    m_Mapping=CreateFileMapping(m_File,NULL,PAGE_READWRITE,0,0,NULL);
    if (m_Mapping==NULL)
    {
        cerr<<"Failed to map file "<<pFilename<<" to memory\n";
        return false;
    };
    m_pMap=(unsigned char*) MapViewOfFile(m_Mapping,FILE_MAP_WRITE,0,0,FileSize);
    if (!m_pMap)
    {
        cerr<<"Failed to obtain a memory view of  "<<pFilename<<endl;
        return false;
    };
#else
    //get file size
    off64_t FileSize = lseek64(m_File, 0, SEEK_END);
    lseek64(m_File, 0, SEEK_SET);
    m_pMap = (unsigned char *) mmap(NULL, FileSize,PROT_READ | PROT_WRITE,MAP_SHARED, m_File, 0);
#endif
#else
    //get file size
    off64_t FileSize = lseek64(m_File, 0, SEEK_END);
    lseek64(m_File, 0, SEEK_SET);

#endif

    //load the disk header
#ifdef USE_MMAP
    DiskHeader &Header = *(DiskHeader *) m_pMap;
#else
    DiskHeader Header;
    if (read(m_File, &Header, sizeof ( Header)) != sizeof ( Header))
    {
        cerr << "Failed to read disk header from file " << pFilename << endl;
        return false;
    };
#endif
    //the payload data should start after the disk and array header, but aligned to BlockSize boundary
    m_PayloadOffset = sizeof (Header) + ArrayDataSize;
    m_PayloadOffset = ((m_PayloadOffset / BlockSize) + ((m_PayloadOffset % BlockSize) ? 1 : 0)) * m_BlockSize;
    //check the header validity
    if ((Header.MagicNumber != MAGICNUMBER) ||
        (Header.HeaderVersion != DISKHEADERVERSION)) {
        cerr << "Invalid disk header for disk " << pFilename << endl;
        return false;
    };
    //check the disk validity

    if ((Header.BlockSize != BlockSize) || (Header.NumOfBlocks != NumOfBlocks)) {
        cerr << "Disk configuration does not match array configuration for disk " << pFilename << endl;
        return false;
    };
    if (FileSize != m_PayloadOffset + NumOfBlocks * BlockSize) {
        cerr << "File size does not match header data in " << pFilename << endl;
        return false;
    };
    if (Header.ArrayDataSize != ArrayDataSize) {
        cerr << "Array configuration header size mismatch in " << pFilename << endl;
        return false;
    };
    if (Header.DiskID != m_DiskID) {
        cerr << "Disk ID mismatch in " << pFilename << endl;
        return false;
    };
    /*    if (Header.Dirty)
        {
            cerr<<"Disk "<<pFilename<<" was not properly unmounted\n";
            return false;
        };*/
    //load array configuration
#ifdef USE_MMAP
    memcpy(m_pArrayData, m_pMap + sizeof(Header), m_ArrayDataSize);
#else
    if (read(m_File, m_pArrayData, m_ArrayDataSize) != m_ArrayDataSize)
    {
        cerr << "Failed to read array configuration in " << pFilename << endl;
        return false;
    };
#endif
    if (Header.Valid) {
        m_DiskState = dsOffline;
        m_LastUnmount = Header.LastUnmount;
    } else m_DiskState = dsInvalid; //the disk was not properly initialized before

    return true;
};


///close the file and deallocate memory

CDisk::~CDisk() {
    if (m_MountState == msReadWrite) {
        cerr << "Warning, file " << m_pFileName << " was not properly unmounted\n";
    };
#ifdef USE_MMAP
#ifdef WIN32
    if(m_pMap)
    {
        UnmapViewOfFile(m_pMap);
        CloseHandle(m_Mapping);
        CloseHandle(m_File);
    };
#else
    munmap(m_pMap, m_PayloadOffset + m_NumOfBlocks * m_BlockSize);
#endif
#else
    if (m_File >= 0)
        close(m_File);
#endif
    free(m_pArrayData);
    DestroyCS(m_Lock);
};

/**
 Lock the mutex
 */
void CDisk::Lock() {
    LockCS(m_Lock);
};

/**
 * Unlock the mutex
 */
void CDisk::Unlock() {
    UnlockCS(m_Lock);
}

///Try to mount the disk. The disk must be in dsOnline state and not mounted
///@return true on success

bool CDisk::Mount(bool Write) {
    if ((m_DiskState != dsOnline) && (m_MountState != msUnmounted))
        return false;
    m_MountState = (Write) ? msReadWrite : msRead;
    return true;
};
///set the array data block. A copy of data will be made

void CDisk::SetArrayData(const void *pData, unsigned Size) {
    if (m_ArrayDataSize != Size)
        m_pArrayData = realloc(m_pArrayData, Size);
    m_ArrayDataSize = Size;
    memcpy(m_pArrayData, pData, Size);
    m_PayloadOffset = sizeof(DiskHeader) + m_ArrayDataSize;
    m_PayloadOffset = ((m_PayloadOffset / m_BlockSize) + ((m_PayloadOffset % m_BlockSize) ? 1 : 0)) * m_BlockSize;
};


///write an updated header

bool CDisk::WriteHeader() {
    DiskHeader Header = {
        MAGICNUMBER, DISKHEADERVERSION, m_DiskID, m_BlockSize, m_NumOfBlocks, m_LastUnmount,
        m_DiskState == dsOnline, //the disk is assumed to be valid only if it has been taken online
        m_ArrayDataSize
    };
#ifdef USE_MMAP
    //write disk header
    memcpy(m_pMap, &Header, sizeof(Header));
    //write array header
    memcpy(m_pMap + sizeof(Header), m_pArrayData, m_ArrayDataSize);
#else
    lseek64(m_File, 0, SEEK_SET);
    if (write(m_File, &Header, sizeof ( Header)) != sizeof ( Header))
    {
        cerr << "Failed to update disk header for " << m_pFileName << endl;
        m_DiskState = dsInvalid;
        return false;
    };
    //write array header
    if (write(m_File, m_pArrayData, m_ArrayDataSize) != m_ArrayDataSize)
    {
        cerr << "Failed to update array data for " << m_pFileName << endl;
        m_DiskState = dsInvalid;
        return false;
    };
#endif
    return true;
};


///Unmount the disk. The disk must be in dsOnline state and mounted
///@return true on success

bool CDisk::Unmount(time_t Timestamp) {
    switch (m_MountState) {
        case msRead:
            m_MountState = msUnmounted;
            return true;
        case msReadWrite: {
            Lock();
            m_MountState = msUnmounted;
            //update disk header
            m_LastUnmount = Timestamp;
            bool Result = WriteHeader();
            Unlock();
            return Result;
        };
        default:
            return false;
    };
};

///Initialize the disk. The disk must be in dsOffline or dsInvalid state.
///The payload data is filled with zeroes. On success, the disk status is changed to online
///@return true on success

bool CDisk::ResetDisk() {
    if (m_DiskState == dsOnline)
        return false;
    Lock();

#ifdef USE_MMAP
    //close the existing map, intialize file, and re-establish the map
#ifdef WIN32
    if (m_pMap)
    {
        UnmapViewOfFile(m_pMap);
        CloseHandle(m_Mapping);
        CloseHandle(m_File);
    };
    m_File=CreateFile(m_pFileName,GENERIC_READ | GENERIC_WRITE,FILE_SHARE_READ,NULL,CREATE_ALWAYS,0,NULL);
    if (m_File==INVALID_HANDLE_VALUE)
    {
         cerr << "Failed to create file " << m_pFileName << endl;
         return false;
    };
    //resize the file
    unsigned long long NewSize=m_PayloadOffset+m_NumOfBlocks * m_BlockSize;
    SetFilePointerEx(m_File,*(_LARGE_INTEGER*)&NewSize,NULL,FILE_BEGIN);
    SetEndOfFile(m_File);
    //do the mapping
    m_Mapping=CreateFileMapping(m_File,NULL,PAGE_READWRITE,0,0,NULL);
    if (m_Mapping==NULL)
    {
        cerr<<"Failed to map file "<<m_pFileName <<" to memory\n";
        return false;
    };
    m_pMap=(unsigned char*) MapViewOfFile(m_Mapping,FILE_MAP_WRITE,0,0,NewSize);
    if (!m_pMap)
    {
        cerr<<"Failed to obtain a memory view of  "<<m_pFileName<<endl;
        return false;
    };

#else
    if (m_pMap)
        munmap(m_pMap, m_PayloadOffset + m_NumOfBlocks * m_BlockSize);
    if (m_File < 0) {
        //try to create the files
        m_File = open(m_pFileName, O_RDWR | O_CREAT | FILE_IO_OPTIONS, OPEN_FLAGS);
        if (m_File < 0) {
            cerr << "Failed to create file " << m_pFileName << endl;
            Unlock();
            return false;
        };
    };
    //resize the file appropriately
    ftruncate64(m_File, 0);
    size_t TargetSize = m_PayloadOffset + m_NumOfBlocks * m_BlockSize;
    //this will fill the file with zeroes
    if (ftruncate64(m_File, TargetSize)) {
        Unlock();
        return false;
    };
    m_pMap = (unsigned char *) mmap(NULL, TargetSize,PROT_READ | PROT_WRITE,MAP_SHARED, m_File, 0);
    if (!m_pMap) {
        Unlock();
        return false;
    };
#endif
#else
    //we are going to rebuild the file from scratch
    m_DiskState = dsInvalid;
    if (m_File < 0)
    {
        //try to create the files
        m_File = open(m_pFileName, O_RDWR | O_CREAT | FILE_IO_OPTIONS, OPEN_FLAGS);
        if (m_File < 0)
        {
            cerr << "Failed to create file " << m_pFileName << endl;
            Unlock();
            return false;
        };
    };
    //resize the file appropriately
    ftruncate64(m_File, 0);
    size_t TargetSize = m_PayloadOffset + m_NumOfBlocks*m_BlockSize;
    //this will fill the file with zeroes
    if (ftruncate64(m_File, TargetSize))
    {
        Unlock();
        return false;
    };
#endif
    m_LastUnmount = 0;
    //take the disk online
    m_DiskState = dsOnline;
    //write updated header
    if (!WriteHeader()) {
        //something failed on the way
        m_DiskState = dsInvalid;
        Unlock();
        return false;
    };
    Unlock();
    return true;
};

volatile unsigned long long ReadDataCount = 0;

///read a number of payload data blocks. The disk must be mounted
///@return true on success

bool CDisk::ReadData(unsigned long long BlockID, ///the first block to be read
                     unsigned NumOfBlocks, ///the number of data blocks to be read
                     void *pDest ///destination address. Must have size for at least NumOfBlocks*GetBlockSize() bytes
) {
    if (m_MountState == msUnmounted) //invalid disk access
        return false;
    if (BlockID + NumOfBlocks > m_NumOfBlocks) //invalid read request
        return false;
    LOCKEDADD(opRead, NumOfBlocks*m_BlockSize);
#ifdef USE_MMAP
    memcpy(pDest, m_pMap + m_PayloadOffset + BlockID * m_BlockSize, NumOfBlocks * m_BlockSize);
    return true;
#else
    Lock();
    //seek to the data position
    off64_t Pos = m_PayloadOffset + BlockID*m_BlockSize;
    if (lseek64(m_File, Pos, SEEK_SET) != Pos)
    {
        //something is wrong with the disk
        SetDiskState(dsInvalid);
        Unlock();
        return false;
    };
    unsigned DataSize = NumOfBlocks*m_BlockSize;
    if (read(m_File, pDest, DataSize) != DataSize)
    {
        //something is wrong with the disk
        cerr << "Read error " << strerror(errno) << " while reading from disk " << m_pFileName << endl;
        _ASSERT(_CrtCheckMemory());
        _ASSERT(false);
        SetDiskState(dsInvalid);
        Unlock();
        return false;
    }
    else
    {
        Unlock();
        return true;
    };
#endif
};


///write a number of payload data blocks, The disk must be read-write mounted
///@return true on success

bool CDisk::WriteData(unsigned long long BlockID, ///start of the destination area
                      unsigned NumOfBlocks, ///the number of blocks to be written
                      const void *pData ///the data to be written
) {
    if (m_MountState != msReadWrite) //invalid disk access
        return false;
    if (BlockID + NumOfBlocks > m_NumOfBlocks) //invalid write request
        return false;
    LOCKEDADD(opWrite, NumOfBlocks*m_BlockSize);
#ifdef USE_MMAP
    memcpy(m_pMap + m_PayloadOffset + BlockID * m_BlockSize, pData, NumOfBlocks * m_BlockSize);
    return true;
#else
    Lock();
    //seek to the data position
    off64_t Pos = m_PayloadOffset + BlockID*m_BlockSize;
    if (lseek64(m_File, Pos, SEEK_SET) != Pos)
    {
        //something is wrong with the disk
        SetDiskState(dsInvalid);
        Unlock();
        return false;
    };
    unsigned DataSize = NumOfBlocks*m_BlockSize;
    if (write(m_File, pData, DataSize) != DataSize)
    {
        //something is wrong with the disk
        cerr << "Write error " << strerror(errno) << " while reading from disk " << m_pFileName << endl;
        _ASSERT(false);
        SetDiskState(dsInvalid);
        Unlock();
        return false;
    }
    else
    {
        Unlock();
        return true;
    };
#endif
};
