cmake_minimum_required(VERSION 3.29)
project(studentraid)

set(CMAKE_CXX_STANDARD 20)

include_directories(confuse)
include_directories(Include)
add_compile_options(-mssse3)

add_executable(studentraid
    confuse/config_linux.h
    confuse/config_win32.h
    confuse/confuse.c
    confuse/confuse.h
    confuse/lexer.c
    disk/array.cpp
    disk/disk.cpp
    disk/RAIDProcessor.cpp
    Include/arithmetic.h
    Include/array.h
    Include/config.h
    Include/disk.h
    Include/locker.h
    Include/misc.h
    Include/RAID5.h
    Include/RAIDconfig.h
    Include/RAIDProcessor.h
    Include/RS.h
    Include/sync.h
    Include/usecase.h
    RAID/arithmetic.cpp
    RAID/RAID5.cpp
    RAID/RS.cpp
    src/locker.cpp
    src/main.cpp
    src/misc.cpp
    src/usecase.cpp
    Include/RDP.h
    RAID/RDP.cpp)
