/*********************************************************
* main.cpp  - main file for the RAID algorithms testbed
*
* Copyright(C) 2012 Saint-Petersburg State Polytechnic University
*
* Developed in the framework of the "Forward error correction for next generation storage systems" project
*
* Author: P. Trifonov petert@dcn.ftk.spbstu.ru
* ********************************************************/

#include <iostream>
#include <stdlib.h>
#include <string>
#include <string.h>


//suppress warnings about improper use of const char* by confuse 
#pragma GCC diagnostic ignored "-Wwrite-strings"

#include "confuse.h"
#include "misc.h"
#include "usecase.h"
#include "array.h"
#include "RAID5.h"
#include "RS.h"
#include "RDP.h"

//global configuration options

//this is needed to generate config parsing code for all RAID types
#define IMPLEMENT_PARSERS
#include "RAIDconfig.h"

using namespace std;

void Usage() {
    cerr << "Usage: testbed ConfigFile Mode [Options]\n"
            "\tSupported modes (with options):\n"
            "\t\t i  initialize disk array \n"
            "\t\t v  integer write-read-verify cycle ( BlocksPerRequest (0 for the highest possible) )  \n"
            "\t\t s  store a file on the array ( FileName )  \n"
            "\t\t g  get a file from the array ( FileName )  \n"
            "\t\t c  check array consistency\n"
            "\t\t b  run performance benchmarks ( l|r a|n WriteRatio BlockSize ThreadCount Duration )\n"
            "\t\t\t Access mode: l - linear, r - random\n"
            "\t\t\t Access type: a - BlockSize aligned, n - non-aligned\n ";
};

/**Report a configuration file problem
* */
static void conf_error(cfg_t *cfg, const char *fmt, va_list ap) {
    int ret;
    static char buf[1024];

    ret = vsnprintf(buf, sizeof (buf), fmt, ap);
    if (ret < 0 || ret >= sizeof (buf)) {
        cerr << "could not print error message\n";
    } else {
        cerr << buf << endl;
    }
}

///disk configuration record
cfg_opt_t disk_opts[] = {
    CFG_STR("file", NULL, CFGF_NONE),
    CFG_BOOL("online", cfg_true, CFGF_NONE),
    CFG_END()
};

///configuration file format decriptor
cfg_opt_t opts[] = {
    CFG_INT("DiskCapacity", 1024, CFGF_NONE),
    CFG_INT("MaxConcurrentThreads", 4, CFGF_NONE),
    CFG_STR("RAIDType", NULL, CFGF_NONE),
    CFG_SEC("disk", disk_opts, CFGF_MULTI),
    //all RAID types should be listed here
    PARAMCONFIG(RAID5),
    PARAMCONFIG(RS),
    PARAMCONFIG(RDP),
    CFG_END()
};
char *pArrayStates[] = {"Uninitialized", "Failed", "Degraded", "Normal "};

CRAIDProcessor *GetProcessor(cfg_t *pConfig) {
    const char *pType = cfg_getstr(pConfig, "RAIDType");
    if (!pType) {
        cerr << "RAID type was not specified in the configuration file\n";
        return 0;
    };

    for (unsigned i = 0; i < rtEnd; i++) {
        if (strcmp(pType, ppRAIDNames[i]) == 0) {
            CRAIDProcessor *p = Parsers[i](pConfig);
            return p;
        }
    };
    cerr << "Unknown RAID type " << pType << endl;
    return 0;
};

int main(int argc, char **argv) {
#ifdef _M_IX86
    cerr<<"WARNING! THIS PROGRAM MAY WORK INCORRECTLY IF COMPILED IN 32-BIT MODE!\n";
#endif

    if (argc < 3) {
        Usage();
        return 1;
    };
    cfg_t *cfg, *cfg_disk;
    cfg = cfg_init(opts, CFGF_NONE);
    cfg_set_error_function(cfg, conf_error);

    if (cfg_parse(cfg, argv[1]) != CFG_SUCCESS) {
        cerr << "Error parsing configuration file " << argv[1] << endl;
        return 1;
    };
    unsigned DiskCapacity = cfg_getint(cfg, "DiskCapacity");
    unsigned NumOfDisks = cfg_size(cfg, "disk");
    unsigned MaxConcurrentThreads = cfg_getint(cfg, "MaxConcurrentThreads");
    if (!NumOfDisks) {
        cerr << "No disk configuration found in the configuration file " << argv[1] << endl;
        return 1;
    };

    try {
        DiskConf *pDisks = new DiskConf[NumOfDisks];
        for (unsigned i = 0; i < NumOfDisks; i++) {
            cfg_disk = cfg_getnsec(cfg, "disk", i);
            pDisks[i].pFileName = cfg_getstr(cfg_disk, "file");
            pDisks[i].Online = cfg_getbool(cfg_disk, "online") > 0;
        };


        CRAIDProcessor *pProcessor = GetProcessor(cfg);
        if (!pProcessor) {
            cerr << "Failed to initialize RAID processor\n";
            return 1;
        };
        CDiskArray Array(NumOfDisks, pDisks, DiskCapacity, *pProcessor, MaxConcurrentThreads);
        cout << "Array type is " << ppRAIDNames[Array.GetType()] << '*' << Array.GetNumOfSubarrays() << endl;
        cout << "Array state is " << pArrayStates[Array.GetState()] << endl;
        cout << "Disk status ";
        for (unsigned i = 0; i < Array.GetNumOfDisks(); i++)
            cout << Array.IsDiskOnline(i);
        cout << endl;


        int Result = 1;
        char c = argv[2][0];
        switch (c) {
            case 'i':
                Result = InitializeArray(Array);
                break;
            case 'v':
                if (argc == 4) {
                    unsigned BlocksPerRequest = atoi(argv[3]);
                    Result = IntegerReadVerify(Array, BlocksPerRequest);
                } else Usage();
                break;
            case 's':
                if (argc == 4) {
                    Result = StoreFile(Array, argv[3]);
                } else Usage();
                break;
            case 'g':
                if (argc == 4) {
                    Result = ReadFile(Array, argv[3]);
                } else Usage();
                break;
            case 'c':
                Result = Check(Array);
                break;
            case 'b': {
                if (argc == 9) {
                    //run performance benchmarks
                    bool Random, Aligned;
                    switch (argv[3][0]) {
                        case 'l': Random = false;
                            break;
                        case 'r': Random = true;
                            break;
                        default: {
                            cerr << "Access mode can be eight l (linear) or r(random)\n";
                            return 2;
                        }
                    };
                    switch (argv[4][0]) {
                        case 'a': Aligned = true;
                            break;
                        case 'n': Aligned = false;
                            break;
                        default: {
                            cerr << "Access type can be eight a (aligned) or n(non-aligned)\n";
                            return 2;
                        }
                    };
                    double WriteRatio = atof(argv[5]);
                    unsigned BlockSize = atoi(argv[6]);
                    unsigned ThreadCount = atoi(argv[7]);
                    unsigned MaxTime = atoi(argv[8]);
                    Result = Benchmark(Array, Random, BlockSize, Aligned, WriteRatio, ThreadCount, MaxTime);
                } else Usage();
                break;
            }
            default: {
                Usage();
                cfg_free(cfg);
                return 1;
            };
        };
        cfg_free(cfg);
        delete pProcessor;
        delete[]pDisks;
        return Result;
    } catch (const Exception &ex) {
        cerr << ex.what() << endl;
        cfg_free(cfg);
        return 1;
    };

    return 0;
}
