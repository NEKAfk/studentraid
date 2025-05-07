/*********************************************************
 * parseconfig.h  - header file for RAID configurator
 *
 * Copyright(C) 2012 Saint-Petersburg State Polytechnic University
 *
 * Developed in the framework of the "Forward error correction for next generation storage systems" project
 *
 * Author: P. Trifonov petert@dcn.ftk.spbstu.ru
 * ********************************************************/

#include "config.h"
#ifdef IMPLEMENT_PARSERS
#undef PARSECONFIG_H
#endif

#ifndef PARSECONFIG_H
#define PARSECONFIG_H


#ifndef IMPLEMENT_PARSERS

//this is a hack to workaround broken variadic macro support in MS Visual C++
#define DECLAREPARAMS__(count,tuple)DECLAREPARAMS_##count tuple
#define DECLARECONFIG__(count,tuple)DECLARECONFIG_##count tuple
#define INITPARAM__(count,tuple)INITPARAM_##count tuple


#define DECLAREPARAM(type,name) type name;


#define DECLAREPARAMS_0()
#define DECLAREPARAMS_1(type,name) DECLAREPARAM(type,name)
#define DECLAREPARAMS_2(type1,name1,type2,name2) DECLAREPARAM(type1,name1) DECLAREPARAMS_1(type2,name2)
#define DECLAREPARAMS_3(type,name,type2,name2,type3,name3) DECLAREPARAM(type,name) DECLAREPARAMS_2(type2,name2,type3,name3)
#define DECLAREPARAMS_4(type,name,type2,name2,type3,name3,type4,name4) DECLAREPARAM(type,name) DECLAREPARAMS_3(type2,name2,type3,name3,type4,name4)
#define DECLAREPARAMS_5(type,name,type2,name2,type3,name3,type4,name4,type5,name5) DECLAREPARAM(type,name) DECLAREPARAMS_4(type2,name2,type3,name3,type4,name4,type5,name5)
#define DECLAREPARAMS_6(type,name,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6) DECLAREPARAM(type,name) DECLAREPARAMS_5(type2,name2,type3,name3,type4,name4,type5,name5,type6,name6)
#define DECLAREPARAMS_7(type,name,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7) DECLAREPARAM(type,name) DECLAREPARAMS_6(type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7)
#define DECLAREPARAMS_8(type,name,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7,type8,name8) DECLAREPARAM(type,name) DECLAREPARAMS_7(type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7,type8,name8)
#define DECLAREPARAMS_9(type,name,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7,type8,name8,type9,name9) DECLAREPARAM(type,name) DECLAREPARAMS_8(type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7,type8,name8,type9,name9)

#define INITPARAM_0(...)
#define INITPARAM_1(cfg,type,name) ,name(cfg_get##type(cfg,#name))
#define INITPARAM_2(cfg,type,name,type2,name2) ,name(cfg_get##type(cfg,#name))INITPARAM_1(cfg,type2,name2)
#define INITPARAM_3(cfg,type,name,type2,name2,type3,name3) ,name(cfg_get##type(cfg,#name))INITPARAM_2(cfg,type2,name2,type3,name3)
#define INITPARAM_4(cfg,type,name,type2,name2,type3,name3,type4,name4) ,name(cfg_get##type(cfg,#name))INITPARAM_3(cfg,type2,name2,type3,name3,type4,name4)
#define INITPARAM_5(cfg,type,name,type2,name2,type3,name3,type4,name4,type5,name5) ,name(cfg_get##type(cfg,#name))INITPARAM_4(cfg,type2,name2,type3,name3,type4,name4,type5,name5)
#define INITPARAM_6(cfg,type,name,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6) ,name(cfg_get##type(cfg,#name))INITPARAM_5(cfg,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6)
#define INITPARAM_7(cfg,type,name,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7) ,name(cfg_get##type(cfg,#name))INITPARAM_6(cfg,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7)
#define INITPARAM_8(cfg,type,name,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7,type8,name8) ,name(cfg_get##type(cfg,#name))INITPARAM_7(cfg,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7,type8,name8)
#define INITPARAM_9(cfg,type,name,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7,type8,name8,type9,name9) ,name(cfg_get##type(cfg,#name))INITPARAM_8(cfg,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7,type8,name8,type9,name9)

#define DECLARECONFIG_0() CFG_END()
#define DECLARECONFIG_1(type,name) CFG_##type(#name,(cfg_t_##type)0,CFGF_NONE), CFG_END()
#define DECLARECONFIG_2(type,name,type2,name2) CFG_##type(#name,(cfg_t_##type)0,CFGF_NONE),DECLARECONFIG_1(type2,name2)
#define DECLARECONFIG_3(type,name,type2,name2,type3,name3) CFG_##type(#name,(cfg_t_##type)0,CFGF_NONE),DECLARECONFIG_2(type2,name2,type3,name3)
#define DECLARECONFIG_4(type,name,type2,name2,type3,name3,type4,name4) CFG_##type(#name,(cfg_t_##type)0,CFGF_NONE),DECLARECONFIG_3(type2,name2,type3,name3,type4,name4)
#define DECLARECONFIG_5(type,name,type2,name2,type3,name3,type4,name4,type5,name5) CFG_##type(#name,(cfg_t_##type)0,CFGF_NONE),DECLARECONFIG_4(type2,name2,type3,name3,type4,name4,type5,name5)
#define DECLARECONFIG_6(type,name,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6) CFG_##type(#name,(cfg_t_##type)0,CFGF_NONE),DECLARECONFIG_5(type2,name2,type3,name3,type4,name4,type5,name5,type6,name6)
#define DECLARECONFIG_7(type,name,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7) CFG_##type(#name,(cfg_t_##type)0,CFGF_NONE),DECLARECONFIG_6(type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7)
#define DECLARECONFIG_8(type,name,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7,type8,name8) CFG_##type(#name,(cfg_t_##type)0,CFGF_NONE),DECLARECONFIG_7(type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7,type8,name8)
#define DECLARECONFIG_9(type,name,type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7,type8,name8,type9,name9) CFG_##type(#name,(cfg_t_##type)0,CFGF_NONE),DECLARECONFIG_8(type2,name2,type3,name3,type4,name4,type5,name5,type6,name6,type7,name7,type8,name8,type9,name9)

#define cfg_t_unsigned unsigned
#define cfg_t_bool cfg_bool_t

#define FOREACH__(N,tuple) FOREACH_##N tuple
//apply a given macro to each of the arguments
#define FOREACH(N,Macro,...) FOREACH_##N(Macro,__VA_ARGS__)
#define FOREACH_1(Macro,X1)  Macro(X1)
#define FOREACH_2(Macro,X1,X2) Macro(X1)   FOREACH_1(Macro,X2)
#define FOREACH_3(Macro,X1,X2,X3) Macro(X1)   FOREACH_2(Macro,X2,X3)
#define FOREACH_4(Macro,X1,X2,X3,X4) Macro(X1)   FOREACH_3(Macro,X2,X3,X4)
#define FOREACH_5(Macro,X1,X2,X3,X4,X5) Macro(X1)   FOREACH_4(Macro,X2,X3,X4,X5)
#define FOREACH_6(Macro,X1,X2,X3,X4,X5,X6) Macro(X1)   FOREACH_5(Macro,X2,X3,X4,X5,X6)
#define FOREACH_7(Macro,X1,X2,X3,X4,X5,X6,X7) Macro(X1)   FOREACH_6(Macro,X2,X3,X4,X5,X6,X7)
#define FOREACH_8(Macro,X1,X2,X3,X4,X5,X6,X7,X8) Macro(X1)   FOREACH_7(Macro,X2,X3,X4,X5,X6,X7,X8)
#define FOREACH_9(Macro,X1,X2,X3,X4,X5,X6,X7,X8,X9) Macro(X1)   FOREACH_8(Macro,X2,X3,X4,X5,X6,X7,X8,X9)


///generates a descriptor for config section
#define PARAMCONFIG(name) CFG_SEC(#name,name##_opts,CFGF_NONE)
///makes a call to the RAID parameter struct constructor, which loads the data from the configuration section
#define READPARAMS(name,variable,config)   name##Params* variable=new name##Params(cfg_getsec(config,#name));


///generates DECLARERAID(...) calls
#define DECLARE(X) DECLARE##X
///generates DEFINERAID(...) calls
#define DEFINE(X) DEFINE##X
///generates PARSERAID(...) calls
#define PARSE(X) PARSE##X
///generates FUNCTORRAID(...)calls
#define FUNCTOR(X) FUNCTOR##X
///generates NAMERAID(...)calls
#define NAME(X) NAME##X


///generates an entry in eRAIDTypes enum
#define DECLARERAID(name,...) rt##name ,
///generates struct with RAID parameters
#define DEFINERAID(name,count,...) RAIDPARAMS(name,count,__VA_ARGS__);
///a stub for parser generator
#define PARSERAID(name,...)
///RAID name
#define NAMERAID(name,...) #name,

///generates a struct with RAID data
#define RAIDPARAMS(name,count,...) \
struct name##Params:public RAIDParams \
{   \
  DECLAREPARAMS__(count,(__VA_ARGS__));\
  name##Params(cfg_t* cfg);\
};\
CFGOPTIONLIST(name,count,__VA_ARGS__) \
CFGCONSTRUCTORIMPL(name,count,__VA_ARGS__)


///a list of RAID configuration descriptors
///NumOfTypes - is the number of RAID types. It must be followed by NumOfTypes entries
///like RAID(N,type1,name1,type2,name2,...)
///where N is the number of parameters of a specific RAID, type,name - their types and names
#define RAIDLIST(NumOfTypes,...) enum eRAIDTypes{ FOREACH__(NumOfTypes,(DECLARE,__VA_ARGS__ rtEnd)) }; \
          FOREACH__(NumOfTypes,(DEFINE,__VA_ARGS__)) \
          PARSECONFIG(NumOfTypes,__VA_ARGS__)

#define LISTFUCTORS(NumOfTypes,...)
#define LISTNAMES(NumOfTypes,...)

#else

#undef RAIDLIST
#undef PARSERAID
#undef LISTFUCTORS
#undef LISTNAMES
///generates configuration file parser calls
#define RAIDLIST(NumOfTypes,...) PARSECONFIG(NumOfTypes,__VA_ARGS__)
///generates a parser for a single configuration
#define PARSERAID(name,...) CRAIDProcessor* Parse##name(cfg_t* cfg) \
{\
     if (!cfg_size(cfg,#name)) return 0;\
     READPARAMS(name,P,cfg);\
     C##name##Processor* pProcessor=new C##name##Processor(P);\
     return pProcessor;\
     };

///a generic RAID parser
typedef CRAIDProcessor* (*tRAIDParser)(cfg_t* cfg);
///generate an entry in parser function table
#define FUNCTORRAID(name,...) Parse##name,
///generate an array of configuration parsers for different RAID types
#define LISTFUCTORS(NumOfTypes,...) tRAIDParser Parsers[]={FOREACH__(NumOfTypes,(FUNCTOR,__VA_ARGS__)) NULL };
///generate a table of possible RAID names
#define LISTNAMES(NumOfTypes,...) const char* ppRAIDNames[]={FOREACH__(NumOfTypes,(NAME,__VA_ARGS__)) NULL };




#endif

//if config header was included, provide the option list implementation
#ifdef _cfg_h_
///config specification  for a given RAID. It includes the common parameters (RAIDParams)
#define CFGOPTIONLIST(name,count,...) cfg_opt_t name##_opts[] ={ CFG_unsigned("Dimension",0,CFGF_NONE),CFG_unsigned("InterleavingOrder",1,CFGF_NONE),  CFG_unsigned("StripeUnitSize",0,CFGF_NONE), DECLARECONFIG__(count,(__VA_ARGS__)) };
///generates a constructor body from a configuration file section
#define CFGCONSTRUCTORIMPL(name,count,...) name##Params::name##Params(cfg_t* cfg):\
            RAIDParams(rt##name,cfg_getint(cfg,"Dimension"),cfg_getint(cfg,"InterleavingOrder"),cfg_getint(cfg,"StripeUnitSize"))INITPARAM__(count,(cfg,__VA_ARGS__))\
  {}
#define CFG_int CFG_INT
#define CFG_bool CFG_BOOL
#define CFG_unsigned CFG_INT
#define cfg_getunsigned cfg_getint

///define parsers for all RAID configurations
#define PARSECONFIG(NumOfTypes,...) FOREACH__(NumOfTypes,(PARSE,__VA_ARGS__)) ; LISTFUCTORS(NumOfTypes,__VA_ARGS__); LISTNAMES(NumOfTypes,__VA_ARGS__);

#else
//no parser declaration is needed
#define CFGOPTIONLIST(name,count,...)
#define CFGCONSTRUCTORIMPL(name,count,...)
#define PARSECONFIG(NumOfTypes,...)

struct cfg_t;

#endif


///this is the list of supported RAID types and their configuration parameters
///each type is declared using RAID(Name,Number of parameters, type1,param1,type2,param2,...)
///the first parameter of RAIDLIST is the total number of RAID types supported
#pragma pack(push)
#pragma pack(1)
RAIDLIST(3,
         RAID(RAID5,0),
         RAID(RS,1,unsigned,Redundancy),
         RAID(RDP,1,unsigned,PrimeNumber)
)
#pragma pack(pop)


#endif
