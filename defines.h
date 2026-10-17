/* User choices. CMake builds configured with LBFAST_USE_DEFINES_H=OFF
   skip this block and include lbfast_config.h, which CMake generates in
   the build directory from these lines and its cache variables. */
#ifndef LBFAST_CMAKE_CONFIG
#define LATTICE 27
#define HIGHORDER

#define TWOCOMPONENT
#define DENSRATIO
#define CSF

#define DOBENCHMARK
#define DOXDMF
#define noWRITEPRESS

#define noREPULSIVE_FLUX

#define noSMAGORINSKI

#define noINTERNAL_OBSTACLES

#define GPUTILEX 8
#define GPUTILEY 8
#define GPUTILEZ 8

#define ASYNCMPI

#define PRC 8

#define noMIXEDPRC
#define STRPRC 8

#define noVELUNIFORMV
#define noPOISEUILLE
#define noPOISEUILLESTARTREST
#define noTWOPOISEUILLE
#define noTAYLORGREEN
#define noLAMBTEST
#define noLAPLACE
#define noCAPILLARYWAVE

#define USEGNUPLOT
#define noPRINTPHI
#else
#include "lbfast_config.h"
#endif

#ifdef LAPLACE
#define TWOCOMPONENT
#define DENSRATIO
#define PRINTPHI
#define WRITEPRESS
#endif

#ifdef LAMBTEST
#define TWOCOMPONENT
#define DENSRATIO
#define PRINTPHI
#define WRITEPRESS
#endif

#ifdef CAPILLARYWAVE
#define TWOCOMPONENT
#define DENSRATIO
#define PRINTPHI
#endif

#ifdef TWOPOISEUILLE
#define INTERNAL_OBSTACLES
#define TWOCOMPONENT
#define DENSRATIO
#define PRINTPHI
#endif

#ifdef POISEUILLE
#define INTERNAL_OBSTACLES
#undef TWOCOMPONENT
#undef DENSRATIO
#undef INTERFACE_INCOMP
#endif

#ifdef TAYLORGREEN
#undef TWOCOMPONENT
#undef DENSRATIO
#undef INTERFACE_INCOMP
#endif

#ifdef DENSRATIO
#define INTERFACE_INCOMP
#endif

#if !defined(TWOCOMPONENT) && !defined(INTERNAL_OBSTACLES) && !defined(MULTIHIT) && !defined(ELASTIC_FORCE)
#define FUSED_MOMENTS
#endif
