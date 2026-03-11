#define LATTICE 27
#define HIGHORDER

#define TWOCOMPONENT
#define DENSRATIO

#define DOBENCHMARK
#define DOXDMF
#define noWRITEPRESS

#define INTERFACE_INCOMP
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
#define noTAYLORGREEN
#define noUSEGNUPLOT

#ifdef POISEUILLE
#define INTERNAL_OBSTACLES
#define BOUNCE_BACK
#endif
