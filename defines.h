#define LATTICE 27
#define HIGHORDER

#define TWOCOMPONENT
#define DENSRATIO

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
#define noTWOPOISEUILLE
#define noTAYLORGREEN
#define noLAMBTEST
#define noLAPLACE

#define USEGNUPLOT
#define noPRINTPHI

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
#define ELIPSLAMB
#define noMILLER
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
