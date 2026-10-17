#include <sys/resource.h>
#include <sys/time.h>
#include <stdio.h>
/* 
               struct timeval ru_utime; // user CPU time used 
               struct timeval ru_stime; // system CPU time used 
               long   ru_maxrss;        // maximum resident set size
               long   ru_ixrss;         // integral shared memory size 
               long   ru_idrss;         // integral unshared data size 
               long   ru_isrss;         // integral unshared stack size 
               long   ru_minflt;        // page reclaims (soft page faults)
               long   ru_majflt;        // page faults (hard page faults) 
               long   ru_nswap;         // swaps 
               long   ru_inblock;       // block input operations 
               long   ru_oublock;       // block output operations 
               long   ru_msgsnd;        // IPC messages sent 
               long   ru_msgrcv;        // IPC messages received 
               long   ru_nsignals;      // signals received 
               long   ru_nvcsw;         // voluntary context switches 
               long   ru_nivcsw;        // involuntary context switches 
*/

#ifdef __cplusplus
extern "C" {
#endif

/*
 * get_mem: peak resident set size of the calling process, from getrusage(RUSAGE_SELF).
 * Linux reports ru_maxrss in kilobytes, so the returned value ru_maxrss/(1024*1024) is in GiB.
 * Called from Fortran through the bind(C) interface of profiling_m (get_memory, which sums
 * the value over the MPI ranks). When compiled with -DPGI it returns 0.
 */
double get_mem()
{
    struct rusage usg2;       // for memory leaks..
    double mem;
    const double megabyte = 1024 * 1024;
//    long time;
#ifdef PGI
// do nothing
    mem = 0;
#else
    getrusage(RUSAGE_SELF,&usg2);
    mem = ((double)usg2.ru_maxrss) / megabyte;
//    time = usg2.ru_utime.tv_sec;
//    fprintf(stdout," time --> %ld\n ",time);
#endif
    return mem;
}

#ifdef __cplusplus
}
#endif

