#include <sys/sysinfo.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <stdio.h>
/* 
struct sysinfo {
       long uptime;              Seconds since boot  
       unsigned long loads[3];   1, 5, and 15 minute load averages  
       unsigned long totalram;   Total usable main memory size  
       unsigned long freeram;    Available memory size  
       unsigned long sharedram;  Amount of shared memory  
       unsigned long bufferram;  Memory used by buffers  
       unsigned long totalswap;  Total swap space size  
       unsigned long freeswap;   swap space still available  
       unsigned short procs;     Number of current processes  
       unsigned long totalhigh;  Total high memory size  
       unsigned long freehigh;   Available high memory size  
       unsigned int mem_unit;    Memory unit size in bytes  
       char _f[20-2*sizeof(long)-sizeof(int)];  Padding for libc5  
   };
*/

/*
#include <linux/kernel.h>
#include <linux/sys.h>
#include <stdio.h>
#include <sys/sysinfo.h>

int main ()
{
  // Conversion constants. 
 const long minute = 60;
 const long hour = minute * 60;
 const long day = hour * 24;
 const double megabyte = 1024 * 1024;
  // Obtain system statistics. 
 struct sysinfo si;
 sysinfo (&si);
  // Summarize interesting values. 
 printf ("system uptime : %ld days, %ld:%02ld:%02ld\n", 
     si.uptime / day, (si.uptime % day) / hour, 
     (si.uptime % hour) / minute, si.uptime % minute);
 printf ("total RAM   : %5.1f MB\n", si.totalram / megabyte);
 printf ("free RAM   : %5.1f MB\n", si.freeram / megabyte);
 printf ("process count : %d\n", si.procs);
 return 0;
}
*/

#ifdef __cplusplus
extern "C" {
#endif

double get_totalram()
{
    struct sysinfo usg2;       // for memory leaks..
  //unsigned long int myram;
    double mem;
    const double megabyte = 1024 * 1024;
//    long time;
#ifdef PGI
// do nothing
    mem = 0;
#else
    sysinfo (&usg2);
    mem = ((double)usg2.totalram) / megabyte;
//    time = usg2.ru_utime.tv_sec;
//    fprintf(stdout," time --> %ld\n ",time);
#endif
    return mem;
}

double get_freeram()
{
    struct sysinfo usg2;       // for memory leaks..
  //unsigned long int myram;
    double mem;
    const double megabyte = 1024 * 1024;
//    long time;
#ifdef PGI
// do nothing
    mem = 0;
#else
    sysinfo (&usg2);
    mem = ((double)usg2.freeram) / megabyte;
//    time = usg2.ru_utime.tv_sec;
//    fprintf(stdout," time --> %ld\n ",time);
#endif
    return mem;
}

#ifdef __cplusplus
}
#endif
