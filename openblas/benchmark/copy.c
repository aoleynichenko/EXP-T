/***************************************************************************
Copyright (c) 2014, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#ifdef __CYGWIN32__
#include <sys/time.h>
#endif
#include "common.h"


#undef COPY

#ifdef COMPLEX
#ifdef DOUBLE
#define COPY   BLASFUNC(zcopy)
#else
#define COPY   BLASFUNC(ccopy)
#endif
#else
#ifdef DOUBLE
#define COPY   BLASFUNC(dcopy)
#else
#define COPY   BLASFUNC(scopy)
#endif
#endif

#if defined(__WIN32__) || defined(__WIN64__)

#ifndef DELTA_EPOCH_IN_MICROSECS
#define DELTA_EPOCH_IN_MICROSECS 11644473600000000ULL
#endif

int gettimeofday(struct timeval *tv, void *tz){

  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag;

  if (NULL != tv)
    {
      GetSystemTimeAsFileTime(&ft);

      tmpres |= ft.dwHighDateTime;
      tmpres <<= 32;
      tmpres |= ft.dwLowDateTime;

      /*converting file time to unix epoch*/
      tmpres /= 10;  /*convert into microseconds*/
      tmpres -= DELTA_EPOCH_IN_MICROSECS;
      tv->tv_sec = (long)(tmpres / 1000000UL);
      tv->tv_usec = (long)(tmpres % 1000000UL);
    }

  return 0;
}

#endif

#if !defined(__WIN32__) && !defined(__WIN64__) && !defined(__CYGWIN32__) && 0

static void *huge_malloc(BLASLONG size){
  int shmid;
  void *address;

#ifndef SHM_HUGETLB
#define SHM_HUGETLB 04000
#endif

  if ((shmid =shmget(IPC_PRIVATE,
		     (size + HUGE_PAGESIZE) & ~(HUGE_PAGESIZE - 1),
		     SHM_HUGETLB | IPC_CREAT |0600)) < 0) {
    printf( "Memory allocation failed(shmget).\n");
    exit(1);
  }

  address = shmat(shmid, NULL, SHM_RND);

  if ((BLASLONG)address == -1){
    printf( "Memory allocation failed(shmat).\n");
    exit(1);
  }

  shmctl(shmid, IPC_RMID, 0);

  return address;
}

#define malloc huge_malloc

#endif

int main(int argc, char *argv[]){

  FLOAT *x, *y;
  FLOAT alpha[2] = { 2.0, 2.0 };
  blasint m, i;
  blasint inc_x=1,inc_y=1;
  int loops = 1;
  int l;
  char *p;

  int from =   1;
  int to   = 200;
  int step =   1;

  struct timeval start, stop;
  double time1 = 0.0, timeg = 0.0;
  long nanos = 0;
  time_t seconds = 0;
  struct timespec time_start = { 0, 0 }, time_end = { 0, 0 };

  argc--;argv++;

  if (argc > 0) { from     = atol(*argv);		argc--; argv++;}
  if (argc > 0) { to       = MAX(atol(*argv), from);	argc--; argv++;}
  if (argc > 0) { step     = atol(*argv);		argc--; argv++;}

  if ((p = getenv("OPENBLAS_LOOPS")))  loops = atoi(p);
  if ((p = getenv("OPENBLAS_INCX")))   inc_x = atoi(p);
  if ((p = getenv("OPENBLAS_INCY")))   inc_y = atoi(p);

  fprintf(stderr, "From : %3d  To : %3d Step = %3d Inc_x = %d Inc_y = %d Loops = %d\n", from, to, step,inc_x,inc_y,loops);

  if (( x = (FLOAT *)malloc(sizeof(FLOAT) * to * abs(inc_x) * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

  if (( y = (FLOAT *)malloc(sizeof(FLOAT) * to * abs(inc_y) * COMPSIZE)) == NULL){
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }

#ifdef __linux
  srandom(getpid());
#endif

  fprintf(stderr, "   SIZE       Flops\n");

  for(m = from; m <= to; m += step)
  {

   timeg=0;

   fprintf(stderr, " %6d : ", (int)m);
   for(i = 0; i < m * COMPSIZE * abs(inc_x); i++){
       x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
   }

   for(i = 0; i < m * COMPSIZE * abs(inc_y); i++){
       y[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
   }

   for (l=0; l<loops; l++)
   {
       clock_gettime(CLOCK_REALTIME, &time_start);
       COPY (&m, x, &inc_x, y, &inc_y );
       clock_gettime(CLOCK_REALTIME, &time_end);

       nanos = time_end.tv_nsec - time_start.tv_nsec;
       seconds = time_end.tv_sec - time_start.tv_sec;

       time1 = seconds + nanos / 1.e9;
       timeg += time1;
   }

      timeg /= loops;

      fprintf(stderr,
	    " %10.2f MBytes %12.9f sec\n",
	    COMPSIZE * sizeof(FLOAT) * 1. * (double)m / timeg / 1.e6, timeg);

  }

  return 0;
}

// void main(int argc, char *argv[]) __attribute__((weak, alias("MAIN__")));
