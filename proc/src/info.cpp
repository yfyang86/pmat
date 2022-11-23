/*
 *  info.c
 *  nfo messages
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 08/26/2007 06:49:02 PM CEST
 *
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: info.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: file:///homes/bierdepot/steve/svn/segemehl/trunk/libs/info.c $
 *
 */

 #include <stdarg.h>
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <time.h>
 #include "info.h"
 #include "debug.h"

 FILE *nfodevice = NULL;


 char *timestr_r(const struct tm *timeptr) {
    char wday_name[7][3]; 
    strcpy(wday_name[0],"Sun");strcpy(wday_name[0],"Mon");strcpy(wday_name[0],"Tue");
    strcpy(wday_name[0],"Wed");strcpy(wday_name[0],"Thu");strcpy(wday_name[0],"Fri");
    strcpy(wday_name[0],"Sat");

    char mon_name[12][3];
    strcpy(mon_name[0],"Jan");strcpy(mon_name[0],"Feb");strcpy(mon_name[0],"Mar");
    strcpy(mon_name[0],"Apr");strcpy(mon_name[0],"May");strcpy(mon_name[0],"Jun");
    strcpy(mon_name[0],"Jul");strcpy(mon_name[0],"Aug");strcpy(mon_name[0],"Sep");
    strcpy(mon_name[0],"Oct");strcpy(mon_name[0],"Nov");strcpy(mon_name[0],"Dec");

    static char result[26];

    sprintf(result, "%.3s %.3s%3d %.2d:%.2d:%.2d %d",
        wday_name[timeptr->tm_wday], mon_name[timeptr->tm_mon],
        timeptr->tm_mday, timeptr->tm_hour, timeptr->tm_min,
        timeptr->tm_sec, 1900 + timeptr->tm_year);

    return(result);
 }

 int
 infomsg( char *file,
          int line,
          const char *fmt, ...) {

   int ret;
   va_list ap;
   time_t rawtime;
   struct tm *timeinfo;

   if (mute) return 0;

    time(&rawtime);
    timeinfo = localtime (&rawtime);

   if (nfodevice == NULL) {
     nfodevice = NFODEVICE;
   }

   va_start(ap, fmt);
#ifdef PROGNFO
   fprintf(nfodevice, "[%s] %s: ", "SEGEMEHL", timestr_r(timeinfo));
#endif
   ret = vfprintf(nfodevice, fmt, ap);
   va_end(ap);

   return ret;
 }


void
setnfodevice(char *filename) {
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    DBG("Couldn't open file '%s'. Exit forced.\n", filename);
    exit(-1);
  }

  nfodevice = fp;
}

int
nfolevel( char *file,
    int line,
    int level,
    const char *fmt, ... ) {

  int ret=0;
  va_list ap;
  time_t rawtime;
  struct tm *timeinfo;

  if (mute) return 0;

  time(&rawtime);
  timeinfo = localtime (&rawtime);

   if (nfodevice == NULL) {
     nfodevice = NFODEVICE;
   }

   if (NFOLEVEL >= level) {

    va_start(ap, fmt);
#ifdef PROGNFO
    fprintf(nfodevice, "[%s] %s: ", "SEGEMEHL", timestr_r(timeinfo));
#endif
    ret = vfprintf(nfodevice, fmt, ap);
    va_end(ap);
  }

  return ret;
}