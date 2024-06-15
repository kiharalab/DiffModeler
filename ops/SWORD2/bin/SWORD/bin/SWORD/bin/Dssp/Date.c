#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif
#include <stdio.h>
static  char *mon[]={"JAN","FEB","MAR","APR","MAY","JUN",
                "JUL","AUG","SEP","OCT","NOV","DEC"};
/* PROCEDURE DATE(VAR DATESTRING:PACKED ARRAY[1..11] OF CHAR);EXTERN; */
/* activate DATE by removing comment brackets if necessary */
/***/

void Date(string)
char* string;
{
  time_t tt;
  struct tm *t;
  time(&tt);
  t=localtime(&tt);
  sprintf(string,"%d-%s-%4d",t->tm_mday,mon[t->tm_mon],t->tm_year+1900);
}
