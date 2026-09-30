#ifndef _MEMTRACK_H_
#define _MEMTRACK_H_

#ifdef TRACK_MEMORY
#  define realloc(a,s) ssrealloc(a, s, __FUNCTION__, __FILE__, __LINE__)
#  define malloc(s)    ssmalloc(s, __FUNCTION__, __FILE__, __LINE__)
#  define free(a)      ssfree(a)
#endif

void *ssrealloc(void *adrs, size_t newsize, const char *func, char *path, int lnum);
void *ssmalloc(size_t size, const char *func, char *path, int lnum);
void ssfree(void *adrs);
int dumpLiveTrackers(void);

#endif
