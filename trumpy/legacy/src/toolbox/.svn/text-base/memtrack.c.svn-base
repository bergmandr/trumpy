#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TRACK_MEMORY
#  define TRACK_MEMORY_DISABLED
#else
#  undef TRACK_MEMORY
#endif

#include "memtrack.h"
#include "control.h"

typedef struct {
  void *adrs;
  size_t size;
  int indx;
  int lnum;
  char path[256];
  char func[256];
} Tracker;

#define MAX_TRACKERS 0x8000  /* = 32768 */
static int ntrackers = 0;
static Tracker tracker[MAX_TRACKERS];
static int nmalloc = 0;
static int nrealloc = 0;
static int nnullrealloc = 0;
static int nfree = 0;

/* static int _getline(char *path, int lnum, char *line) { */
/*   int i, k; */
/*   char c, fullpath[256]; */
/*   sprintf(fullpath, "%s/%s", TRUMP_HOME, path); */

/*   FILE *src = fopen(fullpath, "r"); */

/*   if ( src == NULL ) { */
/*     perr("unable to open '%s' for reading.\n", fullpath); */
/*     return 1; */
/*   } */

/*   i = 1; */
/*   while ( (c=(char)fgetc(src)) != EOF ) { */
/*     while ( c != '\n' ) */
/*       c = (char)fgetc(src); */

/*     i++; */
/*     if ( i == lnum ) { */
/*       k = 0; */
/*       while ( (line[k++]=(char)fgetc(src)) != '\n' ) */
/* 	continue; */
/*       line[k-1] = '\0'; */
/*       return 0; */
/*     } */
/*   } */

/*   perr("could not find line %d in %s\n", lnum, fullpath); */
/*   fclose(src); */

/*   return 2; */
/* } */

void *ssmalloc(size_t size, const char *func, char *path, int lnum) {
  void *adrs = malloc(size);

  /*
  pout("adding tracker %9d (%9d bytes @ 0x%08X)\n", ntrackers, size, 
       (unsigned int)adrs);
  */

  /*  this part is broken.  it needs some way to test if the LHS variable 
   *    already has memory assigned to it.
  int i;
  char line[256];
  for ( i=0; i<ntrackers; i++ ) {
    if ( strcmp(tracker[i].path, path) == 0 && tracker[i].lnum == lnum ) {
      _getline(path, lnum, line);
      perr("WARNING! pointer redefined before 'free()'\n");
      perr("%s: %4d: %s\n", path, lnum, line);
      break;
    }
  }
  */

  nmalloc++;
  tracker[ntrackers].indx = ntrackers;
  tracker[ntrackers].adrs = adrs;
  tracker[ntrackers].size = size;
  tracker[ntrackers].lnum = lnum;
  strcpy(tracker[ntrackers].path, path);
  strcpy(tracker[ntrackers].func, func);
  ntrackers++;

  if ( ntrackers == MAX_TRACKERS ) {
    perr("too many trackers have been initialized (%d)!\n", ntrackers);
    abort();
  }

  return adrs;
}

void *ssrealloc(void *adrs, size_t newsize, const char *func, char *path, int lnum) {
  if ( adrs != NULL ) {
    int i;
    void *newadrs = realloc(adrs, newsize);
    for ( i=0; i<ntrackers; i++ ) {
      if ( adrs == tracker[i].adrs ) {
	/*
	pout("reallocating tracker %9d (%9d bytes @ 0x%08X)\n", 
	     tracker[i].indx, newsize, (unsigned int)newadrs);
	*/

	nrealloc++;
	tracker[i].adrs = newadrs;
	tracker[i].size = newsize;
	tracker[i].lnum = lnum;
	strcpy(tracker[i].path, path);
	strcpy(tracker[i].func, func);
	break;
      }
    }

    return newadrs;
  }
  else {
    nnullrealloc++;
    return ssmalloc(newsize, func, path, lnum);
  }
}

void ssfree(void *adrs) {
  int i;

  for ( i=0; i<ntrackers; i++ ) {
    if ( adrs == tracker[i].adrs ) {
      /*
      pout("removing tracker %9d  (%9d bytes @ 0x%08X)\n", tracker[i].indx, 
	   tracker[i].size, (unsigned int)tracker[i].adrs);
      */
      free(adrs);
      nfree++;
      ntrackers--;
      for ( ; i<ntrackers; i++ )
	tracker[i] = tracker[i+1];
    }
  }
}

int dumpLiveTrackers() {
#ifdef TRACK_MEMORY_DISABLED
  fprintf(stdout, "Memory tracking is disabled in this build.\n");
#else
  int i;

  fprintf(stdout, "There have been %d calls to memory allocation functions.\n", 
	  nmalloc+nrealloc+nfree);
  fprintf(stdout, "    calls to 'malloc()':  %d\n", nmalloc);
  fprintf(stdout, "    calls to 'realloc()': %d (%d initially 'NULL')\n", 
	  nrealloc, nnullrealloc);
  fprintf(stdout, "    calls to 'free()':    %d\n", nfree);
  fprintf(stdout, "    balance: %d\n\n", nmalloc-nfree);
  fprintf(stdout, "There are %d live trackers.\n", ntrackers);
  for ( i=0; i<ntrackers; i++ ) {
    fprintf(stdout, "  %9d: ID %d: %s line %d : %9d bytes @ 0x%08X\n", i+1, 
	    tracker[i].indx, tracker[i].path, tracker[i].lnum, 
	    (int)tracker[i].size, (unsigned int)tracker[i].adrs);
  }
  fprintf(stdout, "\n");
#endif

  return ntrackers;
}

