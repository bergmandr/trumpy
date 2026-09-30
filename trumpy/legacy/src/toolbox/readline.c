#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "toolbox.h"

#define MAX_STRLEN 256

/*
 *  A utility for reading text files and parsing them in "bash" or "config" 
 * style.  It reads one line of text (string of characters terminated by a 
 * neline ('\n') character), and converts the words in the line to an argument 
 * list array.  Words are defined as sequences of non-whitespace characters 
 * in the line.
 */
int readline(FILE *fp, char *line, char args[][MAX_STRLEN]) {
  int i, k, nc=0, narg=0;
  char c;

  while ( (c=(char)fgetc(fp)) != '\n' ) {
    if ( c == EOF )
      return -1;
    else if ( c == '#' ) {
      /* scan to the end of the line */
      while ( fgetc(fp) != '\n' ) { }
      break;
    }

    line[nc++] = c;
  }
  line[nc] = '\0';

  /*
   *  extract the distinct words in the line
   */
  for ( i=0; i<nc; ) {
    if ( isspace(line[i]) ) {
      /* scan to next non-whitespace character */
      while ( isspace(line[i]) && i < nc )
	i++;
    }
    else {
      /* copy word to args[narg] */
      k = 0;
      while ( ! isspace(line[i]) && i < nc )
	args[narg][k++] = line[i++];
      args[narg][k] = '\0';
      narg++;

      //    printf("  %2d  '%s'\n", narg, args[narg-1]);
    }
  }

  return narg;
}

