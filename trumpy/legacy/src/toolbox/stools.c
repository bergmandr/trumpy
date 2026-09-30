/*
 *  Tools for manipulating strings that may or may not already exist.
 */

#include <stdlib.h>
#include <string.h>
#include "stools.h"

/*
 *  returns '0' if 'str' does not end with 'sfx', '1' if it does, and '-1' if 
 *    'sfx' is longer than 'str'.
 */
int endswith(char *str, const char *sfx) {
/*   int i; */
/*   int n = strlen(str); */
/*   int k = strlen(sfx); */

/*   if ( k > n ) */
/*     return -1; */

/*   for ( i=n-1; i>=k; i-- ) */
/*     if ( str[i] != sfx[i] ) */
/*       return 0; */

/*   return 1; */
  if ( strstr(str, sfx) != NULL )
    return 1;
  else
    return 0;
}

/*
 *  returns '0' if 'str' does not start with 'pfx', '1' if it does, and '-1' 
 *    if 'pfx' is longer than 'str'.
 */
int startswith(char *str, const char *pfx) {
  int i;
  int n = strlen(str);
  int k = strlen(pfx);

  if ( k > n )
    return -1;

  for ( i=0; i<k; i++ )
    if ( str[i] != pfx[i] )
      return 0;

  return 1;
}
