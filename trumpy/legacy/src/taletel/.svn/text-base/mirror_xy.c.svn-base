#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "dst_std_types.h"


void mirror_xy_(integer4 *ipmt, integer4 *irow, integer4 *icol) {
  integer4 id;
  if (*ipmt < 1 || *ipmt > 256) {
    printf(" mirror_xy called for bad ipmt = %d\n",*ipmt);
    abort();
    *irow = 1;
    *icol = 1;
    return;
  }
  id = *ipmt - 1;
  *icol = (id % 16) + 1;
  *irow = (id/16) + 1;
  return;
}