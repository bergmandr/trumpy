#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tlelectronics.h"
// #include "maxmin.h"
// #include "dst_std_types.h"

// #include "parameters.h"
// #include "mir_hit.h"
// #include "mirrs3.h"
// #include "glb_prm.h"

integer4 iang_m_(integer4 *f77i, integer4 *f77j) {

  integer4 return_value;
  
  integer4 iy1,ix1,iy2,ix2;
  integer4 iy,ix;
  
  integer4 i,j;
  
  i = *f77i - 1;
  j = *f77j - 1;
  
  if (i == j) {
    return_value = 0;
    return return_value;
  }
  
  mirror_xy_(&(mir_hit_.id_hit_mir[i]),&iy1,&ix1);
  mirror_xy_(&(mir_hit_.id_hit_mir[j]),&iy2,&ix2);
  iy = abs(iy1 - iy2);
  ix = abs(ix1 - ix2);
  
  if (ix >= 3 || iy >= 3)
    return_value = 3;
  else if (iy == 0)
    return_value = abs(ix1 - ix2);
  else if (iy == 1) {
    return_value = ix;
    if (ix2 >= ix1)
      return_value += (iy2 % 2);
    if (ix2 <= ix1)
      return_value += (iy1 % 2);
  }
  else if (iy == 2)
    return_value = 2 + ix/2;
  else
    return_value = 3;
  
  if (return_value > 3)
    return_value = 3;
  
  return return_value;
}
    