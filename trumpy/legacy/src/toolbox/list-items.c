#include "toolbox.h"

int addListItem (int val, int *n, int max, int a[]) {
  int out, new = 1;

  for (out=0; out<*n; out++) {
    if (a[out] == val) {
      new = 0;
      break;
    }
  }
  if (new) {
    out = *n;
    a[out] = val;
    (*n)++;
    if (*n > max)
      return -1;
  }

  return out;
}
