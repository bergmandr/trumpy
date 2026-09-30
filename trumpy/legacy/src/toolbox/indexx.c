#include "toolbox.h"

void indexx(int n, double *arr, int *indx){
  // this function orders the array arr by generating the index array
  // indx. The result is such that the array arr[indx[i]] is the
  // ordered.        

  int l,j,ir,indxt,i,k;
  double q;             
  if (n<=0) return;          
  indx-=1;            
  arr -=1; // to account for stupid Fortran conventions
  for (j=1;j<=n;j++) indx[j]=j;
  if (n==1) {indx[1]=0; return;}
  l=(n >> 1) + 1;
  ir = n;       
  for (;;) {
    if (l>1) q=arr[(indxt=indx[--l])];
    else {
      q=arr[(indxt=indx[ir])];
      indx[ir]=indx[1];
      if (--ir==1) {
        indx[1]=indxt;
        for (k=1;k<=n;k++) indx[k]--;
        return;
      }      
    }
    i=l;
    j=l << 1;
    while (j<=ir) {
      if (j < ir && arr[indx[j]] < arr[indx[j+1]]) j++;
      if (q < arr[indx[j]]) {
        indx[i]=indx[j];
        j += (i=j);
      } else j=ir+1;
    }
    indx[i]=indxt;
  }
}

