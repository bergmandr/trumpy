#include <stdio.h>
#include "random.h"

static Ranq2 ranq2;

/* This is really a constructor */
#ifdef _OPENMP
void ranq2setup(int i, unsigned long long j) {
    ranq2.v[i] = 4101842887655102017LL;
    ranq2.w[i] = 1;
    ranq2.v[i] ^= j;
    ranq2.w[i] = ranq2int64(i);
    ranq2.v[i] = ranq2int64(i);
}
#else
void ranq2setup(unsigned long long j) {
  ranq2.v = 4101842887655102017LL;
  ranq2.w = 1;
  ranq2.v ^= j;
  ranq2.w = ranq2int64();
  ranq2.v = ranq2int64();
}
#endif



#ifdef _OPENMP
unsigned long long ranq2int64(int i) {

  ranq2.v[i] ^= ranq2.v[i] >> 17; 

  ranq2.v[i] ^= ranq2.v[i] << 31; 

  ranq2.v[i] ^= ranq2.v[i] >> 8;

  ranq2.w[i] = 4294957665U*(ranq2.w[i] & 0xffffffff) + (ranq2.w[i] >> 32);

  return ranq2.v[i] ^ ranq2.w[i];
}
#else
unsigned long long ranq2int64() {

  ranq2.v ^= ranq2.v >> 17; 

  ranq2.v ^= ranq2.v << 31; 

  ranq2.v ^= ranq2.v >> 8;

  ranq2.w = 4294957665U*(ranq2.w & 0xffffffff) + (ranq2.w >> 32);

  return ranq2.v ^ ranq2.w;
}
#endif


#ifdef _OPENMP
double ranq2doub(int i) {return 5.42101086242752217E-20 * ranq2int64(i); }
#else
double ranq2doub() {return 5.42101086242752217E-20 * ranq2int64(); }
#endif

#ifdef _OPENMP
unsigned long ranq2int32(int i) { return (unsigned long) ranq2int64(i);}
#else
unsigned long ranq2int32() { return (unsigned long) ranq2int64();}
#endif

#ifdef _OPENMP
void ranq2printState(int i) {printf("v[%d] = %llu, w[%d] = %llu\n",i,ranq2.v[i],i,ranq2.w[i]);}
#else
void ranq2printState() {printf("v = %llu, w = %llu\n",ranq2.v,ranq2.w);}
#endif

#ifdef _OPENMP
void ranq2setState(int i, unsigned long long v, unsigned long long w) {
  ranq2.v[i] = v; ranq2.w[i] = w;}
#else
void ranq2setState(unsigned long long v, unsigned long long w) {
  ranq2.v = v; ranq2.w = w;}
#endif

#ifdef _OPENMP
void ranq2getState(int i, unsigned long long *v, unsigned long long *w) {
  *v = ranq2.v[i]; *w = ranq2.w[i];}
#else
void ranq2getState(unsigned long long *v, unsigned long long *w) {
  *v = ranq2.v; *w = ranq2.w;}
#endif
