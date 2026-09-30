#ifndef _RANDOM_H_
#define _RANDOM_H_


#ifdef _OPENMP
#include <omp.h>
#define MAXNTHREAD 24
#define RANDOM_NUMBER ranq2doub(omp_get_thread_num())

#else


//#define RANDOM_NUMBER getRandomNumber()
#define RANDOM_NUMBER ranq2doub()
#endif

/* in "random.c" */
#ifdef __cplusplus
extern "C"{
#endif
void getSeeds(int *hi, int *lo);
void setSeeds(int hi, int lo);
double getRandomNumber();
int prand(double mu);
double grand(void);
double trand(void);

/* in "ran2.c" */
float ran2(long *idum);

/* from "ranss.c" */
double ranss(int *h, int *l);
#ifdef __cplusplus
}
#endif

/* Implentation of Ranq2 from NR3: ranq2.c */
#ifdef _OPENMP
typedef struct {
  unsigned long long v[MAXNTHREAD];
  unsigned long long w[MAXNTHREAD];
} Ranq2;
#else
typedef struct {
  unsigned long long v;
  unsigned long long w;
} Ranq2;
#endif

#ifdef __cplusplus
extern "C"{
#endif

#ifdef _OPENMP
void ranq2setup(int i, unsigned long long j);
unsigned long long ranq2int64(int i);
double ranq2doub(int i);
unsigned long ranq2int32(int i);
void ranq2printState(int i);
void ranq2setState(int i, unsigned long long v, unsigned long long w);
void ranq2getState(int i, unsigned long long *v, unsigned long long *w);
#else
void ranq2setup(unsigned long long j);
unsigned long long ranq2int64();
double ranq2doub();
unsigned long ranq2int32();
void ranq2printState();
void ranq2setState(unsigned long long v, unsigned long long w);
void ranq2getState(unsigned long long *v, unsigned long long *w);
#endif

#ifdef __cplusplus
}
#endif

#endif
