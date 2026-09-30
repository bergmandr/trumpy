#ifndef _PELIST_H_
#define _PELIST_H_

typedef struct {
  int cam;    // Camera in which PE's hit
  int tube;   // Tube in which PE's hit
  int npe;    // Number of PE's which hit
  double tpe; // Time (ns) at which PE's hit
} PERecord;

typedef struct {
  int nper;
  PERecord* *peRecord;
} PEList;

/* in "pelist.c" */
//int makePEList(const AcceptedTrack* at, PEList* pelist);
void clearPEList(PEList* pelist);
int cmpPERecord(const void *pe1, const void *pe2);
void sortPEList(PEList* pelist);
int fillPETimes(const PEList* pelist, PETimes* petimes);
void clearPETimes(PETimes* pe);
void zeroPETimes(PETimes* pe);

#endif
