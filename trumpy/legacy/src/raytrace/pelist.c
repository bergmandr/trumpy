#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "control.h"
#include "event.h"
#include "fdconstants.h"
#include "toolbox.h"

#include "airshower.h"
#include "tacalibration.h"
#include "taelectronics.h"
#include "raytrace.h"
#include "showerlib.h"
#include "track.h"
#include "fdsite.h"
#include "acpttrack.h"
#include "pelist.h"

int makePEList(const AcceptedTrack* at, PEList* pelist) {
  int i, j, k, m;
  int na, nr;

  na = 0;
  for (i=0; i<GEOFD_MAXMIR; i++)
    for (j=0; j<GEOFD_MIRTUBE; j++)
      na += at->len[i][j];

  if (na == 0) {
    pelist->peRecord = NULL;
    pelist->nper = 0;
    return 0;
  }

  pelist->peRecord = (PERecord**)malloc(na*sizeof(PERecord*));
  nr = 0;

  for (i=0; i<GEOFD_MAXMIR; i++) {
    for (m=0; m<at->nmir; m++) {
      if (at->mir[m] == i) {
	for (j=0; j<GEOFD_MIRTUBE; j++) {
	  for (k=0; k<at->len[m][j]; k++) {
	    pelist->peRecord[nr] = (PERecord*)malloc(sizeof(PERecord));
	    pelist->peRecord[nr]->cam  = at->mir[m];
	    pelist->peRecord[nr]->tube = j;
	    pelist->peRecord[nr]->npe  = at->npe[m][j][k];
	    pelist->peRecord[nr]->tpe  = (double)at->tpe[m][j][k];
	    nr++;
	  }
	}
	break;
      }
    }
  }
  pelist->nper = nr;

  return pelist->nper;
}

void clearPEList(PEList* pe) {
  int i;
  for (i=0; i<pe->nper; i++)
    free(pe->peRecord[i]);
  free(pe->peRecord);
}

int cmpPERecord(const void *vpe1, const void *vpe2) {
  PERecord* pe1 = *(PERecord**)vpe1;
  PERecord* pe2 = *(PERecord**)vpe2;
  if (pe1->cam != pe2->cam)
    return (pe1->cam < pe2->cam)? -1 : +1;
  if (pe1->tube != pe2->tube)
    return (pe1->tube < pe2->tube)? -1 : +1;
  if (pe1->tpe != pe2->tpe)
    return (pe1->tpe < pe2->tpe)? -1 : +1;
  else
    return 0;
}

void sortPEList(PEList* pelist) {
  qsort(pelist->peRecord,pelist->nper,sizeof(PERecord*),cmpPERecord);
}

int fillPETimes(const PEList* ls, PETimes* tm) {
  int iper = 0;
  int nper = ls->nper;
  int iper2;
  int cam,tube,len;

  int newcam, icam;

  if (nper > 0) {
    tm->t1 = ls->peRecord[iper]->tpe;
    tm->t2 = ls->peRecord[iper]->tpe;
  }
  else {
    tm->t1 = 0.;
    tm->t2 = 0.;
  }

  tm->nmir = 0;

  while (iper < nper) {
    /*
     * Find block of PE's from a given camera/tube
     * Make sure not to go past end of array
     */
    iper2 = iper + 1;
    while (iper2 < nper &&
	   ls->peRecord[iper]->cam == ls->peRecord[iper2]->cam &&
	   ls->peRecord[iper]->tube == ls->peRecord[iper2]->tube)
      iper2++;

    cam = ls->peRecord[iper]->cam;

    newcam = TRUE;
    for (icam=0; icam<tm->nmir; icam++) {
      if (tm->mir[icam] == cam) {
	newcam = FALSE;
	break;
      }
    }
    if (newcam) {
      icam = tm->nmir;
      tm->mir[icam] = cam;
      tm->nmir++;
    }
    
    /* Alocate space in t */
    len = iper2 - iper;
    tube = ls->peRecord[iper]->tube;
    tm->n[icam][tube] = (int*)malloc(len*sizeof(int));
    tm->t[icam][tube] = (int*)malloc(len*sizeof(int));

    /* Fill new arrays in t */
    int ip,j;
    double lastTPE = -1.;
    for (j=-1,ip=iper; ip<iper2; ip++) {
      if ( fabs(round(ls->peRecord[ip]->tpe) - lastTPE) < 1.0e-10 && j >= 0) {
	/*
	 * Combine PERecords for a given camera/tube which round to the 
	 * same ns.  In this case the time will be an integer number of ns
	 */
	tm->n[icam][tube][j] += ls->peRecord[ip]->npe;
	tm->t[icam][tube][j] = lastTPE;
      }
      else {
	j++;
	tm->n[icam][tube][j] = ls->peRecord[ip]->npe;
	tm->t[icam][tube][j] = ls->peRecord[ip]->tpe;
	tm->t1 = min(tm->t1,ls->peRecord[ip]->tpe);
	tm->t2 = max(tm->t2,ls->peRecord[ip]->tpe);
	lastTPE = round(ls->peRecord[ip]->tpe);
      }
    }

    len = j+1;
    tm->len[icam][tube] = len;    
    tm->n[icam][tube] = (int*)realloc(tm->n[icam][tube],len*sizeof(int));
    tm->t[icam][tube] = (int*)realloc(tm->t[icam][tube],len*sizeof(int));

    /* Move to next set in ls->peRecord */
    iper = iper2;
  }

  return 0;
}

void zeroPETimes(PETimes* pe) {
  int j,k;
  for (j=0;j<GEOFD_MAXMIR;j++)
    for (k=0;k<GEOFD_MIRTUBE;k++) {
      pe->n[j][k] = NULL;
      pe->t[j][k] = NULL;
      pe->len[j][k] = 0;
    }
  pe->t1 = 0.;
  pe->t2 = 0.;
}

void clearPETimes(PETimes* pe) {
  int j,k;
  for (j=0;j<GEOFD_MAXMIR;j++)
    for (k=0;k<GEOFD_MIRTUBE;k++) {
      /*
       * pe->n and pe-> are allocated in fillPETimes from a PEList.
       * They are only allocated when a given camera/tube appear, so
       * need to check before free'ing
       */
      if (pe->n[j][k] != NULL) free(pe->n[j][k]);
      if (pe->t[j][k] != NULL) free(pe->t[j][k]);
      pe->n[j][k] = NULL;
      pe->t[j][k] = NULL;
      pe->len[j][k] = 0;
    }
  pe->t1 = 0.;
  pe->t2 = 0.;
}
