#ifndef _SSHBOOK_H_
#define _SSHBOOK_H_

void openHistogramFile(char *name);

void closeHistogramFile(void);

void hbook1(int id, char *chtitl, int ncx, float xlow, float xup);

void hbook2(int id, char *chtitl, int ncx, float xlow, float xup, int ncy, float ylow, float yup);

void hbprof(int id, char *chtitl, int ncx, float xlow, float xup, float ymin, float ymax);

void hfill(int id, float x, float y, float w);

void hf1(int id, float x, float w);

void hf2(int id, float x, float y, float w);

void hfn(int id, float *data);

void hbkrow(int id, char *title, int nentries);

void rowtag(int id, char *name);

int hbk_init_(char *path);

int hbk_deftag_(void);

int hbk_nt_tag_(int *icol,char *name);

int hbk_row_(int *id, char *name, int *ncol);

int hbk_1d_(int *id, char *name, int *nx, float *xmi, float *xma);

int hbk_2d_(int *id, char *name,int *nx, float *xmi, float *xma, int *ny, float *ymi, float *yma);

int hbk_prof_(int *id, char *name, int *nx, float *xmi, float *xma, float *ymi, float *yma);

int hbk_end_(void);

void hf1_(int *id, float *x, float *w);

void hf2_(int *id, float *x, float *y, float *w);

void hfn_(int *id, float *x);

void hfill_(int *id, float *x, float *y, float *w);

#endif
