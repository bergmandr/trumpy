#ifndef _TLELECTRONICS_H_
#define _TLELECTRONICS_H_


#include "dst_std_types.h"







#include "trump.h"

// typedef struct {
//   double xmu_noise;
// //   double tau_
//   double ped_fadc[GEOTL_MAXMIR][nchan_mir];
//   double gain_fadc[GEOTL_MAXMIR][nchan_mir];
// 
// } TLElectronicsState;

// extern TLElectronicsState tes;

// // #include "_control.h"
// #ifndef __CONTROL_H_
// #define __CONTROL_H_
// typedef struct {
//   integer4 ifluc;
//   integer4 inoise_mir;
//   integer4 inoise_eye;
//   integer4 nev_plot;
//   integer4 ifilter_r;
//   integer4 ifilter_r_m;
//   integer4 icerenkov;
//   integer4 ireconstruct;
//   integer4 isim_fadc;
//   integer4 isim_utah;
//   integer4 inew;
// } control;
// 
// extern control control_;
// #endif



// #include "db_mcru.h"
// #ifndef _DB_MCRU_H_
// #define _DB_MCRU_H_
// typedef struct {
//   integer4 db_date[3];
//   integer4 db_events;
//   real4 runtime_night;
//   real4 db_dac1;
//   real4 db_dac2;
//   real4 db_thresh1;
//   real4 db_thresh2;
//   integer4 db_trig_version;
//   integer4 db_adj_opt;
//   integer4 db_prescale;
//   real4 db_avnoise;
//   real4 db_rmsnoise;
//   integer4 db_deadmir;
//   integer4 db_dead[42];
// } db_info;

// extern db_info db_info_;

// typedef struct {
//   integer4 use_data_set;
//   integer4 use_database;
// } param_src;
// 
// extern param_src param_src_;

// typedef struct {
//   integer4 setnumber;
// } which_set;
// 
// extern which_set which_set_;

// typedef struct {
//   integer4 def_sec_trig;
//   integer4 def_adj_opt;
// } trig_def;
// 
// extern trig_def trig_def_;

// integer4 adjinit[14][6] = {
//   { 1, 2, 3, 4, 0, 0},
//   { 2, 1, 3, 4, 0, 0},
//   { 3, 4, 1, 2, 5, 6},
//   { 4, 3, 1, 2, 5, 6},
//   { 5, 6, 3, 4, 8, 7},
//   { 6, 5, 3, 4, 8, 7},
//   { 7, 8, 5, 6, 9,10},
//   { 8, 7, 5, 6, 9,10},
//   { 9,10, 8, 7,11,12},
//   {10, 9, 8, 7,11,12},
//   {11,12, 9,10,13,14},
//   {12,11, 9,10,13,14},
//   {13,14,11,12, 0, 0},
//   {14,13,11,12, 0, 0}
// };

// #endif



// #include "fadcgain.h"
// #ifndef _FADC_GAIN_H
// #define _FADC_GAIN_H

// typedef struct
// {
//   real4 ped_fadc[GEOTL_MAXMIR][nchan_mir];
// } pedfadc;
// 
// extern pedfadc pedfadc_;

// typedef struct
// {
//   real4 gain_fadc[GEOTL_MAXMIR][nchan_mir];
// } gainfadc;
// 
// extern gainfadc gainfadc_;

// typedef struct
// {
//   real4 tau_fadc;
//   real4 tau_16;
//   real4 tau_neg;
// } taufadc;
// 
// extern taufadc taufadc_;
// #endif



// #include "filter_m.h"
#ifndef _FILTER_M_H_
#define _FILTER_M_H_
typedef struct
{
  integer4 lcut_m[nhit_mir_max];
  real4 dth_m[nhit_mir_max];
  integer4 iedge_m[nhit_mir_max];
} filter_m;

extern filter_m filter_m_;
#endif




// #include "glb_prm.h"
// #ifndef _GLB_PRM_H_
// #define _GLB_PRM_H_

// typedef struct
// {
//   real4 xmu_noise;
//   real4 dt_adc;
//   real4 sigt_adc;
//   real4 thresh_1;
//   real4 thresh_2;
//   real4 side;
// } glb_prm;

// extern glb_prm glb_prm_;
// #endif




// #include "mir_hit.h"
#ifndef _MIR_HIT_H_
#define _MIR_HIT_H_

typedef struct
{
  integer4 nhit_mir;
  integer4 id_hit_mir[nhit_mir_max];
  integer4 it0_hit_mir[nhit_mir_max];
  integer4 nt_hit_mir[nhit_mir_max];
  real4 t_hit_mir[nhit_mir_max];
  real4 dt_hit_mir[nhit_mir_max];
  integer4 npe_hit_mir[nhit_mir_max];
  integer4 ichan_trig_mir[nhit_mir_max];
  integer4 it0_trig_mir[nhit_mir_max];
  real4 sig_hit_mir[nhit_mir_max];
  integer4 ichan_hit_mir[nhit_mir_max];
  real4 sig_trig_mir[nhit_mir_max];
  integer4 max_fadc_hit_mir[nhit_mir_max];
  integer4 nfadc_hit_mir[nhit_mir_max];
} mir_hit;

extern mir_hit mir_hit_;
#endif




// #include "mir_trig.h"
#ifndef _MIR_TRIG_H_
#define _MIR_TRIG_H_

typedef struct {
  integer4 nmir_trig_1;
  integer4 nmir_trig_2;
  integer4 nmir_trig_utah;
  integer4 idmir_trig_1[GEOTL_MAXMIR];
  integer4 idmir_trig_2[GEOTL_MAXMIR];
  integer4 idmir_trig_utah[GEOTL_MAXMIR];
} mir_trig;

extern mir_trig mir_trig_;

#endif



// #include "mirrs3.h"
#ifndef _MIRRS3_H_
#define _MIRRS3_H_

typedef struct {
  integer4 mirs;
  integer4 ntbmir[GEOTL_MAXMIR];
  integer4 mrtbid[ntube_eye];
  integer4 mirtub[GEOTL_MAXMIR][256];
} mirrs;

extern mirrs mirrs_;


typedef struct {
  real4 vmir[GEOTL_MAXMIR][3];
  integer4 tbv[ntube_eye][3];
} centre;

extern centre centre_;
#endif




// #include "numcon.h"
#ifndef _NUMCON_H_
#define _NUMCON_H_
typedef struct {
  real4 pi;
  real4 radian;
  real4 root2;
  real4 root3;
  real4 twopi;
  real4 root_2pi;
  real4 halfpi;
  real4 cinv;
} numcon;

extern numcon numcon_;
#endif




// #include "petime.h"
#ifndef _PETIME_H_
#define _PETIME_H_

typedef struct {
  double t0_pmt[nhit_max];
  double tl_pmt[nhit_max];
  integer4 npe_pmt[nhit_max];
  double tav_pmt[nhit_max];
} petime;

extern petime petime_;

typedef struct {
  integer4 id_temp[nhit_max];
  double tpe_temp[nhit_max][npe_max];
  integer4 incr_npe_temp[nhit_max][npe_max];
  integer4 ntpe_temp[nhit_max];
  double tgen_av[nhit_max];
} petemp;

extern petemp petemp_;

typedef struct {
  integer4 nmir_view;
  integer4 imir_view[GEOTL_MAXMIR];
  integer4 jmir[GEOTL_MAXMIR];
  double pemir[GEOTL_MAXMIR];
  double pemir_sc[GEOTL_MAXMIR];
  double t0_mir[GEOTL_MAXMIR];
  double t1_mir[GEOTL_MAXMIR];
  integer4 nt_mir[GEOTL_MAXMIR];
} pemir;

extern pemir pemir_;

typedef struct {
  integer4 ntba;
  integer4 itba[nhit_max];
  integer4 jhit[ntube_eye];
} idtemp;

extern idtemp idtemp_;
#endif




// #include "skynoise.h"
#ifndef _SKYNOISE_H_
#define _SKYNOISE_H_

#define nnbins1 100
#define nnbins2 100
#define nnbins3 100

typedef struct {
  real4 skynoise2a[nnbins2];
  real4 skynoise2t[nnbins2];
  real4 skynoise2d[nnbins2];
  real4 skynoise3a[nnbins3];
  real4 skynoise3t[nnbins3];
  real4 skynoise3d[nnbins3];
  real4 skynoise1[nnbins1];
  real4 norm1;
  real4 norm2t;
  real4 norm2a;
  real4 norm2d;
  real4 norm3a;
  real4 norm3t;
  real4 norm3d;
} skynoise;

extern skynoise skynoise_;

typedef struct {
  real4 amb_noise[42];
} bg_noise;

extern bg_noise bg_noise_;

#endif




// #include "trigger_az.h"
#ifndef _TRIGGER_AZ_H_
#define _TRIGGER_AZ_H_

integer4 fraw1_copy_second;
integer4 fraw1_copy_clkcnt;

integer2 fraw1_copy_mir_num;
integer2 fraw1_copy_num_chan;

integer2 fraw1_copy_channel[320];
integer2 fraw1_copy_it0_chan[320];

integer2 fraw1_copy_nt_chan[320];
integer1 fraw1_copy_m_fadc[320][nt_chan_max];

integer2 ftrg1_copy_trig_code;
integer2 ftrg1_copy_mir_num;

integer2 ftrg1_copy_nhit_dsp;
integer2 ftrg1_copy_ichan_dsp[320];
integer2 ftrg1_copy_it0_dsp[320];

integer2 ftrg1_copy_nt_dsp[320];
integer2 ftrg1_copy_tav_dsp[320];

integer2 ftrg1_copy_sigt_dsp[320];
integer2 ftrg1_copy_nadc_dsp[320];

integer2 ftrg1_copy_t_pld_start;

integer2 ftrg1_copy_t_pld_end;

integer2 ftrg1_copy_row_pattern;

integer2 ftrg1_copy_col_pattern;

integer2 ftrg1_copy_wp[16];
integer2 ftrg1_copy_delay;
#endif

typedef struct {
  integer4 nmir_try;
} rawcheck;

extern rawcheck rawcheck_;

typedef struct {
  integer4 itbflag;
  integer4 itry;
} jim;

extern jim jim_;

typedef struct {
  integer4 nsingle;
  integer4 nprims;
} prescale;

extern prescale prescale_;

// typedef struct {
//   real4 calib_corr;
// } calib;
// 
// extern calib calib_;

#ifdef __cplusplus
extern "C"{
#endif

// void const_init_(void);
void init_noise_(void);

void electronics_fadc_(TLCalibration *tes);
int sec_trigger_(TLCalibration *tes, int *mirid1, int *itrig1);
int setTLPmtID(void);


double convertPETimesToTLPE(PETimes *pt);
void initTLElectronics(geofd_dst_common *geo);
int simTLTriggerResponse(PETimes *pt, fraw1_dst_common *fraw1, TLCalibration *tes);


void trig_digit_(TLCalibration *pes, integer4 *mir_id, integer4 *istat);
integer2 fadc_trig_pld_(integer4 *mirid, integer2 *v_pattern, integer2 *h_pattern,
                       integer2 *start, integer2 *stop, integer2 *code);

void mirror_xy_(integer4 *ipmt, integer4 *irow, integer4 *icol);
integer4 jstar1_(integer1 *i1);
integer4 jstar2_(integer2 *i2);
void TLdsp_scan_(double *tau, int madc[1024], int *nbin, int *mped,
               int *pulse, int *itime, int *i1, int *nt,
               double *area, double *time);
// void filter_1_m_(void);

void filterIsolatedHits(void);
int iang_m_(int *f77i, int *f77j);

int tentativeRayleighFilter(void);

// void time_order_m_(void);
void sortHitsByTime(void);
// real4 dot_(real4 a[3],real4 b[3]);
// void unit_norm_(real4 v[3]/*, integer4 *iloc*/);

void initTLWaveformGenerator(void);
void generateTLWaveform(int *ntpe, double tpe[100000], int incr_npe[100000], 
            double *tstart, int *nt, double v[max_tim], double v16[max_tim],
                        TLCalibration *tes);

// void rnpssn_(real4 *a, integer4 *b, integer4 *c);
// void ranlux_(real4 a[], integer4 *b);


#ifdef __cplusplus
}
#endif
#endif
