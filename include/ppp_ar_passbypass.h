/*------------------------------------------------------------------------------
 * ppp_ar_passbypass.h : Pass-by-Pass PPP ambiguity resolution header
 *-----------------------------------------------------------------------------*/
#ifndef PPP_AR_PASSBYPASS_H
#define PPP_AR_PASSBYPASS_H

#include "rtklib.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MAXARC
#define MAXARC 64
#endif

typedef struct {
    int    sat;
    int    day;
    gtime_t ts, te;
    int    nobs;
    double N_IF;
    double var_IF;
    double N_WL;
    double var_WL;
    double wl_wsum;
    double wl_wxsum;
    int    fixed_WL;
    int    fixed_NL;
} ambarc_t;

typedef struct {
    int      n;
    ambarc_t arc[MAXARC];
} satamb_t;

typedef struct {
    int    sat1, sat2;
    int    arc1, arc2;
    double DD_WL;
    double var_DD_WL;
    double DD_WL_fix;
    int    fixed_WL;
    double DD_NL;
    double var_DD_NL;
    double DD_NL_fix;
    int    fixed_NL;
    double DD_IF;
    double var_DD_IF;
    double DD_IF_fix;
} ddamb_t;

extern satamb_t satamb[MAXSAT];
extern int      n_ddamb;
extern ddamb_t  ddamb[MAXSAT * MAXSAT];
extern int      refsat;
extern int      pbp_base_day_id;

extern int pbp_day_tag;
extern int pbp_collect_flag;
extern int pbp_apply_flag;
extern int pbp_resolve_flag;
extern int pbp_neq_accum_flag;
extern int pbp_epoch_collected;
extern int pbp_current_day;

extern int pbp_epoch_fix_count;
extern int pbp_epoch_total;
extern int pbp_constraint_sum;

extern gtime_t pbp_day_start_win[2];
extern gtime_t pbp_day_end_win[2];
extern int     pbp_epoch_offset[2];
extern int     pbp_day_epoch_n[2];

extern void init_arc_data(void);
extern void print_arc_summary(void);
extern int  collect_ambiguities(const rtk_t *rtk, const obsd_t *obs,
                                int n, int day, satamb_t *satamb);
extern int  compute_dd_ambiguities(const satamb_t *satamb, int refsat,
                                   ddamb_t *ddamb, int *n_dd);
extern int  fix_wl_nl_ambiguities(ddamb_t *ddamb, int n_dd);
extern int  fix_ambiguity(double N_float, double std, double threshold,
                          int *fixed);
extern int  pbp_bds_is_sidereal(int sat);

extern int  pbp_neq_init(gtime_t t0, gtime_t t1, double ti,
                         const prcopt_t *opt);
extern int  pbp_neq_add_epoch(rtk_t *rtk, const obsd_t *obs, int n,
                              const double *v, const double *H,
                              const double *R, int nv,
                              const double *x_lin);
extern int  pbp_build_arc_columns(void);
extern int  pbp_store_fixed_constraints(const ddamb_t *dd, int n_dd,
                                        double Pb);
extern int  pbp_finalize_final_neq(void);
extern void pbp_clear_fixed_constraints(void);
extern int  pbp_has_fixed_constraints(void);
extern void pbp_set_day_window(int day, gtime_t ts, gtime_t te, double ti);

extern int  pbp_write_day1_fixed_clock_file(const char *path);
extern int  pbp_write_day1_neqfloat_clock_file(const char *path);
extern int  pbp_get_fixed_clock(gtime_t t, int sys_idx, double *clk);
extern int  pbp_get_fixed_arc_bias(gtime_t t, int sat, double *bias,
                                   double *var);

extern int  pbp_apply_session_pseudoobs(rtk_t *rtk);
extern int  apply_ar_fixed(rtk_t *rtk, const ddamb_t *ddamb, int n_dd);
extern void pbp_print_epoch_stats(void);

extern int  ppp_ar_48h(const prcopt_t *popt, rtk_t *rtk, const obs_t *obs);
extern int  collect_ambiguities_epoch(const rtk_t *rtk, const obsd_t *obs,
                                      int n, int day);

extern double pbp_varerr(int sat, int sys, double el, double snr_rover,
                         int f, const prcopt_t *opt, const obsd_t *obs);

#ifdef __cplusplus
}
#endif

#endif /* PPP_AR_PASSBYPASS_H */
