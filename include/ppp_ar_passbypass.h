/*------------------------------------------------------------------------------
 * ppp_ar_passbypass.h : Pass-by-Pass PPP ambiguity resolution header
 *
 * Type definitions: ambarc_t, satamb_t, ddamb_t
 * Function declarations for ppp_ar_passbypass.c / ppp_ar_integration.c
 *-----------------------------------------------------------------------------*/
#ifndef PPP_AR_PASSBYPASS_H
#define PPP_AR_PASSBYPASS_H

#include "rtklib.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ── Constants ─────────────────────────────────────────────────────────── */
#ifndef MAXARC
#define MAXARC 64   /* max arcs per satellite */
#endif

/* ── Arc ambiguity data ────────────────────────────────────────────────── */
typedef struct {
    int    sat;           /* satellite number                             */
    int    day;           /* day index (0 or 1)                           */
    gtime_t ts, te;       /* arc start/end time                           */
    int    nobs;          /* number of observations in arc                */
    double N_IF;          /* IF ambiguity estimate (metres)               */
    double var_IF;        /* variance of IF ambiguity                     */
    double N_WL;          /* WL ambiguity estimate (cycles)               */
    double var_WL;        /* variance of WL ambiguity                     */
    double wl_wsum;       /* WL inverse-variance weight sum               */
    double wl_wxsum;      /* WL weighted value sum                        */
    int    fixed_WL;      /* WL fixed flag                                */
    int    fixed_NL;      /* NL fixed flag                                */
} ambarc_t;

typedef struct {
    int      n;           /* number of arcs                               */
    ambarc_t arc[MAXARC]; /* arc data array                               */
} satamb_t;

/* ── DD ambiguity data ─────────────────────────────────────────────────── */
typedef struct {
    int    sat1, sat2;    /* reference and rover satellite numbers        */
    int    arc1, arc2;    /* arc indices (day0)                           */
    double DD_WL;         /* DD WL ambiguity (cycles)                     */
    double var_DD_WL;     /* DD WL variance                               */
    double DD_WL_fix;     /* DD WL fixed integer value                    */
    int    fixed_WL;      /* DD WL fixed flag                             */
    double DD_NL;         /* DD NL ambiguity (cycles)                     */
    double var_DD_NL;     /* DD NL variance                               */
    double DD_NL_fix;     /* DD NL fixed integer value                    */
    int    fixed_NL;      /* DD NL fixed flag                             */
    double DD_IF;         /* DD IF ambiguity (metres)                     */
    double var_DD_IF;     /* DD IF variance                               */
    double DD_IF_fix;     /* DD IF fixed value (metres)                   */
} ddamb_t;

/* ── Global variables (defined in ppp_ar_passbypass.c) ─────────────────── */
extern satamb_t satamb[MAXSAT];
extern int      n_ddamb;
extern ddamb_t  ddamb[MAXSAT * MAXSAT];
extern int      refsat;
extern int      pbp_base_day_id;

/* Driver flags (defined in ppp_ar_integration.c) */
extern int pbp_day_tag;
extern int pbp_collect_flag;
extern int pbp_apply_flag;

/* NEQ/resolve flags (defined in ppp_ar_passbypass.c) */
extern int pbp_resolve_flag;
extern int pbp_neq_accum_flag;
extern int pbp_epoch_collected;
extern int pbp_current_day;

/* Epoch counters (defined in ppp_ar_passbypass.c) */
extern int pbp_epoch_fix_count;
extern int pbp_epoch_total;
extern int pbp_constraint_sum;

/* Day windows (defined in ppp_ar_passbypass.c) */
extern gtime_t pbp_day_start_win[2];
extern gtime_t pbp_day_end_win[2];
extern int     pbp_epoch_offset[2];
extern int     pbp_day_epoch_n[2];

/* ── Functions: ppp_ar_passbypass.c ────────────────────────────────────── */

/* Arc data management */
extern void init_arc_data(void);
extern void print_arc_summary(void);
extern int  collect_ambiguities(const rtk_t *rtk, const obsd_t *obs,
                                 int n, int day, satamb_t *satamb);

/* DD ambiguity computation and fixing */
extern int  compute_dd_ambiguities(const satamb_t *satamb, int refsat,
                                    ddamb_t *ddamb, int *n_dd);
extern int  fix_wl_nl_ambiguities(ddamb_t *ddamb, int n_dd);
extern int  fix_ambiguity(double N_float, double std, double threshold,
                           int *fixed);

/* BDS satellite sidereal filter */
extern int  pbp_bds_is_sidereal(int sat);

/* NEQ system */
extern int  pbp_neq_init(gtime_t t0, gtime_t t1, double ti,
                           const prcopt_t *opt);
extern int  pbp_neq_add_epoch(rtk_t *rtk, const obsd_t *obs, int n,
                                const double *v, const double *H,
                                const double *R, int nv);
extern int  pbp_build_arc_columns(void);
extern int  pbp_store_fixed_constraints(const ddamb_t *dd, int n_dd,
                                         double Pb);
extern int  pbp_finalize_final_neq(void);
extern void pbp_clear_fixed_constraints(void);
extern int  pbp_has_fixed_constraints(void);
extern void pbp_set_day_window(int day, gtime_t ts, gtime_t te, double ti);

/* Results */
extern int  pbp_write_day1_fixed_clock_file(const char *path);
extern int  pbp_get_fixed_clock(gtime_t t, int sys_idx, double *clk);
extern int  pbp_get_fixed_arc_bias(gtime_t t, int sat, double *bias,
                                    double *var);

/* Legacy stubs */
extern int  pbp_apply_session_pseudoobs(rtk_t *rtk);
extern int  apply_ar_fixed(rtk_t *rtk, const ddamb_t *ddamb, int n_dd);

/* Diagnostics */
extern void pbp_print_epoch_stats(void);

/* ── Functions: ppp_ar_integration.c ───────────────────────────────────── */
extern int  ppp_ar_48h(const prcopt_t *popt, rtk_t *rtk, const obs_t *obs);
extern int  collect_ambiguities_epoch(const rtk_t *rtk, const obsd_t *obs,
                                       int n, int day);

/* ── Function: ppp.c (varerr wrapper) ──────────────────────────────────── */
extern double pbp_varerr(int sat, int sys, double el, double snr_rover,
                          int f, const prcopt_t *opt, const obsd_t *obs);

#ifdef __cplusplus
}
#endif

#endif /* PPP_AR_PASSBYPASS_H */
