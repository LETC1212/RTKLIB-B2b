/*------------------------------------------------------------------------------
* ppp_ar_passbypass.h : Pass-by-Pass PPP ambiguity resolution types and API
*
* history : 2025/01/xx  1.0  original
*           2025/03/xx  2.0  add wl_wsum / wl_wxsum to ambarc_t for
*                             inverse-variance weighted WL arc estimation
*-----------------------------------------------------------------------------*/
#ifndef PPP_AR_PASSBYPASS_H
#define PPP_AR_PASSBYPASS_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MAXARC
#define MAXARC  32          /* max continuous arcs per satellite           */
#endif

/* ── ambarc_t : single continuous ambiguity arc ──────────────────────────────
 *
 * [CHANGED v2.0]
 *   N_IF / var_IF  : now set to the LAST epoch's float PPP state value
 *                    (not a weighted mean over the arc).
 *   N_WL / var_WL  : inverse-variance weighted mean of per-epoch HMW values.
 *                    Accumulated via wl_wsum / wl_wxsum (two new fields).
 *
 * [REMOVED]  IF weighted-mean accumulators (implicit via 1/var encoding).
 * [ADDED]    wl_wsum  : Σ(1 / var_wl_epoch)
 * [ADDED]    wl_wxsum : Σ(N_WL_epoch / var_wl_epoch)
 * ──────────────────────────────────────────────────────────────────────────*/
typedef struct {
    int     sat;       /* satellite number                                  */
    int     day;       /* day index (0 = day-0, 1 = day-1)                 */
    gtime_t ts, te;   /* arc start / end epoch                             */
    int     nobs;      /* number of epochs in arc                          */

    /* IF ambiguity: value of the LAST epoch in this arc [m] / [m²]        */
    double  N_IF;
    double  var_IF;

    /* WL ambiguity: inverse-variance weighted mean over arc [cycles/cyc²] */
    double  N_WL;
    double  var_WL;

    /* [NEW] Accumulators for WL inverse-variance weighted mean             */
    double  wl_wsum;   /* Σ(1 / var_wl_epoch)                             */
    double  wl_wxsum;  /* Σ(N_WL_epoch / var_wl_epoch)                    */

    int     fixed_WL;  /* WL integer fixed flag                            */
    int     fixed_NL;  /* NL integer fixed flag                            */
} ambarc_t;

/* ── satamb_t : arc collection per satellite ─────────────────────────────── */
typedef struct {
    int       n;
    ambarc_t  arc[MAXARC];
} satamb_t;

/* ── ddamb_t : double-differenced ambiguity ──────────────────────────────── */
typedef struct {
    int    sat1, sat2;       /* ref / non-ref satellite numbers             */

    double DD_WL;            /* float DD WL [cycles]                        */
    double var_DD_WL;        /* float DD WL variance [cycles²]              */
    double DD_WL_fix;        /* fixed DD WL (integer) stored as double      */
    int    fixed_WL;         /* WL fixed flag                               */

    double DD_IF;            /* float DD IF [m]                             */
    double var_DD_IF;        /* float DD IF variance [m²]                   */

    double DD_NL;            /* float DD NL [cycles]                        */
    double var_DD_NL;        /* float DD NL variance [cycles²]              */
    double DD_NL_fix;        /* fixed DD NL (integer) stored as double      */
    int    fixed_NL;         /* NL fixed flag                               */

    double DD_IF_fix;        /* fixed DD IF [m]                             */

    int    arc1, arc2;       /* day-0 arc indices for ref / sat             */
} ddamb_t;

/* ── Globals defined in ppp_ar_passbypass.c ─────────────────────────────── */
extern satamb_t satamb[MAXSAT];
extern int      n_ddamb;
extern ddamb_t  ddamb[MAXSAT * MAXSAT];
extern int      refsat;
extern int      pbp_resolve_flag;   /* 1 = pass-2 constrained re-processing */
extern int      pbp_base_day_id;    /* day-id of first collected epoch       */

/* ── ppp.c wrapper for static varerr() ──────────────────────────────────── */
extern double pbp_varerr(int sat, int sys, double el, double snr_rover,
                          int f, const prcopt_t *opt, const obsd_t *obs);

/* ── ppp_ar_passbypass.c API ─────────────────────────────────────────────── */
extern void init_arc_data(void);
extern void print_arc_summary(void);
extern int  collect_ambiguities(const rtk_t *rtk, const obsd_t *obs, int n,
                                 int day, satamb_t *satamb);
extern int  compute_dd_ambiguities(const satamb_t *satamb, int refsat,
                                    ddamb_t *ddamb, int *n_dd);
extern int  fix_wl_nl_ambiguities(ddamb_t *ddamb, int n_dd);
extern int  fix_ambiguity(double N_float, double std, double threshold,
                           int *fixed);

/* Ambiguity-subblock NEQ re-solve (paper Eq.18-20) */
extern void pbp_clear_fixed_constraints(void);
extern int  pbp_store_fixed_constraints(const ddamb_t *ddamb, int n_dd,
                                         double Pb);
extern int  pbp_get_fixed_arc_bias(gtime_t t, int sat, double *bias,
                                    double *var);
extern int  pbp_has_fixed_constraints(void);

/* Legacy / link-compatibility stubs */
extern int  apply_ar_fixed(rtk_t *rtk, const ddamb_t *ddamb, int n_dd);
extern int  pbp_apply_session_pseudoobs(rtk_t *rtk);

/* Diagnostics */
extern void pbp_print_epoch_stats(void);
extern int  pbp_epoch_fix_count;
extern int  pbp_epoch_total;
extern int  pbp_constraint_sum;

/* ── ppp_ar_integration.c API ────────────────────────────────────────────── */
extern int  ppp_ar_48h(const prcopt_t *popt, rtk_t *rtk, const obs_t *obs);
extern int  collect_ambiguities_epoch(const rtk_t *rtk, const obsd_t *obs,
                                       int n, int day);

/* driver flags */
extern int  pbp_day_tag;
extern int  pbp_collect_flag;
extern int  pbp_apply_flag;

#ifdef __cplusplus
}
#endif

#endif /* PPP_AR_PASSBYPASS_H */
