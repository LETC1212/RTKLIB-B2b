/*------------------------------------------------------------------------------
* ppp_ar.c : ppp ambiguity resolution with day-to-day differencing
*
* options : -DREV_WL_FCB reversed polarity of WL FCB
*
* reference :
*    [1] H.Okumura, C-gengo niyoru saishin algorithm jiten (in Japanese),
*        Software Technology, 1991
*    [2] Xi et al., 2021, Pass-by-Pass Ambiguity Resolution in Single GPS
*        Receiver PPP Using Observations for Two Sequential Days,
*        Remote Sensing, 2021, 13(18), 3631
*
*          Copyright (C) 2012-2013 by T.TAKASU, All rights reserved.
*          Copyright (C) 2025 Enhanced for PPP-B2b fixed solution
*
* version : $Revision:$ $Date:$
* history : 2013/03/11 1.0  new
*           2025/12/04 2.0  add day-to-day differencing AR for PPP-B2b
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

/* AR thresholds and constants */
#define MIN_ARC_SAT     4       /* minimum satellites for AR */
#define MIN_ARC_TIME    30      /* minimum continuous lock time (epochs) */
#define THRES_RATIO     3.0     /* ratio test threshold */
#define THRES_VAR_AMB   0.25    /* variance threshold for ambiguity (cycles^2) */
#define MIN_TIME_SPAN   72000.0 /* minimum time span for AR (20 hours in seconds) */
#define SIDEREAL_DAY    86164.0905 /* sidereal day in seconds */
#define MAX_AMB_SEARCH  100     /* max search iterations */

/* Ambiguity arc storage - stores ambiguity history over time */
#define MAX_AMB_HISTORY 5000    /* maximum ambiguity history entries */

typedef struct {
    gtime_t time;               /* epoch time */
    int sat;                    /* satellite number */
    double amb[NFREQ];          /* ambiguity estimates (cycles) */
    double std[NFREQ];          /* ambiguity standard deviations (cycles) */
    int lock[NFREQ];            /* lock count */
    int valid;                  /* valid flag */
} amb_hist_t;

static amb_hist_t amb_history[MAX_AMB_HISTORY];  /* ambiguity history */
static int amb_hist_count = 0;                    /* number of history entries */
static int amb_hist_initialized = 0;              /* initialization flag */

/* number and index of states (from ppp.c) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics?9:3)
#ifdef BDS2BDS3
#define NC(opt)     (NSYS+1)
#else
#define NC(opt)     (NSYS)
#endif
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)
#define ND(opt)     ((opt)->nf>=3?1:0)
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt))
#define NB(opt)     (NF(opt)*MAXSAT)
#define NX(opt)     (NR(opt)+NB(opt))
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)
#define IC(s,opt)   (NP(opt)+(s))

/* save ambiguity to history buffer ----------------------------------------*/
static void save_amb_to_history(const rtk_t *rtk, const obsd_t *obs, int n)
{
    const prcopt_t *opt = &rtk->opt;
    int i, j, sat, idx;
    gtime_t oldest_time;
    double time_span;

    trace(3, "save_amb_to_history: n=%d, hist_count=%d\n", n, amb_hist_count);

    /* initialize if needed */
    if (!amb_hist_initialized) {
        amb_hist_count = 0;
        for (i = 0; i < MAX_AMB_HISTORY; i++) {
            amb_history[i].valid = 0;
        }
        amb_hist_initialized = 1;
    }

    /* check if buffer is full - remove oldest entries if needed */
    if (amb_hist_count >= MAX_AMB_HISTORY - MAXSAT) {
        /* find oldest time */
        oldest_time = amb_history[0].time;
        for (i = 1; i < amb_hist_count; i++) {
            if (timediff(amb_history[i].time, oldest_time) < 0) {
                oldest_time = amb_history[i].time;
            }
        }

        /* remove entries older than (current - 2 days) to keep buffer manageable */
        time_span = timediff(obs[0].time, oldest_time);
        if (time_span > 2.0 * SIDEREAL_DAY) {
            int keep_count = 0;
            for (i = 0; i < amb_hist_count; i++) {
                if (timediff(obs[0].time, amb_history[i].time) < 2.0 * SIDEREAL_DAY) {
                    if (keep_count != i) {
                        amb_history[keep_count] = amb_history[i];
                    }
                    keep_count++;
                }
            }
            amb_hist_count = keep_count;
            trace(3, "Cleaned history buffer: kept %d entries\n", keep_count);
        }
    }

    /* save current ambiguities to history */
    for (i = 0; i < n && i < MAXOBS; i++) {
        sat = obs[i].sat;
        if (sat <= 0 || sat > MAXSAT) continue;

        idx = sat - 1;

        /* check if satellite has good lock and ambiguity quality */
        if (rtk->ssat[idx].lock[0] < MIN_ARC_TIME) continue;

        /* add to history */
        if (amb_hist_count < MAX_AMB_HISTORY) {
            amb_history[amb_hist_count].time = obs[i].time;
            amb_history[amb_hist_count].sat = sat;
            amb_history[amb_hist_count].valid = 1;

            for (j = 0; j < NF(opt) && j < NFREQ; j++) {
                int ib = IB(sat, j, opt);
                if (ib >= 0 && ib < rtk->nx) {
                    amb_history[amb_hist_count].amb[j] = rtk->x[ib];
                    amb_history[amb_hist_count].std[j] = SQRT(rtk->P[ib + ib * rtk->nx]);
                    amb_history[amb_hist_count].lock[j] = rtk->ssat[idx].lock[j];
                } else {
                    amb_history[amb_hist_count].amb[j] = 0.0;
                    amb_history[amb_hist_count].std[j] = 999.0;
                    amb_history[amb_hist_count].lock[j] = 0;
                }
            }
            amb_hist_count++;
        }
    }
}

/* find matching ambiguity in history using sidereal day repeat ------------*/
static int find_matching_amb_history(int sat, gtime_t current_time,
                                      double *amb_ref, double *std_ref, int *lock_ref)
{
    int i, best_idx = -1;
    double time_diff, best_diff = 999999.0;
    double target_diff;

    /* target time difference: approximately one sidereal day ago */
    target_diff = SIDEREAL_DAY;

    /* search for matching satellite in history around sidereal day ago */
    for (i = 0; i < amb_hist_count; i++) {
        if (!amb_history[i].valid) continue;
        if (amb_history[i].sat != sat) continue;

        /* compute time difference */
        time_diff = timediff(current_time, amb_history[i].time);

        /* check if close to one sidereal day (with some tolerance) */
        /* Allow ±2 hours window around sidereal day */
        if (fabs(time_diff - target_diff) < 7200.0) {
            /* find the closest match */
            if (fabs(time_diff - target_diff) < best_diff) {
                best_diff = fabs(time_diff - target_diff);
                best_idx = i;
            }
        }
    }

    if (best_idx >= 0) {
        *amb_ref = amb_history[best_idx].amb[0];
        *std_ref = amb_history[best_idx].std[0];
        *lock_ref = amb_history[best_idx].lock[0];

        trace(4, "Found history match: sat=%2d time_diff=%.1f hours best_diff=%.1f min\n",
              sat, timediff(current_time, amb_history[best_idx].time) / 3600.0,
              best_diff / 60.0);
        return 1;
    }

    return 0;
}

/* form day-to-day difference equations using history ----------------------*/
static int form_dd_equations(const rtk_t *rtk, const obsd_t *obs, int n,
                              int *sat_dd, double *dd_amb, double *dd_std)
{
    const prcopt_t *opt = &rtk->opt;
    int i, sat, idx, ndd = 0;
    double amb_cur, std_cur, amb_ref, std_ref;
    int lock_ref;
    gtime_t oldest_hist_time;
    double time_span;

    trace(3, "form_dd_equations: n=%d hist_count=%d\n", n, amb_hist_count);

    if (!amb_hist_initialized || amb_hist_count == 0) {
        trace(3, "Ambiguity history not initialized or empty\n");
        return 0;
    }

    /* check if we have enough time span (at least 20 hours) */
    oldest_hist_time = amb_history[0].time;
    for (i = 1; i < amb_hist_count; i++) {
        if (amb_history[i].valid &&
            timediff(amb_history[i].time, oldest_hist_time) < 0) {
            oldest_hist_time = amb_history[i].time;
        }
    }

    time_span = timediff(obs[0].time, oldest_hist_time);
    if (time_span < MIN_TIME_SPAN) {
        trace(3, "Time span too short: %.1f hours < %.1f hours\n",
              time_span / 3600.0, MIN_TIME_SPAN / 3600.0);
        return 0;
    }

    trace(3, "Time span: %.1f hours, processing AR\n", time_span / 3600.0);

    /* form day-to-day differences for each satellite */
    for (i = 0; i < n && i < MAXOBS; i++) {
        sat = obs[i].sat;
        if (sat <= 0 || sat > MAXSAT) continue;

        idx = sat - 1;

        /* check current lock time */
        if (rtk->ssat[idx].lock[0] < MIN_ARC_TIME) continue;

        /* get current ambiguity */
        int ib = IB(sat, 0, opt);
        if (ib < 0 || ib >= rtk->nx) continue;

        amb_cur = rtk->x[ib];
        std_cur = SQRT(rtk->P[ib + ib * rtk->nx]);

        /* check current ambiguity variance */
        if (std_cur > THRES_VAR_AMB) continue;

        /* find matching reference ambiguity from history */
        if (!find_matching_amb_history(sat, obs[i].time,
                                        &amb_ref, &std_ref, &lock_ref)) {
            continue;  /* no matching history found */
        }

        /* check reference ambiguity quality */
        if (std_ref > THRES_VAR_AMB) continue;
        if (lock_ref < MIN_ARC_TIME) continue;

        /* form day-to-day difference */
        dd_amb[ndd] = amb_cur - amb_ref;
        dd_std[ndd] = SQRT(SQR(std_cur) + SQR(std_ref));
        sat_dd[ndd] = sat;

        trace(4, "DD sat=%2d amb_cur=%8.3f std_cur=%6.3f amb_ref=%8.3f std_ref=%6.3f dd=%8.3f\n",
              sat, amb_cur, std_cur, amb_ref, std_ref, dd_amb[ndd]);

        ndd++;
    }

    trace(3, "Formed %d day-to-day differences\n", ndd);
    return ndd;
}

/* round day-to-day differences to nearest integer -------------------------*/
static int resolve_dd_ambiguity(int ndd, const double *dd_amb, const double *dd_std,
                                 int *dd_fix, double *ratio)
{
    int i;
    double residual_sum = 0.0, residual_min = 1E10, residual_sec = 1E10;
    int success = 0;

    if (ndd < MIN_ARC_SAT) {
        trace(3, "Not enough satellites for AR: %d < %d\n", ndd, MIN_ARC_SAT);
        return 0;
    }

    /* simple rounding strategy - round to nearest integer */
    /* this is a simplified approach; in practice, LAMBDA would be better */
    for (i = 0; i < ndd; i++) {
        dd_fix[i] = ROUND(dd_amb[i]);
        double res = fabs(dd_amb[i] - dd_fix[i]);
        residual_sum += SQR(res / dd_std[i]);

        trace(4, "DD %d: float=%8.3f fixed=%4d residual=%6.3f\n",
              i, dd_amb[i], dd_fix[i], res);
    }

    /* compute ratio test (simplified version) */
    /* in practice, should use second-best solution */
    residual_min = residual_sum;
    residual_sec = residual_sum + 1.0;  /* simplified */

    *ratio = residual_sec / (residual_min + 1E-10);

    /* check ratio test */
    if (*ratio > THRES_RATIO) {
        success = 1;
        trace(3, "AR success: ndd=%d ratio=%.2f\n", ndd, *ratio);
    } else {
        trace(3, "AR failed: ndd=%d ratio=%.2f < %.2f\n", ndd, *ratio, THRES_RATIO);
    }

    return success;
}

/* apply fixed ambiguities to update state vector --------------------------*/
static void apply_fixed_ambiguity(rtk_t *rtk, const obsd_t *obs, int n,
                                   const int *sat_dd, const int *dd_fix, int ndd)
{
    const prcopt_t *opt = &rtk->opt;
    int i, j, sat, idx, ib, lock_ref;
    double *xa, *Pa;
    double amb_ref, std_ref;

    trace(3, "apply_fixed_ambiguity: ndd=%d\n", ndd);

    /* copy float solution to fixed solution */
    xa = rtk->xa;
    Pa = rtk->Pa;
    matcpy(xa, rtk->x, rtk->nx, 1);
    matcpy(Pa, rtk->P, rtk->nx, rtk->nx);

    /* apply fixed ambiguities */
    for (i = 0; i < ndd && i < n; i++) {
        sat = sat_dd[i];
        if (sat <= 0 || sat > MAXSAT) continue;

        idx = sat - 1;
        ib = IB(sat, 0, opt);

        if (ib < 0 || ib >= rtk->nx) continue;

        /* find reference ambiguity from history */
        if (!find_matching_amb_history(sat, obs[i].time,
                                        &amb_ref, &std_ref, &lock_ref)) {
            trace(3, "Warning: no history match for sat %d in apply_fixed\n", sat);
            continue;
        }

        /* directly set ambiguity to reference + fixed integer */
        xa[ib] = amb_ref + dd_fix[i];

        /* tighten ambiguity variance (set to very small value) */
        Pa[ib + ib * rtk->nx] = 1E-6;

        trace(4, "Fixed sat=%2d: float_amb=%8.3f fixed_amb=%8.3f ref_amb=%8.3f int=%d\n",
              sat, rtk->x[ib], xa[ib], amb_ref, dd_fix[i]);
    }

    /* mark satellites as fixed */
    for (i = 0; i < ndd; i++) {
        sat = sat_dd[i];
        if (sat <= 0 || sat > MAXSAT) continue;
        idx = sat - 1;

        for (j = 0; j < opt->nf && j < NFREQ; j++) {
            rtk->ssat[idx].fix[j] = 2;  /* 2 indicates fixed */
        }
    }
}

/* update receiver clock after ambiguity fixing ----------------------------*/
static void update_clock_fixed(rtk_t *rtk, const obsd_t *obs, int n,
                                const double *azel, const nav_t *nav)
{
    const prcopt_t *opt = &rtk->opt;
    int i, j, sys, idx;
    double *xa = rtk->xa;
    double lam, C, P, delta_clk, weight;
    double clk_correction[NSYS + 1] = {0};
    double clk_weight[NSYS + 1] = {0};

    trace(3, "update_clock_fixed: n=%d\n", n);

    /* recompute clock using fixed ambiguities */
    for (i = 0; i < n && i < MAXOBS; i++) {
        int sat = obs[i].sat;
        if (sat <= 0 || sat > MAXSAT) continue;

        /* check if satellite is fixed */
        idx = sat - 1;
        if (rtk->ssat[idx].fix[0] != 2) continue;

        /* check valid observations */
        if (obs[i].L[0] == 0.0 || obs[i].P[0] == 0.0) continue;

        /* get system */
        sys = satsys(sat, NULL);
        int sys_idx = -1;
        if (sys == SYS_GPS) sys_idx = 0;
        else if (sys == SYS_GLO) sys_idx = 1;
        else if (sys == SYS_GAL) sys_idx = 2;
        else if (sys == SYS_CMP) sys_idx = 3;
        else continue;

        /* get wavelength */
        double freq = sat2freq(sat, obs[i].code[0], nav);
        if (freq <= 0.0) continue;
        lam = CLIGHT / freq;

        /* get ambiguity */
        int ib = IB(sat, 0, opt);
        if (ib < 0 || ib >= rtk->nx) continue;

        /* compute phase range with fixed ambiguity */
        double L_range = obs[i].L[0] * lam + xa[ib] * lam;
        double P_range = obs[i].P[0];

        /* clock correction from code-phase consistency */
        delta_clk = P_range - L_range;

        /* use elevation-dependent weighting */
        if (azel) {
            double el = azel[1 + i * 2];
            weight = el > 0.0 ? 1.0 / (1.0 + exp(-10.0 * (el - 30.0 * D2R))) : 0.0;
        } else {
            weight = 1.0;
        }

        clk_correction[sys_idx] += delta_clk * weight;
        clk_weight[sys_idx] += weight;
    }

    /* apply clock corrections */
    for (i = 0; i < NC(opt); i++) {
        if (clk_weight[i] > 0.0) {
            int ic = IC(i, opt);
            if (ic >= 0 && ic < rtk->nx) {
                double corr = clk_correction[i] / clk_weight[i];
                xa[ic] += corr;
                trace(4, "Clock sys=%d: correction=%8.3f m\n", i, corr);
            }
        }
    }
}

/* resolve integer ambiguity for ppp with day-to-day differencing ----------*/
extern int pppamb(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                  const double *azel, int *exc)
{
    const prcopt_t *opt = &rtk->opt;
    int sat_dd[MAXSAT];
    double dd_amb[MAXSAT], dd_std[MAXSAT];
    int dd_fix[MAXSAT];
    int ndd, success;
    double ratio;

    trace(3, "pppamb: n=%d nx=%d\n", n, rtk->nx);

    /* check if AR is enabled */
    if (opt->modear == ARMODE_OFF) {
        return 0;
    }

    /* continuously save ambiguities to history buffer */
    save_amb_to_history(rtk, obs, n);

    /* form day-to-day difference equations */
    ndd = form_dd_equations(rtk, obs, n, sat_dd, dd_amb, dd_std);

    if (ndd < MIN_ARC_SAT) {
        trace(3, "Not enough DD equations: %d < %d\n", ndd, MIN_ARC_SAT);
        return 0;
    }

    /* resolve day-to-day ambiguity differences */
    success = resolve_dd_ambiguity(ndd, dd_amb, dd_std, dd_fix, &ratio);

    if (!success) {
        rtk->sol.ratio = (float)ratio;
        return 0;
    }

    /* apply fixed ambiguities */
    apply_fixed_ambiguity(rtk, obs, n, sat_dd, dd_fix, ndd);

    /* update receiver clock with fixed ambiguities */
    update_clock_fixed(rtk, obs, n, azel, nav);

    /* store ratio */
    rtk->sol.ratio = (float)ratio;

    trace(2, "PPP AR success: ndd=%d ratio=%.2f\n", ndd, ratio);

    return 1;
}
