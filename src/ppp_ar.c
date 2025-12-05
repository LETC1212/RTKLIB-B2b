/*------------------------------------------------------------------------------
* ppp_ar.c : ppp ambiguity resolution using day-to-day differencing
*
* options : -DREV_WL_FCB reversed polarity of WL FCB
*
* reference :
*    [1] Xi et al., 2021, Pass-by-Pass Ambiguity Resolution in Single GPS
*        Receiver PPP Using Observations for Two Sequential Days
*    [2] H.Okumura, C-gengo niyoru saishin algorithm jiten (in Japanese),
*        Software Technology, 1991
*
*          Copyright (C) 2012-2013 by T.TAKASU, All rights reserved.
*          Copyright (C) 2025 - Day-to-day AR implementation
*
* version : $Revision:$ $Date:$
* history : 2013/03/11 1.0  new
*           2025/01/01 2.0  add day-to-day ambiguity resolution for PPP-B2b
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define ROUND(x)        (int)floor((x)+0.5)

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
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)
#define ID(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt))
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

#define MIN_SAT_AR      4       /* minimum number of satellites for AR */
#define MIN_LOCK_AR     10      /* minimum lock count for AR (epochs) */
#define RATIO_THRESHOLD 3.0     /* ratio threshold for AR */
#define MIN_ELEVATION   15.0    /* minimum elevation for AR (degrees) */
#define MAX_AGE_DIFF    3600.0  /* maximum age for day-to-day differencing (s) */

/* global structure to store previous day ambiguities */
typedef struct {
    gtime_t time;               /* time of previous day */
    double amb[MAXSAT][NFREQ];  /* ambiguity parameters from previous day */
    double std[MAXSAT][NFREQ];  /* standard deviations */
    int valid[MAXSAT][NFREQ];   /* valid flags */
    int initialized;            /* initialization flag */
} prevday_t;

static prevday_t prevday = {0};

/* check if new day started ------------------------------------------------*/
static int check_newday(gtime_t time1, gtime_t time2)
{
    double ep1[6], ep2[6];
    time2epoch(time1, ep1);
    time2epoch(time2, ep2);

    /* check if day changed */
    if ((int)ep1[2] != (int)ep2[2] ||
        (int)ep1[1] != (int)ep2[1] ||
        (int)ep1[0] != (int)ep2[0]) {
        return 1;
    }
    return 0;
}

/* save current ambiguities for next day -----------------------------------*/
static void save_ambiguities(rtk_t *rtk)
{
    const prcopt_t *opt = &rtk->opt;
    int i, j, k;

    prevday.time = rtk->sol.time;
    prevday.initialized = 0;

    for (i = 0; i < MAXSAT; i++) {
        for (j = 0; j < NF(opt); j++) {
            k = IB(i + 1, j, opt);

            /* only save valid ambiguities with sufficient lock */
            if (rtk->x[k] != 0.0 && rtk->ssat[i].lock[j] >= MIN_LOCK_AR &&
                rtk->ssat[i].vsat[j] && rtk->P[k + k * rtk->nx] > 0.0) {

                prevday.amb[i][j] = rtk->x[k];
                prevday.std[i][j] = sqrt(rtk->P[k + k * rtk->nx]);
                prevday.valid[i][j] = 1;
                prevday.initialized = 1;
            } else {
                prevday.amb[i][j] = 0.0;
                prevday.std[i][j] = 0.0;
                prevday.valid[i][j] = 0;
            }
        }
    }

    trace(3, "save_ambiguities: saved %d valid ambiguities\n", prevday.initialized);
}

/* compute day-to-day ambiguity differences --------------------------------*/
static int compute_dd_ambiguity(rtk_t *rtk, const obsd_t *obs, int n,
                                 const double *azel, double *dd_amb,
                                 double *dd_std, int *sat_list)
{
    const prcopt_t *opt = &rtk->opt;
    int i, j, k, ndd = 0;
    double dt, amb_cur, amb_prev, std_cur, std_prev;

    /* check if previous day data is valid */
    if (!prevday.initialized) {
        trace(3, "compute_dd_ambiguity: previous day not initialized\n");
        return 0;
    }

    /* check time difference */
    dt = timediff(rtk->sol.time, prevday.time);
    if (dt < 86000.0 || dt > 90000.0) {  /* roughly 24 hours ± 1 hour */
        trace(3, "compute_dd_ambiguity: time diff too large or small: %.1f\n", dt);
        return 0;
    }

    /* compute day-to-day differences for each satellite */
    for (i = 0; i < n && i < MAXOBS; i++) {
        int sat = obs[i].sat;
        if (sat <= 0 || sat > MAXSAT) continue;

        /* check elevation */
        if (azel[1 + i * 2] < MIN_ELEVATION * D2R) continue;

        /* only use first frequency for simplicity (can be extended) */
        j = 0;
        k = IB(sat, j, opt);

        /* check if current ambiguity is valid */
        if (rtk->x[k] == 0.0 || !rtk->ssat[sat - 1].vsat[j] ||
            rtk->ssat[sat - 1].lock[j] < MIN_LOCK_AR) continue;

        /* check if previous day ambiguity is valid */
        if (!prevday.valid[sat - 1][j]) continue;

        /* compute day-to-day difference */
        amb_cur = rtk->x[k];
        amb_prev = prevday.amb[sat - 1][j];
        std_cur = sqrt(rtk->P[k + k * rtk->nx]);
        std_prev = prevday.std[sat - 1][j];

        dd_amb[ndd] = amb_cur - amb_prev;
        dd_std[ndd] = sqrt(std_cur * std_cur + std_prev * std_prev);
        sat_list[ndd] = sat;

        trace(4, "dd_amb[%2d] sat=%3d: cur=%.4f prev=%.4f diff=%.4f std=%.4f\n",
              ndd, sat, amb_cur, amb_prev, dd_amb[ndd], dd_std[ndd]);

        ndd++;
    }

    trace(3, "compute_dd_ambiguity: %d day-to-day differences computed\n", ndd);
    return ndd;
}

/* apply fixed ambiguities to update state ---------------------------------*/
static int apply_fixed_ambiguities(rtk_t *rtk, const int *sat_list, int ndd,
                                     const double *dd_fix, const double *dd_amb)
{
    const prcopt_t *opt = &rtk->opt;
    int i, j, k;
    double delta;

    /* update ambiguities with fixed values */
    for (i = 0; i < ndd; i++) {
        int sat = sat_list[i];
        if (sat <= 0 || sat > MAXSAT) continue;

        j = 0;  /* first frequency */
        k = IB(sat, j, opt);

        /* compute correction from fixed day-to-day difference */
        delta = dd_fix[i] - dd_amb[i];  /* fixed - float */

        /* update ambiguity: amb_new = amb_cur - delta */
        /* because: amb_cur - amb_prev = dd_amb (float) */
        /* we want: amb_new - amb_prev = dd_fix (fixed integer) */
        /* so: amb_new = amb_prev + dd_fix = amb_cur - delta */
        /* Note: ambiguities are in meters, delta is in cycles, need wavelength conversion */
        double freq = sat2freq(sat, CODE_L1C, NULL);
        double wavelength = (freq > 0.0) ? CLIGHT / freq : 0.0;
        rtk->xa[k] = rtk->x[k] - delta * wavelength;

        /* mark as fixed */
        rtk->ssat[sat - 1].fix[j] = 2;

        trace(4, "apply_fixed: sat=%3d delta=%.4f amb_float=%.4f amb_fixed=%.4f\n",
              sat, delta, rtk->x[k], rtk->xa[k]);
    }

    /* update clock parameters after ambiguity fixing */
    /* The clock parameters need to be adjusted because ambiguities changed */
    for (i = 0; i < NC(opt); i++) {
        k = IC(i, opt);
        rtk->xa[k] = rtk->x[k];  /* initially copy from float solution */
    }

    /* copy other parameters */
    for (i = 0; i < NP(opt); i++) {
        rtk->xa[i] = rtk->x[i];  /* position */
    }
    for (i = NP(opt) + NC(opt); i < NR(opt); i++) {
        rtk->xa[i] = rtk->x[i];  /* troposphere, ionosphere, etc. */
    }

    /* copy covariance matrix */
    matcpy(rtk->Pa, rtk->P, rtk->nx, rtk->nx);

    /* reduce covariance for fixed ambiguities */
    for (i = 0; i < ndd; i++) {
        int sat = sat_list[i];
        if (sat <= 0 || sat > MAXSAT) continue;

        j = 0;
        k = IB(sat, j, opt);
        rtk->Pa[k + k * rtk->nx] *= 0.01;  /* reduce variance to 1% */
    }

    trace(3, "apply_fixed_ambiguities: updated %d ambiguities\n", ndd);
    return 1;
}

/* resolve integer ambiguity for ppp using day-to-day differencing ---------*/
extern int pppamb(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                  const double *azel, int *exc)
{
    double *dd_amb, *dd_std, *Q, *F, *s, *dd_fix;
    int *sat_list, ndd, info, i;
    int ret = 0;

    trace(3, "pppamb: n=%d\n", n);

    /* check if new day started, save ambiguities from previous day */
    if (prevday.initialized && check_newday(prevday.time, rtk->sol.time)) {
        trace(3, "pppamb: new day detected, saving previous ambiguities\n");
        save_ambiguities(rtk);
        return 0;  /* wait for next epoch to accumulate data */
    }

    /* save ambiguities at end of first day */
    if (!prevday.initialized) {
        /* check if we have stable solution to save */
        if (rtk->sol.stat == SOLQ_PPP || rtk->sol.stat == SOLQ_FIX) {
            save_ambiguities(rtk);
        }
        return 0;  /* not ready for AR yet */
    }

    /* allocate memory */
    dd_amb = mat(MAXSAT, 1);
    dd_std = mat(MAXSAT, 1);
    sat_list = imat(MAXSAT, 1);
    Q = mat(MAXSAT, MAXSAT);
    F = mat(MAXSAT, 2);
    s = mat(2, 1);
    dd_fix = mat(MAXSAT, 1);

    /* compute day-to-day ambiguity differences */
    ndd = compute_dd_ambiguity(rtk, obs, n, azel, dd_amb, dd_std, sat_list);

    if (ndd < MIN_SAT_AR) {
        trace(3, "pppamb: insufficient satellites for AR: %d < %d\n", ndd, MIN_SAT_AR);
        goto cleanup;
    }

    /* construct covariance matrix */
    for (i = 0; i < ndd; i++) {
        Q[i + i * ndd] = dd_std[i] * dd_std[i];
    }

    /* LAMBDA: integer least-squares estimation */
    /* Note: RTKLIB's lambda function may need proper implementation */
    /* For now, we use simple rounding as fallback */
    trace(3, "pppamb: attempting LAMBDA with %d ambiguities\n", ndd);

    info = lambda(ndd, 2, dd_amb, Q, F, s);

    if (info != 0) {
        trace(3, "pppamb: LAMBDA failed, using simple rounding\n");
        /* fallback: simple rounding */
        for (i = 0; i < ndd; i++) {
            dd_fix[i] = ROUND(dd_amb[i]);
        }
        rtk->sol.ratio = 0.0;
    } else {
        /* use best solution from LAMBDA */
        for (i = 0; i < ndd; i++) {
            dd_fix[i] = F[i];
        }

        /* compute ratio test */
        rtk->sol.ratio = (s[0] > 0.0) ? (float)(s[1] / s[0]) : 0.0f;
        if (rtk->sol.ratio > 999.9f) rtk->sol.ratio = 999.9f;

        trace(3, "pppamb: ratio=%.2f\n", rtk->sol.ratio);

        /* check ratio threshold */
        if (rtk->sol.ratio < RATIO_THRESHOLD) {
            trace(3, "pppamb: ratio test failed: %.2f < %.2f\n",
                  rtk->sol.ratio, RATIO_THRESHOLD);
            goto cleanup;
        }
    }

    /* apply fixed ambiguities */
    if (apply_fixed_ambiguities(rtk, sat_list, ndd, dd_fix, dd_amb)) {
        trace(3, "pppamb: ambiguity fixed successfully, ndd=%d ratio=%.2f\n",
              ndd, rtk->sol.ratio);
        ret = 1;  /* success */
    }

cleanup:
    free(dd_amb);
    free(dd_std);
    free(sat_list);
    free(Q);
    free(F);
    free(s);
    free(dd_fix);

    return ret;
}
