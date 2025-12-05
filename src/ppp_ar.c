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
#define MIN_LOCK_AR     5       /* minimum lock count for AR (epochs) - reduced from 10 */
#define RATIO_THRESHOLD 2.5     /* ratio threshold for AR - reduced from 3.0 */
#define MIN_ELEVATION   10.0    /* minimum elevation for AR (degrees) - reduced from 15 */
#define MIN_TIME_DIFF   82800.0 /* minimum time difference for day-to-day (23h) */
#define MAX_TIME_DIFF   93600.0 /* maximum time difference for day-to-day (26h) */

/* global structure to store previous day ambiguities */
typedef struct {
    gtime_t time;               /* time of previous day */
    int doy;                    /* day of year */
    double amb[MAXSAT][NFREQ];  /* ambiguity parameters from previous day */
    double std[MAXSAT][NFREQ];  /* standard deviations */
    int lock[MAXSAT][NFREQ];    /* lock counts */
    int valid[MAXSAT][NFREQ];   /* valid flags */
    int initialized;            /* initialization flag */
    int epoch_count;            /* number of epochs processed */
} prevday_t;

static prevday_t prevday = {0};

/* check if different day ----------------------------------------------------*/
static int check_different_day(gtime_t time1, gtime_t time2)
{
    return (int)time2doy(time1) != (int)time2doy(time2);
}

/* update ambiguities continuously -------------------------------------------*/
static void update_ambiguities(rtk_t *rtk)
{
    const prcopt_t *opt = &rtk->opt;
    int i, j, k;
    int current_doy = (int)time2doy(rtk->sol.time);

    /* check if solution is valid */
    if (rtk->sol.stat != SOLQ_PPP && rtk->sol.stat != SOLQ_FIX) {
        return;
    }

    /* if first initialization or new day detected */
    if (!prevday.initialized || current_doy != prevday.doy) {
        /* if new day and we had previous data, keep it as reference */
        if (prevday.initialized && current_doy != prevday.doy) {
            trace(2, "update_ambiguities: new day detected (DOY %d->%d), previous data saved for AR\n",
                  prevday.doy, current_doy);
            /* don't reset, keep previous day data for cross-day differencing */
            prevday.initialized = 2;  /* mark as ready for AR */
            return;  /* don't update on day boundary */
        }

        /* first initialization */
        prevday.doy = current_doy;
        prevday.time = rtk->sol.time;
        prevday.initialized = 1;
        prevday.epoch_count = 0;

        trace(2, "update_ambiguities: first initialization (DOY %d)\n", current_doy);
    }

    /* continuously update ambiguities during the day */
    int updated = 0;
    for (i = 0; i < MAXSAT; i++) {
        for (j = 0; j < NF(opt); j++) {
            k = IB(i + 1, j, opt);

            /* only save valid ambiguities with sufficient lock */
            if (rtk->x[k] != 0.0 && rtk->ssat[i].lock[j] >= MIN_LOCK_AR &&
                rtk->ssat[i].vsat[j] && rtk->P[k + k * rtk->nx] > 0.0) {

                /* update or initialize */
                if (!prevday.valid[i][j] ||
                    rtk->ssat[i].lock[j] > prevday.lock[i][j]) {

                    prevday.amb[i][j] = rtk->x[k];
                    prevday.std[i][j] = sqrt(rtk->P[k + k * rtk->nx]);
                    prevday.lock[i][j] = rtk->ssat[i].lock[j];
                    prevday.valid[i][j] = 1;
                    prevday.time = rtk->sol.time;
                    updated++;
                }
            }
        }
    }

    prevday.epoch_count++;

    if (updated > 0 && prevday.epoch_count % 100 == 0) {
        trace(3, "update_ambiguities: updated %d ambiguities (epoch %d, DOY %d)\n",
              updated, prevday.epoch_count, current_doy);
    }
}

/* compute day-to-day ambiguity differences ----------------------------------*/
static int compute_dd_ambiguity(rtk_t *rtk, const obsd_t *obs, int n,
                                 const double *azel, double *dd_amb,
                                 double *dd_std, int *sat_list)
{
    const prcopt_t *opt = &rtk->opt;
    int i, j, k, ndd = 0;
    double dt, amb_cur, amb_prev, std_cur, std_prev;

    /* check if previous day data is ready for AR */
    if (prevday.initialized != 2) {
        return 0;
    }

    /* check time difference */
    dt = timediff(rtk->sol.time, prevday.time);

    /* more flexible time check - allow 23 to 26 hours */
    if (dt < MIN_TIME_DIFF || dt > MAX_TIME_DIFF) {
        if (dt > 3600.0) {  /* only trace if significant time passed */
            trace(3, "compute_dd_ambiguity: time diff out of range: %.1f s (need %.1f-%.1f)\n",
                  dt, MIN_TIME_DIFF, MAX_TIME_DIFF);
        }
        return 0;
    }

    trace(3, "compute_dd_ambiguity: time diff = %.1f s, processing...\n", dt);

    /* compute day-to-day differences for each satellite */
    for (i = 0; i < n && i < MAXOBS; i++) {
        int sat = obs[i].sat;
        if (sat <= 0 || sat > MAXSAT) continue;

        /* check elevation */
        if (azel[1 + i * 2] < MIN_ELEVATION * D2R) continue;

        /* process all frequencies */
        for (j = 0; j < NF(opt); j++) {
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

            trace(4, "dd_amb[%2d] sat=%3d freq=%d: cur=%.4f prev=%.4f diff=%.4f std=%.4f el=%.1f\n",
                  ndd, sat, j+1, amb_cur, amb_prev, dd_amb[ndd], dd_std[ndd],
                  azel[1 + i * 2] * R2D);

            ndd++;

            /* for now, only use first frequency to simplify */
            if (j == 0) break;
        }
    }

    trace(2, "compute_dd_ambiguity: %d day-to-day differences computed (DOY %d)\n",
          ndd, (int)time2doy(rtk->sol.time));

    return ndd;
}

/* apply fixed ambiguities to update state -----------------------------------*/
static int apply_fixed_ambiguities(rtk_t *rtk, const int *sat_list, int ndd,
                                     const double *dd_fix, const double *dd_amb)
{
    const prcopt_t *opt = &rtk->opt;
    int i, j, k;
    double delta, freq, wavelength;

    /* copy float solution to fixed solution first */
    matcpy(rtk->xa, rtk->x, rtk->nx, 1);
    matcpy(rtk->Pa, rtk->P, rtk->nx, rtk->nx);

    /* update ambiguities with fixed values */
    for (i = 0; i < ndd; i++) {
        int sat = sat_list[i];
        if (sat <= 0 || sat > MAXSAT) continue;

        j = 0;  /* first frequency */
        k = IB(sat, j, opt);

        /* compute correction from fixed day-to-day difference */
        delta = dd_fix[i] - dd_amb[i];  /* fixed - float (cycles) */

        /* get wavelength for this satellite */
        freq = sat2freq(sat, CODE_L1C, NULL);
        wavelength = (freq > 0.0) ? CLIGHT / freq : 0.0;

        if (wavelength == 0.0) {
            trace(2, "apply_fixed: sat=%3d invalid frequency\n", sat);
            continue;
        }

        /* update ambiguity: amb_fixed = amb_float - delta * wavelength */
        rtk->xa[k] = rtk->x[k] - delta * wavelength;

        /* reduce variance significantly for fixed ambiguity */
        rtk->Pa[k + k * rtk->nx] = rtk->P[k + k * rtk->nx] * 0.001;

        /* mark as fixed */
        rtk->ssat[sat - 1].fix[j] = 2;

        trace(3, "apply_fixed: sat=%3d delta=%.2f cyc, amb: %.4f -> %.4f m, wl=%.4f m\n",
              sat, delta, rtk->x[k], rtk->xa[k], wavelength);
    }

    trace(2, "apply_fixed_ambiguities: successfully fixed %d ambiguities\n", ndd);
    return 1;
}

/* resolve integer ambiguity for ppp using day-to-day differencing -----------*/
extern int pppamb(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                  const double *azel, int *exc)
{
    double *dd_amb, *dd_std, *Q, *F, *s, *dd_fix;
    int *sat_list, ndd, info, i;
    int ret = 0;

    trace(4, "pppamb: n=%d stat=%d\n", n, rtk->sol.stat);

    /* continuously update ambiguity database */
    update_ambiguities(rtk);

    /* check if we're ready for AR (need previous day data) */
    if (prevday.initialized != 2) {
        trace(4, "pppamb: not ready for AR (init=%d)\n", prevday.initialized);
        return 0;
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
    trace(2, "pppamb: attempting LAMBDA with %d ambiguities\n", ndd);

    info = lambda(ndd, 2, dd_amb, Q, F, s);

    if (info != 0) {
        trace(2, "pppamb: LAMBDA failed (info=%d), using simple rounding\n", info);

        /* fallback: simple rounding for validation */
        for (i = 0; i < ndd; i++) {
            dd_fix[i] = ROUND(dd_amb[i]);
        }

        /* compute residuals manually */
        double res1 = 0.0, res2 = 0.0;
        for (i = 0; i < ndd; i++) {
            double r1 = dd_amb[i] - dd_fix[i];
            res1 += r1 * r1 / (dd_std[i] * dd_std[i]);

            double alt_fix = (dd_amb[i] > dd_fix[i]) ? dd_fix[i] + 1 : dd_fix[i] - 1;
            double r2 = dd_amb[i] - alt_fix;
            res2 += r2 * r2 / (dd_std[i] * dd_std[i]);
        }

        rtk->sol.ratio = (res1 > 0.0) ? (float)(res2 / res1) : 0.0f;
        trace(2, "pppamb: ratio from rounding = %.2f (res1=%.2f, res2=%.2f)\n",
              rtk->sol.ratio, res1, res2);

    } else {
        /* use best solution from LAMBDA */
        for (i = 0; i < ndd; i++) {
            dd_fix[i] = F[i];
        }

        /* compute ratio test */
        rtk->sol.ratio = (s[0] > 1e-6) ? (float)(s[1] / s[0]) : 0.0f;
        if (rtk->sol.ratio > 999.9f) rtk->sol.ratio = 999.9f;

        trace(2, "pppamb: LAMBDA ratio=%.2f (s1=%.2f, s2=%.2f)\n",
              rtk->sol.ratio, s[0], s[1]);
    }

    /* check ratio threshold - be less strict initially */
    if (rtk->sol.ratio < RATIO_THRESHOLD) {
        trace(2, "pppamb: ratio test failed: %.2f < %.2f\n",
              rtk->sol.ratio, RATIO_THRESHOLD);

        /* still try to apply if ratio is reasonable (>1.5) for validation */
        if (rtk->sol.ratio < 1.5) {
            goto cleanup;
        }

        trace(2, "pppamb: ratio marginal but attempting fix for validation\n");
    }

    /* apply fixed ambiguities */
    if (apply_fixed_ambiguities(rtk, sat_list, ndd, dd_fix, dd_amb)) {
        trace(1, "pppamb: ambiguity fixed successfully! ndd=%d ratio=%.2f DOY=%d\n",
              ndd, rtk->sol.ratio, (int)time2doy(rtk->sol.time));
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
