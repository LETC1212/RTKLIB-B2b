/*------------------------------------------------------------------------------
* ppp_ar_passbypass.c : PPP ambiguity resolution using pass-by-pass method
*
* reference :
*    Pass-by-Pass Ambiguity Resolution in Single GPS Receiver PPP Using
*    Observations for Two Sequential Days: An Exploratory Study
*
* version : $Revision:$ $Date:$
* history : 2025/01/xx  1.0  original
*           2025/03/xx  2.0  paper-style one-shot NEQ re-solve
*           2025/03/xx  3.0  arc IF/WL estimation rework
*             [CHANGED] IF:  arc value = LAST epoch's float PPP state
*             [CHANGED] WL:  arc value = inverse-variance weighted HMW mean
*             [ADDED]   pbp_var_wl_epoch() via varerr() error propagation
*             [REMOVED] IF inverse-variance weighted mean accumulation
*             [REMOVED] WL simple mean + SIGMA_MW_CYC^2/nobs
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include "ppp_ar_passbypass.h"

/* ── Pass-2 epoch diagnostics ─────────────────────────────────────────────── */
int pbp_epoch_fix_count = 0;
int pbp_epoch_total     = 0;
int pbp_constraint_sum  = 0;

static void pbp_atexit_summary(void)
{
    if (pbp_epoch_total == 0) return;
    fprintf(stderr,
        "\n======== PBP Pass-2 Epoch Fix Summary (atexit) ========\n"
        "  Total epochs processed   : %d\n"
        "  Epochs with AR fix       : %d  (%.1f%%)\n"
        "  Epochs without AR fix    : %d  (%.1f%%)\n"
        "  Total constraints applied: %d\n"
        "  Avg constraints/epoch    : %.2f\n"
        "  Avg constraints/fix_epoch: %.2f\n",
        pbp_epoch_total,
        pbp_epoch_fix_count,
        100.0*pbp_epoch_fix_count/pbp_epoch_total,
        pbp_epoch_total-pbp_epoch_fix_count,
        100.0*(pbp_epoch_total-pbp_epoch_fix_count)/pbp_epoch_total,
        pbp_constraint_sum,
        (double)pbp_constraint_sum/pbp_epoch_total,
        pbp_epoch_fix_count > 0
            ? (double)pbp_constraint_sum/pbp_epoch_fix_count : 0.0);
    if (pbp_epoch_fix_count == 0)
        fprintf(stderr, "  *** WARNING: 0 fix epochs. Check pbp_ar_debug.log.\n");
    fprintf(stderr, "=======================================================\n\n");
    fflush(stderr);
}
extern void pbp_print_epoch_stats(void) { pbp_atexit_summary(); }

static FILE *pbp_dbg_fp = NULL;
static void pbp_dbg_open(void)
{
    if (!pbp_dbg_fp) {
        pbp_dbg_fp = fopen("pbp_ar_debug.log", "w");
        if (pbp_dbg_fp)
            fprintf(pbp_dbg_fp,
                "# PBP AR diagnostic log\n"
                "# Columns: epoch_time  n_active  nv_raw  nv_after_gate  "
                "applied  nv_constraints  res_vals...\n");
    }
}

/* ── Constants ──────────────────────────────────────────────────────────── */
#ifndef PBP_EPOCH_GAP_SEC
#define PBP_EPOCH_GAP_SEC 60.0
#endif

static int pbp_has_slip(const obsd_t *obs);   /* forward decl */

#define CLIGHT       299792458.0
#define FREQ_GPS_L1  1.57542E9
#define FREQ_GPS_L2  1.22760E9
#ifndef FREQ_BDS_B1I
#define FREQ_BDS_B1I 1561.098e6
#endif
#ifndef FREQ_BDS_B3I
#define FREQ_BDS_B3I 1268.520e6
#endif

/* Fallback WL variance for epochs where varerr returns 0 */
#define SIGMA_MW_CYC 0.35

/* ── PPP state index helpers (must match ppp.c layout) ───────────────────── */
static inline int pbp_NP(const prcopt_t *opt) { return opt->dynamics ? 9 : 3; }
static inline int pbp_NC(const prcopt_t *opt)
{
#ifdef BDS2BDS3
    return NSYS + 1;
#else
    return NSYS;
#endif
}
static inline int pbp_NT(const prcopt_t *opt)
{
    return (opt->tropopt < TROPOPT_EST) ? 0
         : (opt->tropopt == TROPOPT_EST) ? 1 : 3;
}
static inline int pbp_NI(const prcopt_t *opt)
{
    return (opt->ionoopt == IONOOPT_EST) ? MAXSAT : 0;
}
static inline int pbp_ND(const prcopt_t *opt) { return (opt->nf >= 3) ? 1 : 0; }
static inline int pbp_NR(const prcopt_t *opt)
{
    return pbp_NP(opt)+pbp_NC(opt)+pbp_NT(opt)+pbp_NI(opt)+pbp_ND(opt);
}
static inline int pbp_IB(int sat, int f, const prcopt_t *opt)
{
    return pbp_NR(opt) + MAXSAT * f + (sat - 1);
}

/* ── Dual-frequency pair ─────────────────────────────────────────────────── */
static int get_dual_freq_pair(int sat, double *f1, double *f2)
{
    int sys = satsys(sat, NULL);
    if (!f1 || !f2) return 0;
    if (sys == SYS_GPS) { *f1 = FREQ_GPS_L1; *f2 = FREQ_GPS_L2; return 1; }
    if (sys == SYS_CMP) { *f1 = FREQ_BDS_B1I; *f2 = FREQ_BDS_B3I; return 1; }
    return 0;
}

/* ── HMW wide-lane in cycles (system-aware) ──────────────────────────────── */
static double mw_wl_cycles(const obsd_t *obs)
{
    double f1, f2, lam_wl, mw_m;
    if (!obs) return 0.0;
    if (!get_dual_freq_pair(obs->sat, &f1, &f2)) return 0.0;
	int frq2 = (obs->L[1] == 0.0) ? 2 : 1;
    if (obs->L[0]==0.0||obs->L[frq2]==0.0||obs->P[0]==0.0||obs->P[frq2]==0.0)
        return 0.0;
    lam_wl = CLIGHT / (f1 - f2);
    mw_m = (obs->L[0] - obs->L[frq2]) * CLIGHT / (f1 - f2)
          - (f1 * obs->P[0] + f2 * obs->P[frq2]) / (f1 + f2);
    return mw_m / lam_wl;
}

/* ── pbp_var_wl_epoch : single-epoch HMW variance via varerr() ──────────────
 *
 * [NEW in v3.0]
 *
 * HMW combination (matching mw_wl_cycles above):
 *   mw_cyc = (L1_cyc - L2_cyc) - (f1·P1 + f2·P2) / ((f1+f2)·lam_wl)
 *
 * Rearranging in terms of the four raw observations (L1_m, L2_m, P1, P2):
 *   mw_cyc  =  a·L1_m/lam_wl  +  b·L2_m/lam_wl
 *           −  c·P1/lam_wl    −  d·P2/lam_wl
 * where
 *   a =  f1/(f1−f2)    L1 phase coefficient
 *   b = -f2/(f1−f2)    L2 phase coefficient
 *   c =  f1/(f1+f2)    P1 code coefficient
 *   d =  f2/(f1+f2)    P2 code coefficient
 *
 * Error propagation (L1/L2/P1/P2 assumed independent):
 *   var_wl_ep = (a²·var_L1 + b²·var_L2 + c²·var_P1 + d²·var_P2) / lam_wl²
 *
 * varerr() argument f:
 *   0 → L1 phase  (frq=0, code=0)    snr = snr0
 *   1 → P1 code   (frq=0, code=1)    snr = snr0
 *   2 → L2 phase  (frq=1, code=0)    snr = snr1
 *   3 → P2 code   (frq=1, code=1)    snr = snr1
 *
 * Returns var_wl in cycles², or 0.0 on any invalid input.
 * ────────────────────────────────────────────────────────────────────────── */
static double pbp_var_wl_epoch(int sat, int sys, double el,
                                double snr0, double snr1,
                                const prcopt_t *opt, const obsd_t *obs)
{
    double f1, f2, lam_wl;
    double var_L1, var_L2, var_P1, var_P2;
	int frq2 = (obs->L[1] == 0.0) ? 2 : 1;
    if (!obs || !opt) return 0.0;
    if (el <= 0.0) return 0.0;
    if (obs->L[0]==0.0||obs->L[frq2]==0.0||obs->P[0]==0.0||obs->P[frq2]==0.0)
        return 0.0;
    if (!get_dual_freq_pair(sat, &f1, &f2)) return 0.0;

    lam_wl = CLIGHT / (f1 - f2);

    /* pbp_varerr wraps the static varerr() in ppp.c (declared in header) */
    var_L1 = pbp_varerr(sat, sys, el, snr0, 0, opt, obs);  /* L1 phase */
    var_P1 = pbp_varerr(sat, sys, el, snr0, 1, opt, obs);  /* P1 code  */
    var_L2 = pbp_varerr(sat, sys, el, snr1, 2, opt, obs);  /* L2 phase */
    var_P2 = pbp_varerr(sat, sys, el, snr1, 3, opt, obs);  /* P2 code  */

    if (var_L1<=0.0 || var_L2<=0.0 || var_P1<=0.0 || var_P2<=0.0) return 0.0;

    const double a =  f1 / (f1 - f2);   /* L1 coefficient */
    const double b = -f2 / (f1 - f2);   /* L2 coefficient */
    const double c =  f1 / (f1 + f2);   /* P1 coefficient */
    const double d =  f2 / (f1 + f2);   /* P2 coefficient */

    return (a*a*var_L1 + b*b*var_L2 + c*c*var_P1 + d*d*var_P2)
           / (lam_wl * lam_wl);
}

static int NINT(double x) { return (int)floor(x + 0.5); }


static double fix_decision_prob(double frac_abs, double sigma)
{
    if (sigma <= 0.0) return (frac_abs < 0.2) ? 1.0 : 0.0;
    const double s2 = sqrt(2.0) * sigma;
    double xi = 1.0;
    for (int k = 1; k <= 200; k++) {
        double lo = (k - frac_abs) / s2;
        double hi = (k + frac_abs) / s2;
        double term = erfc(lo) - erfc(hi);
        xi -= term;
    }
    return xi;
}

extern int fix_ambiguity(double N_float, double std, double threshold, int *fixed)
{
    const double FRAC_THRESH = 0.2;    /* paper criterion: |frac| < 0.2 cyc */
    const double XI_THRESH   = 0.999;  /* fixable probability threshold 99.9% */

    const int    N_fix = NINT(N_float);
    const double frac  = fabs(N_float - (double)N_fix);

    *fixed = 0;
    if (frac >= FRAC_THRESH) return N_fix;

    /* compute decision function ξ (Eq. 14) */
    double xi = fix_decision_prob(frac, std);
    if (xi >= XI_THRESH) *fixed = 1;

    return N_fix;
}

/* ── Global ambiguity database ───────────────────────────────────────────── */
extern satamb_t satamb[MAXSAT];
extern int      n_ddamb ;
extern ddamb_t  ddamb[MAXSAT * MAXSAT];
extern int      refsat  ;
int      pbp_base_day_id = -1;

/* ── init_arc_data ───────────────────────────────────────────────────────── */
extern void init_arc_data(void)
{
    for (int i = 0; i < MAXSAT; i++) satamb[i].n = 0;
    pbp_base_day_id = -1;
    pbp_epoch_fix_count = 0;
    pbp_epoch_total     = 0;
    pbp_constraint_sum  = 0;
    static int atexit_registered = 0;
    if (!atexit_registered) { atexit(pbp_atexit_summary); atexit_registered = 1; }
    pbp_dbg_open();
}

/* ── print_arc_summary ───────────────────────────────────────────────────── */
extern void print_arc_summary(void)
{
    int i, j, total_arcs = 0, sats_with_both_days = 0;
    int arcs_day0 = 0, arcs_day1 = 0, arcs_other = 0;
    char satid[8];
    printf("\n========== Ambiguity Arc Summary ==========\n");
    for (i = 0; i < MAXSAT; i++)
        for (j = 0; j < satamb[i].n; j++) {
            if      (satamb[i].arc[j].day == 0) arcs_day0++;
            else if (satamb[i].arc[j].day == 1) arcs_day1++;
            else                                 arcs_other++;
        }
    printf("Arc distribution: day0=%d, day1=%d, other=%d\n",
           arcs_day0, arcs_day1, arcs_other);
    for (i = 0; i < MAXSAT; i++) {
        if (satamb[i].n == 0) continue;
        satno2id(i + 1, satid);
        int has_day0 = 0, has_day1 = 0;
        for (j = 0; j < satamb[i].n; j++) {
            if (satamb[i].arc[j].day==0 && satamb[i].arc[j].nobs>=10) has_day0=1;
            if (satamb[i].arc[j].day==1 && satamb[i].arc[j].nobs>=10) has_day1=1;
        }
        if (has_day0 && has_day1) sats_with_both_days++;
        printf("%s: %d arcs [%s]\n", satid, satamb[i].n,
               (has_day0&&has_day1)?"day0+day1":(has_day0?"day0 only":"day1 only"));
        for (j = 0; j < satamb[i].n; j++)
            printf("  Arc %d: day=%d nobs=%d N_IF=%.3f±%.3f N_WL=%.3f±%.4f\n",
                   j, satamb[i].arc[j].day, satamb[i].arc[j].nobs,
                   satamb[i].arc[j].N_IF,  sqrt(satamb[i].arc[j].var_IF),
                   satamb[i].arc[j].N_WL,  sqrt(satamb[i].arc[j].var_WL));
        total_arcs += satamb[i].n;
    }
    printf("Total arcs: %d\n", total_arcs);
    printf("Satellites with both day0 and day1 (>=10 obs): %d\n", sats_with_both_days);
    printf("==========================================\n\n");
}

/* ── collect_ambiguities ─────────────────────────────────────────────────────
 *
 * v3.0 changes:
 *
 * A. IF ambiguity (N_IF / var_IF)                 [CHANGED]
 *    The arc value is ALWAYS overwritten with the CURRENT EPOCH's float PPP
 *    state.  When the loop exits, the arc stores the LAST EPOCH's estimate,
 *    which is the Kalman filter's best (most propagated) value for this arc.
 *    [REMOVED] Inverse-variance weighted mean accumulation.
 *
 * B. WL ambiguity (N_WL / var_WL)                [CHANGED]
 *    Inverse-variance weighted mean of per-epoch HMW values, accumulated
 *    in ambarc_t.wl_wsum / ambarc_t.wl_wxsum (two new fields).
 *    Single-epoch HMW variance computed by pbp_var_wl_epoch() which calls
 *    varerr() through the pbp_varerr() wrapper in ppp.c.
 *    [REMOVED] Simple running mean + fixed SIGMA_MW_CYC^2/nobs.
 *
 * Arc-break logic (LLI slip or time gap > PBP_EPOCH_GAP_SEC) unchanged.
 *-----------------------------------------------------------------------------*/
extern int collect_ambiguities(const rtk_t *rtk, const obsd_t *obs, int n,
                                int day, satamb_t *satamb)
{
    const double ARC_GAP = PBP_EPOCH_GAP_SEC;
    int count = 0;
    if (!rtk || !obs || n <= 0 || !satamb) return 0;

    for (int i = 0; i < n; i++) {
        int sat = obs[i].sat;
        gtime_t time = obs[i].time;
        if (sat <= 0 || sat > MAXSAT) continue;

        int idx_IF = pbp_IB(sat, 0, &rtk->opt);
        if (idx_IF < 0 || idx_IF >= rtk->nx) continue;

        double N_IF   = rtk->x[idx_IF];
        double var_IF = rtk->P[idx_IF + idx_IF * rtk->nx];
        if (N_IF == 0.0) continue;

        double N_WL = mw_wl_cycles(&obs[i]);
        if (N_WL == 0.0) continue;

        /* ── arc break detection ─────────────────────────────────────────── */
        const int has_slip = pbp_has_slip(&obs[i]);
        int idx_arc = -1;
        int is_new  = 0;

        /* fast path: try last arc */
        if (satamb[sat-1].n > 0) {
            int last = satamb[sat-1].n - 1;
            ambarc_t *a = &satamb[sat-1].arc[last];
            if (!has_slip && a->day == day) {
                double dt = timediff(time, a->te);
                if (dt >= 0.0 && dt <= ARC_GAP) idx_arc = last;
            }
        }
        /* full search */
        if (idx_arc < 0) {
            for (int j = 0; j < satamb[sat-1].n; j++) {
                ambarc_t *a = &satamb[sat-1].arc[j];
                if (a->day != day || has_slip) continue;
                double dt = timediff(time, a->te);
                if (dt >= 0.0 && dt <= ARC_GAP) { idx_arc = j; break; }
            }
        }

        /* create new arc */
        if (idx_arc < 0) {
            if (satamb[sat-1].n >= MAXARC) continue;
            idx_arc = satamb[sat-1].n++;
            ambarc_t *a = &satamb[sat-1].arc[idx_arc];
            memset(a, 0, sizeof(*a));
            a->sat = sat; a->day = day;
            a->ts  = time; a->te = time;
            a->nobs = 0;
            /* [CHANGED] WL accumulators start at zero */
            a->wl_wsum  = 0.0;
            a->wl_wxsum = 0.0;
            /* fallback values overwritten immediately below */
            a->N_WL   = N_WL;
            a->var_WL = SIGMA_MW_CYC * SIGMA_MW_CYC;
            a->fixed_WL = 0;
            a->fixed_NL = 0;
            is_new = 1;
        }

        ambarc_t *a = &satamb[sat-1].arc[idx_arc];
        a->te = time;
        a->nobs++;

        /* ── A. IF: always use LAST epoch's float PPP state ─────────────────
         * [CHANGED from v2: was inverse-variance weighted mean]
         * Overwrites the arc value every epoch; at arc end this holds the
         * final Kalman filter estimate for this arc's IF ambiguity.
         * ─────────────────────────────────────────────────────────────────── */
        a->N_IF = N_IF;
        if (var_IF > 0.0) a->var_IF = var_IF;
        /* if var_IF <= 0 (degenerate): keep previous valid variance */

        /* ── B. WL: inverse-variance weighted mean of HMW epochs ────────────
         * [CHANGED from v2: was simple mean + SIGMA_MW_CYC^2/nobs]
         * pbp_var_wl_epoch() propagates the varerr() noise model through the
         * HMW combination to obtain a physically-based per-epoch WL variance.
         * ─────────────────────────────────────────────────────────────────── */
        {
            double el   = rtk->ssat[sat-1].azel[1];
            double snr0 = SNR_UNIT * (double)rtk->ssat[sat-1].snr_rover[0];
            double snr1 = SNR_UNIT * (double)rtk->ssat[sat-1].snr_rover[1];
            int sys_sat = satsys(sat, NULL);

            double vwl = pbp_var_wl_epoch(sat, sys_sat, el, snr0, snr1,
                                           &rtk->opt, &obs[i]);
            if (vwl > 0.0) {
                double w = 1.0 / vwl;
                a->wl_wsum  += w;
                a->wl_wxsum += w * N_WL;
                a->N_WL   = a->wl_wxsum / a->wl_wsum;
                a->var_WL = 1.0 / a->wl_wsum;
            } else if (is_new) {
                /* first epoch, varerr invalid: use fallback simple assignment */
                a->N_WL   = N_WL;
                a->var_WL = SIGMA_MW_CYC * SIGMA_MW_CYC;
            }
		else {
		a->nobs--;
		}
        }

        count++;
    }
    return count;
}

/* ── Helpers ─────────────────────────────────────────────────────────────── */
static double pbp_intdist(double x) { return fabs(x - floor(x + 0.5)); }

static int pbp_has_slip(const obsd_t *obs)
{
    if (!obs) return 0;
    for (int k = 0; k < NFREQ; k++)
        if (obs->LLI[k] & 1) return 1;
    return 0;
}

#ifndef PBP_SIDEREAL_SHIFT_SEC
#define PBP_SIDEREAL_SHIFT_SEC 85920.0
#endif

static double pbp_default_shift_sec(int sat)
{
    (void)sat;
    return PBP_SIDEREAL_SHIFT_SEC;
}

static double pbp_overlap_sec(gtime_t a0, gtime_t a1, gtime_t b0, gtime_t b1)
{
    double A0 = 0.0, A1 = timediff(a1, a0);
    double B0 = -timediff(a0, b0), B1 = -timediff(a0, b1);
    return fmax(0.0, fmin(A1, B1) - fmax(A0, B0));
}

/* ── select_best_arc_pair ─────────────────────────────────────────────────── */
static int select_best_arc_pair(const satamb_t *satamb, int sat,
                                 int *idx0, int *idx1)
{
    const satamb_t *sa;
    double best_overlap = -1.0;
    int best0 = -1, best1 = -1, best_minobs = -1;

    if (!satamb || sat <= 0 || sat > MAXSAT) return 0;
    sa = &satamb[sat-1];
    if (sa->n <= 0) return 0;

    const double shift = pbp_default_shift_sec(sat);
    for (int i0 = 0; i0 < sa->n; i0++) {
        const ambarc_t *a0 = &sa->arc[i0];
        if (a0->day != 0 || a0->nobs < 10) continue;
        gtime_t ts0s = timeadd(a0->ts, shift);
        gtime_t te0s = timeadd(a0->te, shift);
        for (int i1 = 0; i1 < sa->n; i1++) {
            const ambarc_t *a1 = &sa->arc[i1];
            if (a1->day != 1 || a1->nobs < 10) continue;
            double ov = pbp_overlap_sec(ts0s, te0s, a1->ts, a1->te);
            if (ov < 3.0 * 30.0) continue;
            int minobs = a0->nobs < a1->nobs ? a0->nobs : a1->nobs;
            if (ov > best_overlap + 1E-6 ||
                (fabs(ov - best_overlap) <= 1E-6 && minobs > best_minobs)) {
                best_overlap = ov; best0 = i0; best1 = i1; best_minobs = minobs;
            }
        }
    }
    if (best0 >= 0 && best1 >= 0) {
        if (idx0) *idx0 = best0; if (idx1) *idx1 = best1; return 1;
    }

    double best_frac = 1e9; best0 = best1 = -1; best_minobs = -1;
    for (int i0 = 0; i0 < sa->n; i0++) {
        const ambarc_t *a0 = &sa->arc[i0];
        if (a0->day != 0 || a0->nobs < 10) continue;
        for (int i1 = 0; i1 < sa->n; i1++) {
            const ambarc_t *a1 = &sa->arc[i1];
            if (a1->day != 1 || a1->nobs < 10) continue;
            double frac  = pbp_intdist(a1->N_WL - a0->N_WL);
            int minobs   = a0->nobs < a1->nobs ? a0->nobs : a1->nobs;
            if (frac > 0.25 || minobs < 20) continue;
            if (frac < best_frac - 1e-12 ||
                (fabs(frac - best_frac) <= 1e-12 && minobs > best_minobs)) {
                best_frac = frac; best0 = i0; best1 = i1; best_minobs = minobs;
            }
        }
    }
    if (best0 < 0 || best1 < 0) return 0;
    if (idx0) *idx0 = best0; if (idx1) *idx1 = best1;
    return 1;
}

/* ── compute_sd_across_days ─────────────────────────────────────────────── */
static int compute_sd_across_days(const satamb_t *satamb, int sat,
                                   double *SD_WL, double *var_SD_WL,
                                   double *SD_IF,  double *var_SD_IF,
                                   int *i0_out, int *i1_out)
{
    int i0 = -1, i1 = -1;
    if (!SD_WL || !var_SD_WL || !SD_IF || !var_SD_IF) return 0;
    if (!select_best_arc_pair(satamb, sat, &i0, &i1)) return 0;
    if (i0_out) *i0_out = i0; if (i1_out) *i1_out = i1;

    const ambarc_t *a0 = &satamb[sat-1].arc[i0];
    const ambarc_t *a1 = &satamb[sat-1].arc[i1];

    *SD_WL     = a1->N_WL   - a0->N_WL;
    *var_SD_WL = a0->var_WL + a1->var_WL;
    *SD_IF     = a1->N_IF   - a0->N_IF;

    /* arc-level var_IF floor (3.2 mm)^2 per arc */
   // const double var_arc_fl = 1e-5;
   // double vIF0 = (a0->var_IF > 0 && a0->var_IF < var_arc_fl) ? var_arc_fl : a0->var_IF;
   // double vIF1 = (a1->var_IF > 0 && a1->var_IF < var_arc_fl) ? var_arc_fl : a1->var_IF;
   // *var_SD_IF = vIF0 + vIF1;
   *var_SD_IF = a0->var_IF + a1->var_IF;
    return 1;
}

/* ── compute_dd_ambiguities ─────────────────────────────────────────────── */
extern int compute_dd_ambiguities(const satamb_t *satamb, int refsat,
                                   ddamb_t *ddamb, int *n_dd)
{
    int count = 0;
    if (!satamb || !ddamb || !n_dd || refsat <= 0 || refsat > MAXSAT) return 0;
    *n_dd = 0;
    int sys_ref = satsys(refsat, NULL);
    if (sys_ref == 0) return 0;

    double SD_WL_ref=0, varSD_WL_ref=0, SD_IF_ref=0, varSD_IF_ref=0;
    int ref_arc0=-1, ref_arc1=-1;
    if (!compute_sd_across_days(satamb, refsat,
                                 &SD_WL_ref, &varSD_WL_ref,
                                 &SD_IF_ref,  &varSD_IF_ref,
                                 &ref_arc0,   &ref_arc1)) return 0;

    for (int sat = 1; sat <= MAXSAT; sat++) {
        if (sat == refsat || satamb[sat-1].n <= 0) continue;
        if (satsys(sat, NULL) != sys_ref) continue;
        double SD_WL_sat=0, varSD_WL_sat=0, SD_IF_sat=0, varSD_IF_sat=0;
        int sat_arc0=-1, sat_arc1=-1;
        if (!compute_sd_across_days(satamb, sat,
                                     &SD_WL_sat, &varSD_WL_sat,
                                     &SD_IF_sat,  &varSD_IF_sat,
                                     &sat_arc0,   &sat_arc1)) continue;

        ddamb[count].sat1      = refsat;
        ddamb[count].sat2      = sat;
        ddamb[count].DD_WL     = SD_WL_sat - SD_WL_ref;
        ddamb[count].var_DD_WL = varSD_WL_sat + varSD_WL_ref;
        ddamb[count].fixed_WL  = 0;
        ddamb[count].fixed_NL  = 0;
        ddamb[count].DD_IF     = SD_IF_sat - SD_IF_ref;
        ddamb[count].var_DD_IF = varSD_IF_sat + varSD_IF_ref;
        ddamb[count].arc1 = ref_arc0;
        ddamb[count].arc2 = sat_arc0;
        count++;
    }
    *n_dd = count;
    return (count > 0) ? 1 : 0;
}

/* ── fix_wl_nl_ambiguities ───────────────────────────────────────────────── */
extern int fix_wl_nl_ambiguities(ddamb_t *ddamb, int n_dd)
{
    int n_fixed = 0;
    const double th_wl = 0.20, th_nl = 0.20;

    for (int i = 0; i < n_dd; i++) {
        int fixed = 0;
        double std_wl = sqrt(ddamb[i].var_DD_WL);
        ddamb[i].DD_WL_fix = (double)fix_ambiguity(ddamb[i].DD_WL, std_wl, th_wl, &fixed);
        ddamb[i].fixed_WL  = fixed;
        if (!ddamb[i].fixed_WL) continue;

        double f1=0, f2=0;
        if (!get_dual_freq_pair(ddamb[i].sat1, &f1, &f2)) continue;
        const double f1_2=f1*f1, f2_2=f2*f2, den=f1_2-f2_2;
        if (fabs(den) < 1e-6) continue;
        const double lam1=CLIGHT/f1, lam2=CLIGHT/f2;
        const double alpha=f1_2/den, beta=f2_2/den;
        const double C = alpha*lam1 - beta*lam2;
        if (fabs(C) < 1e-12) continue;

        ddamb[i].DD_NL = (ddamb[i].DD_IF - beta*lam2*ddamb[i].DD_WL_fix) / C;
        const double k_if = 1.0/C, k_wl = (-beta*lam2)/C;
        const double var_nl = k_if*k_if*ddamb[i].var_DD_IF
                            + k_wl*k_wl*ddamb[i].var_DD_WL;
        ddamb[i].var_DD_NL = var_nl;

        double std_nl = (var_nl > 0.0) ? sqrt(var_nl) : 1.0;
        ddamb[i].DD_NL_fix = (double)fix_ambiguity(ddamb[i].DD_NL, std_nl, th_nl, &fixed);
        ddamb[i].fixed_NL  = fixed;
        if (!ddamb[i].fixed_NL) continue;

        ddamb[i].DD_IF_fix = C*ddamb[i].DD_NL_fix + beta*lam2*ddamb[i].DD_WL_fix;
        n_fixed++;
    }
    return n_fixed;
}

/* ── Paper-style ambiguity NEQ re-solve (Eq.18-20) ──────────────────────── */

typedef struct {
    int sat, arc, day;
    gtime_t ts, te;
    double b_float, var_float, b_fix, var_fix;
    int used;
} pbp_arcfix_t;

static pbp_arcfix_t *pbp_arcfix   = NULL;
static int           pbp_n_arcfix = 0;
static double        pbp_pb_weight = 1.0e10;
int                  pbp_resolve_flag = 0;

static void pbp_free_arcfix(void)
{
    if (pbp_arcfix) free(pbp_arcfix);
    pbp_arcfix = NULL; pbp_n_arcfix = 0;
}

static int pbp_find_arcfix(int sat, int arc)
{
    for (int i = 0; i < pbp_n_arcfix; i++)
        if (pbp_arcfix[i].used && pbp_arcfix[i].sat==sat && pbp_arcfix[i].arc==arc)
            return i;
    return -1;
}

static int pbp_build_arcfix_list(void)
{
    int n = 0, k = 0;
    pbp_free_arcfix();
    for (int s = 1; s <= MAXSAT; s++)
        for (int j = 0; j < satamb[s-1].n; j++) {
            const ambarc_t *a = &satamb[s-1].arc[j];
            if ((a->day!=0&&a->day!=1)||a->nobs<10||a->var_IF<=0.0||a->N_IF==0.0)
                continue;
            n++;
        }
    if (n <= 0) return 0;
    pbp_arcfix = (pbp_arcfix_t *)calloc((size_t)n, sizeof(pbp_arcfix_t));
    if (!pbp_arcfix) return 0;
    for (int s = 1; s <= MAXSAT; s++)
        for (int j = 0; j < satamb[s-1].n; j++) {
            const ambarc_t *a = &satamb[s-1].arc[j];
            if ((a->day!=0&&a->day!=1)||a->nobs<10||a->var_IF<=0.0||a->N_IF==0.0)
                continue;
            pbp_arcfix[k].sat = s; pbp_arcfix[k].arc = j; pbp_arcfix[k].day = a->day;
            pbp_arcfix[k].ts = a->ts; pbp_arcfix[k].te = a->te;
            pbp_arcfix[k].b_float = a->N_IF; pbp_arcfix[k].var_float = a->var_IF;
            pbp_arcfix[k].b_fix   = a->N_IF; pbp_arcfix[k].var_fix   = a->var_IF;
            pbp_arcfix[k].used = 1; k++;
        }
    pbp_n_arcfix = k;
    return pbp_n_arcfix;
}

extern void pbp_clear_fixed_constraints(void) { pbp_pb_weight=1e10; pbp_free_arcfix(); }
extern int  pbp_has_fixed_constraints(void) { return pbp_n_arcfix > 0; }

extern int pbp_get_fixed_arc_bias(gtime_t t, int sat, double *bias, double *var)
{
    if (!pbp_resolve_flag || pbp_n_arcfix<=0 || !bias) return 0;
    for (int i = 0; i < pbp_n_arcfix; i++) {
        if (!pbp_arcfix[i].used || pbp_arcfix[i].sat!=sat) continue;
        if (timediff(t, pbp_arcfix[i].ts) < -DTTOL) continue;
        if (timediff(t, pbp_arcfix[i].te) >  DTTOL) continue;
        *bias = pbp_arcfix[i].b_fix;
        if (var) *var = pbp_arcfix[i].var_fix;
        return 1;
    }
    return 0;
}

extern int pbp_store_fixed_constraints(const ddamb_t *dd, int n_dd, double Pb)
{
    double *N=NULL, *w=NULL, *x=NULL, *Q=NULL;
    pbp_clear_fixed_constraints();
    if (!dd || n_dd<=0) return 0;
    if (Pb > 0.0) pbp_pb_weight = Pb;
    if (pbp_build_arcfix_list() <= 0) return 0;

    N = zeros(pbp_n_arcfix, pbp_n_arcfix);
    w = zeros(pbp_n_arcfix, 1);
    x = zeros(pbp_n_arcfix, 1);
    if (!N || !w || !x) { free(N); free(w); free(x); pbp_clear_fixed_constraints(); return 0; }

    for (int i = 0; i < pbp_n_arcfix; i++) {
        double var = pbp_arcfix[i].var_float;
        if (var < 1e-8) var = 1e-8;
        N[i + i*pbp_n_arcfix] += 1.0/var;
        w[i] += pbp_arcfix[i].b_float / var;
    }

    for (int i = 0; i < n_dd; i++) {
        int r0, r1, s0, s1, ir0, ir1, is0, is1;
        double bc, d[4] = {+1.0, -1.0, -1.0, +1.0};
        int idx[4];
        if (!dd[i].fixed_WL || !dd[i].fixed_NL) continue;
        if (!select_best_arc_pair(satamb, dd[i].sat1, &r0, &r1)) continue;
        if (!select_best_arc_pair(satamb, dd[i].sat2, &s0, &s1)) continue;
        ir0=pbp_find_arcfix(dd[i].sat1,r0); ir1=pbp_find_arcfix(dd[i].sat1,r1);
        is0=pbp_find_arcfix(dd[i].sat2,s0); is1=pbp_find_arcfix(dd[i].sat2,s1);
        if (ir0<0||ir1<0||is0<0||is1<0) continue;
        bc = dd[i].DD_IF_fix;
        idx[0]=ir0; idx[1]=ir1; idx[2]=is0; idx[3]=is1;
        for (int a=0; a<4; a++) {
            w[idx[a]] += pbp_pb_weight * d[a] * bc;
            for (int b2=0; b2<4; b2++)
                N[idx[a]+idx[b2]*pbp_n_arcfix] += pbp_pb_weight * d[a] * d[b2];
        }
    }

    Q = mat(pbp_n_arcfix, pbp_n_arcfix);
    matcpy(Q, N, pbp_n_arcfix, pbp_n_arcfix);
    if (matinv(Q, pbp_n_arcfix)) {
        free(N); free(w); free(x); free(Q); pbp_clear_fixed_constraints(); return 0;
    }
    matmul("NN", pbp_n_arcfix, 1, pbp_n_arcfix, Q, w, x);
    for (int i = 0; i < pbp_n_arcfix; i++) {
        pbp_arcfix[i].b_fix   = x[i];
        pbp_arcfix[i].var_fix = Q[i + i*pbp_n_arcfix];
        if (pbp_arcfix[i].var_fix < 0.01*0.01) pbp_arcfix[i].var_fix = 0.01*0.01;
        if (pbp_arcfix[i].var_fix > pbp_arcfix[i].var_float)
            pbp_arcfix[i].var_fix = pbp_arcfix[i].var_float;
    }
    free(N); free(w); free(x); free(Q);
    return pbp_n_arcfix;
}

/* ── Legacy link stubs ───────────────────────────────────────────────────── */
extern int pbp_apply_session_pseudoobs(rtk_t *rtk) { (void)rtk; return 0; }
extern int apply_ar_fixed(rtk_t *rtk, const ddamb_t *ddamb, int n_dd)
{
    (void)rtk; (void)ddamb; (void)n_dd; return 0;
}
