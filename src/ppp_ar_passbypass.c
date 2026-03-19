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

#ifndef PBP_MAX2
#define PBP_MAX2(a,b) (( (a) > (b) ) ? (a) : (b))
#endif

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

            /* ── A. IF: arc-level inverse-variance weighted fusion ───────────────
         * Using only the last epoch can keep arc variance too large and make
         * cross-day DD NL fixing over-conservative. Here we fuse epoch IF
         * estimates by 1/var weights to get a stable arc-level IF and variance.
         */
        if (var_IF > 0.0) {
            if (a->nobs <= 1 || a->var_IF <= 0.0) {
                a->N_IF   = N_IF;
                a->var_IF = var_IF;
            }
            else {
                const double w_old = 1.0 / a->var_IF;
                const double w_new = 1.0 / var_IF;
                a->N_IF   = (a->N_IF * w_old + N_IF * w_new) / (w_old + w_new);
                a->var_IF = 1.0 / (w_old + w_new);
            }
        }
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
    if(!pbp_bds_is_sidereal(refsat)&&satsys(refsat,NULL)==SYS_CMP)  return 0;

    double SD_WL_ref=0, varSD_WL_ref=0, SD_IF_ref=0, varSD_IF_ref=0;
    int ref_arc0=-1, ref_arc1=-1;
    if (!compute_sd_across_days(satamb, refsat,
                                 &SD_WL_ref, &varSD_WL_ref,
                                 &SD_IF_ref,  &varSD_IF_ref,
                                 &ref_arc0,   &ref_arc1)) return 0;

    for (int sat = 1; sat <= MAXSAT; sat++) {
        if (sat == refsat || satamb[sat-1].n <= 0) continue;
        if (satsys(sat, NULL) != sys_ref) continue;
        if(!pbp_bds_is_sidereal(sat)&&satsys(sat,NULL)==SYS_CMP)  continue;
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
        if (ddamb[i].var_DD_WL <= 0.0 || ddamb[i].var_DD_IF <= 0.0) continue;
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
        const double var_nl = k_if*k_if*ddamb[i].var_DD_IF;
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
int                  pbp_neq_accum_flag = 0;
int                  pbp_current_day = -1;
gtime_t              pbp_day_start_win[2]={{0}};
gtime_t              pbp_day_end_win[2]={{0}};
int                  pbp_epoch_offset[2]={0};
int                  pbp_ztd_offset[2]={0};
int                  pbp_day_epoch_n[2]={0};
int                  pbp_day_ztd_n[2]={0};


/* [FIX] Reduced from 512 to 200: arc columns are only 3% of n_total (epoch clocks dominate).
 * Combined with n_clk_sys=1: n_total=6014, N matrix=289 MB. Both 512 and 200 work with
 * n_clk_sys=1; the critical fix is n_clk_sys=1 below, not this value. */
#define PBP_MAX_ARC_PARAM 400
#define PBP_MAX_DD_CONSTR 256
#define PBP_MAX_CLKSYS    8

typedef struct {
    int sat, day, arc_id, amb_col;
    gtime_t ts, te;
} pbp_arc_col_t;

typedef struct {
    int sat1, sat2;
    double bc, weight;
} pbp_ddcon_t;

typedef struct {
    int n_xyz, n_amb, n_clk_sys, n_clk_epoch, n_ztd, n_total;
    gtime_t t0, t1; double ti;
    int n_epoch;
    double *N, *w, *xhat;
    gtime_t *epoch_time;
    int *clock_col; /* [n_epoch * n_clk_sys] */
    int *ztd_col;   /* [n_ztd * ntrop]        */
    double *fixed_clk; /* solved clocks for each epoch/system */
    double *day1_fixed_clock; /* solved day1 GPS clock series */
    /* clk_ref[e]: EKF linearisation point for the GPS clock at epoch e (metres).
     * The NEQ solves for CORRECTIONS dx to this reference, so the absolute
     * fixed clock is: x_clk_fixed[e] = clk_ref[e] + xhat[clock_col[e]]       */
    double *clk_ref;
    int day1_epoch_start, day1_epoch_end, day1_epoch_count;
    pbp_arc_col_t arc_cols[PBP_MAX_ARC_PARAM];
    int n_arc_used;
    pbp_ddcon_t ddc[PBP_MAX_DD_CONSTR];
    int n_ddc;
    int ready;
} pbp_neq_t;

static pbp_neq_t g_pbp_neq = {0};

static void pbp_neq_free(void)
{
    free(g_pbp_neq.N); g_pbp_neq.N=NULL;
    free(g_pbp_neq.w); g_pbp_neq.w=NULL;
    free(g_pbp_neq.xhat); g_pbp_neq.xhat=NULL;
    free(g_pbp_neq.epoch_time); g_pbp_neq.epoch_time=NULL;
    free(g_pbp_neq.clock_col); g_pbp_neq.clock_col=NULL;
    free(g_pbp_neq.ztd_col); g_pbp_neq.ztd_col=NULL;
    free(g_pbp_neq.fixed_clk); g_pbp_neq.fixed_clk=NULL;
    free(g_pbp_neq.day1_fixed_clock); g_pbp_neq.day1_fixed_clock=NULL;
    free(g_pbp_neq.clk_ref); g_pbp_neq.clk_ref=NULL;
    memset(&g_pbp_neq,0,sizeof(g_pbp_neq));
}

static int pbp_neq_ntrop(const prcopt_t *opt)
{
    return pbp_NT(opt)>0 ? 1 : 0; /* residual wet delay only */
}

extern void pbp_clear_fixed_constraints(void)
{
    pbp_pb_weight = 1e10;
    pbp_neq_free();
    pbp_current_day = -1;
    memset(pbp_day_start_win,0,sizeof(pbp_day_start_win));
    memset(pbp_day_end_win,0,sizeof(pbp_day_end_win));
    memset(pbp_epoch_offset,0,sizeof(pbp_epoch_offset));
    memset(pbp_ztd_offset,0,sizeof(pbp_ztd_offset));
    memset(pbp_day_epoch_n,0,sizeof(pbp_day_epoch_n));
    memset(pbp_day_ztd_n,0,sizeof(pbp_day_ztd_n));
}
extern int pbp_has_fixed_constraints(void) { return g_pbp_neq.ready; }

extern void pbp_set_day_window(int day, gtime_t ts, gtime_t te, double ti)
{
    if (day < 0 || day > 1) return; if (ti<=0.0){fprintf(stderr,"[PBP-NEQ] day%d ti default 30s\n",day);ti=30.0;}
    pbp_current_day = day;
    pbp_day_start_win[day] = ts;
    pbp_day_end_win[day]   = te;
    pbp_day_epoch_n[day] = (int)floor(timediff(te,ts)/ti + 0.5) + 1;
    pbp_day_ztd_n[day]   = (int)floor(timediff(te,ts)/3600.0 + 1.0) + 1;
    if (day == 0) {
        pbp_epoch_offset[0] = 0;
        pbp_ztd_offset[0] = 0;
    }
    else {
        pbp_epoch_offset[1] = pbp_day_epoch_n[0];
        pbp_ztd_offset[1] = pbp_day_ztd_n[0];
    }
    fprintf(stderr,"[PBP-NEQ] day%d window: %s -> %s, epoch_n=%d ztd_n=%d off_e=%d off_z=%d\n",
            day, time_str(ts,0), time_str(te,0), pbp_day_epoch_n[day], pbp_day_ztd_n[day],
            pbp_epoch_offset[day], pbp_ztd_offset[day]);
}

extern int pbp_neq_init(gtime_t t0, gtime_t t1, double ti, const prcopt_t *opt)
{
    int e, s, ntrop;
    pbp_arc_col_t arc_keep[PBP_MAX_ARC_PARAM];
    pbp_ddcon_t ddc_keep[PBP_MAX_DD_CONSTR];
    int n_arc_keep = g_pbp_neq.n_arc_used;
    int n_ddc_keep = g_pbp_neq.n_ddc;
    if (n_arc_keep>0) memcpy(arc_keep,g_pbp_neq.arc_cols,sizeof(pbp_arc_col_t)*n_arc_keep);
    if (n_ddc_keep>0) memcpy(ddc_keep,g_pbp_neq.ddc,sizeof(pbp_ddcon_t)*n_ddc_keep);
    pbp_neq_free();
    if(!opt){fprintf(stderr,"[PBP-NEQ] opt NULL\n");return 0;} if(timediff(t1,t0)<0.0){fprintf(stderr,"[PBP-NEQ] t1<t0\n");return 0;} if(ti<=0.0){ti=30.0;}
    g_pbp_neq.n_arc_used = n_arc_keep;
    g_pbp_neq.n_ddc = n_ddc_keep;
    if (n_arc_keep>0) memcpy(g_pbp_neq.arc_cols,arc_keep,sizeof(pbp_arc_col_t)*n_arc_keep);
    if (n_ddc_keep>0) memcpy(g_pbp_neq.ddc,ddc_keep,sizeof(pbp_ddcon_t)*n_ddc_keep);
    g_pbp_neq.t0 = t0; g_pbp_neq.t1 = t1; g_pbp_neq.ti = ti;
    g_pbp_neq.n_xyz = 3;
    /* [FIX] Use only 1 receiver clock column per epoch (GPS reference clock).
     * With NSYS=6 and 2-day 30s data: n_total≈35125 → N≈10 GB → zeros() returns NULL.
     * With n_clk_sys=1: n_total≈6326 → N≈320 MB → feasible on cloud server.
     * BDS and other system inter-system biases are white-noise nuisance parameters;
     * keeping only the GPS clock is sufficient to extract position + arc ambiguities.
     * Their H-matrix columns (index NP+1 .. NP+NC-1) remain gmap=-1 and are not
     * accumulated into the NEQ — this is equivalent to pre-eliminating them. */
    g_pbp_neq.n_clk_sys = 1;  /* was: pbp_NC(opt) */
    g_pbp_neq.n_epoch = (int)floor(timediff(t1,t0)/ti + 0.5) + 1;
    g_pbp_neq.n_clk_epoch = g_pbp_neq.n_epoch * g_pbp_neq.n_clk_sys;
    g_pbp_neq.n_ztd = (int)floor(timediff(t1,t0)/3600.0 + 1.0) + 1;
    ntrop = pbp_neq_ntrop(opt);
    g_pbp_neq.epoch_time = (gtime_t*)calloc((size_t)g_pbp_neq.n_epoch,sizeof(gtime_t));
    g_pbp_neq.clock_col = (int*)calloc((size_t)g_pbp_neq.n_clk_epoch,sizeof(int));
    g_pbp_neq.ztd_col   = (int*)calloc((size_t)PBP_MAX2(1,g_pbp_neq.n_ztd*PBP_MAX2(1,ntrop)),sizeof(int));
    g_pbp_neq.fixed_clk = (double*)calloc((size_t)g_pbp_neq.n_clk_epoch,sizeof(double));
    g_pbp_neq.day1_fixed_clock = (double*)calloc((size_t)g_pbp_neq.n_epoch,sizeof(double));
    /* clk_ref: EKF linearisation point for GPS clock at each epoch (metres).
     * Saved in pbp_neq_add_epoch so we can recover absolute clock after NEQ solve. */
    g_pbp_neq.clk_ref = (double*)calloc((size_t)g_pbp_neq.n_clk_epoch,sizeof(double));
    if (!g_pbp_neq.epoch_time || !g_pbp_neq.clock_col || !g_pbp_neq.ztd_col
        || !g_pbp_neq.fixed_clk || !g_pbp_neq.day1_fixed_clock || !g_pbp_neq.clk_ref) {
        fprintf(stderr,"[PBP-NEQ] calloc failed\n"); pbp_neq_free(); return 0;
    }
    for (e=0;e<g_pbp_neq.n_epoch;e++) g_pbp_neq.epoch_time[e]=timeadd(t0,e*ti);
    for (e=0;e<g_pbp_neq.n_epoch;e++) for (s=0;s<g_pbp_neq.n_clk_sys;s++)
        g_pbp_neq.clock_col[e*g_pbp_neq.n_clk_sys+s] = g_pbp_neq.n_xyz + PBP_MAX_ARC_PARAM + e*g_pbp_neq.n_clk_sys + s;
    for (e=0;e<g_pbp_neq.n_ztd*PBP_MAX2(1,ntrop);e++)
        g_pbp_neq.ztd_col[e] = g_pbp_neq.n_xyz + PBP_MAX_ARC_PARAM + g_pbp_neq.n_clk_epoch + e;
    g_pbp_neq.n_total = g_pbp_neq.n_xyz + PBP_MAX_ARC_PARAM + g_pbp_neq.n_clk_epoch + g_pbp_neq.n_ztd*PBP_MAX2(1,ntrop);
    g_pbp_neq.N = zeros(g_pbp_neq.n_total,g_pbp_neq.n_total);
    g_pbp_neq.w = zeros(g_pbp_neq.n_total,1);
    g_pbp_neq.xhat = zeros(g_pbp_neq.n_total,1);
    if(!g_pbp_neq.N||!g_pbp_neq.w||!g_pbp_neq.xhat){fprintf(stderr,"[PBP-NEQ] zeros failed\n");pbp_neq_free();return 0;}
    fprintf(stderr,"[PBP-NEQ] init: epochs=%d clk_sys=%d ztd_blk=%d total=%d\n",
            g_pbp_neq.n_epoch,g_pbp_neq.n_clk_sys,g_pbp_neq.n_ztd,g_pbp_neq.n_total);
    return 1;
}

static int pbp_find_arc_col(int sat, int day, int arc_id)
{
    for (int i=0;i<g_pbp_neq.n_arc_used;i++) {
        if (g_pbp_neq.arc_cols[i].sat==sat && g_pbp_neq.arc_cols[i].day==day && g_pbp_neq.arc_cols[i].arc_id==arc_id)
            return g_pbp_neq.arc_cols[i].amb_col;
    }
    return -1;
}

extern int pbp_build_arc_columns(void)
{
    int new_k=g_pbp_neq.n_arc_used, matched=0, added=0;
    for (int sat=1;sat<=MAXSAT;sat++) {
        for (int j=0;j<satamb[sat-1].n;j++) {
            const ambarc_t *a=&satamb[sat-1].arc[j];
            if ((a->day!=0&&a->day!=1)||a->nobs<10||a->var_IF<=0.0||a->N_IF==0.0) continue;
            int best_k=-1; double best_ov=-1.0;
            for (int k=0;k<g_pbp_neq.n_arc_used;k++) {
                pbp_arc_col_t *ac=&g_pbp_neq.arc_cols[k];
                if (ac->sat!=sat||ac->day!=a->day||ac->arc_id>=0) continue;
                double ov=pbp_overlap_sec(a->ts,a->te,ac->ts,ac->te);
                if (ov>best_ov){best_ov=ov;best_k=k;}
            }
            if (best_k>=0&&best_ov>0.0) {
                g_pbp_neq.arc_cols[best_k].arc_id=j;
                g_pbp_neq.arc_cols[best_k].ts=a->ts;
                g_pbp_neq.arc_cols[best_k].te=a->te;
                matched++;
            } else {
                if (new_k>=PBP_MAX_ARC_PARAM){fprintf(stderr,"[PBP-NEQ] arc overflow\n");continue;}
                g_pbp_neq.arc_cols[new_k].sat=sat; g_pbp_neq.arc_cols[new_k].day=a->day;
                g_pbp_neq.arc_cols[new_k].arc_id=j;
                g_pbp_neq.arc_cols[new_k].ts=a->ts; g_pbp_neq.arc_cols[new_k].te=a->te;
                g_pbp_neq.arc_cols[new_k].amb_col=g_pbp_neq.n_xyz+new_k;
                new_k++; added++;
            }
        }
    }
    g_pbp_neq.n_arc_used=new_k;
    fprintf(stderr,"[PBP-NEQ] arc_cols: %d matched %d new total=%d\n",matched,added,new_k);
    return (new_k>0)?1:0;
}

static int pbp_epoch_id(gtime_t t)
{
    if (pbp_current_day >= 0 && pbp_current_day <= 1 &&
        (pbp_day_start_win[pbp_current_day].time || pbp_day_start_win[pbp_current_day].sec!=0.0)) {
        int e = (int)floor(timediff(t,pbp_day_start_win[pbp_current_day])/g_pbp_neq.ti + 0.5) + pbp_epoch_offset[pbp_current_day];
        if (e<0) e=0;
        if (e>=g_pbp_neq.n_epoch) e=g_pbp_neq.n_epoch-1;
        return e;
    }
    {
        int e = (int)floor(timediff(t,g_pbp_neq.t0)/g_pbp_neq.ti + 0.5);
        if (e<0 || e>=g_pbp_neq.n_epoch) return -1;
        return e;
    }
}

static int pbp_ztd_id(gtime_t t)
{
    if (pbp_current_day >= 0 && pbp_current_day <= 1 &&
        (pbp_day_start_win[pbp_current_day].time || pbp_day_start_win[pbp_current_day].sec!=0.0)) {
        int h = (int)floor(timediff(t,pbp_day_start_win[pbp_current_day])/3600.0) + pbp_ztd_offset[pbp_current_day];
        if (h<0) h=0;
        if (h>=g_pbp_neq.n_ztd) h=g_pbp_neq.n_ztd-1;
        return h;
    }
    {
        int h = (int)floor(timediff(t,g_pbp_neq.t0)/3600.0);
        if (h<0) h=0;
        if (h>=g_pbp_neq.n_ztd) h=g_pbp_neq.n_ztd-1;
        return h;
    }
}

static int pbp_find_arc_by_time(int sat, gtime_t t)
{
    for (int i=0;i<g_pbp_neq.n_arc_used;i++) {
        if (g_pbp_neq.arc_cols[i].sat != sat) continue;
        if (timediff(t,g_pbp_neq.arc_cols[i].ts) < -DTTOL) continue;
        if (timediff(t,g_pbp_neq.arc_cols[i].te) >  DTTOL) continue;
        return g_pbp_neq.arc_cols[i].amb_col;
    }
    return -1;
}

/* pbp_lazy_get_arc_col: find or lazily create arc column during Pass1.
 * No rtk->x[ib]==0 skip: skipping zeroed-out sats zeros N[clk,arc]
 * making DD fixing invisible in the fixed clock series. */
static int pbp_lazy_get_arc_col(int sat, gtime_t t, int has_slip)
{
    for (int k=0;k<g_pbp_neq.n_arc_used;k++) {
        pbp_arc_col_t *ac=&g_pbp_neq.arc_cols[k];
        if (ac->sat!=sat||ac->day!=pbp_current_day) continue;
        double dt=timediff(t,ac->te);
        if (dt<0.0||dt>PBP_EPOCH_GAP_SEC||has_slip) continue;
        ac->te=t; return ac->amb_col;
    }
    if (g_pbp_neq.n_arc_used>=PBP_MAX_ARC_PARAM) return -1;
    int k=g_pbp_neq.n_arc_used;
    g_pbp_neq.arc_cols[k].sat=sat; g_pbp_neq.arc_cols[k].day=pbp_current_day;
    g_pbp_neq.arc_cols[k].arc_id=-1;
    g_pbp_neq.arc_cols[k].ts=t; g_pbp_neq.arc_cols[k].te=t;
    g_pbp_neq.arc_cols[k].amb_col=g_pbp_neq.n_xyz+k;
    g_pbp_neq.n_arc_used++;
    return g_pbp_neq.arc_cols[k].amb_col;
}

extern int pbp_neq_add_epoch(rtk_t *rtk, const obsd_t *obs, int n,
                              const double *v, const double *H,
                              const double *R, int nv)
{
    const prcopt_t *opt;
    double *Ri=NULL;
    int *gmap=NULL;
    int e, ntrop;
    if (!rtk||!obs||n<=0||!v||!H||!R||nv<=0||!pbp_neq_accum_flag) return 0;
    if (!g_pbp_neq.N) return 0;
    opt=&rtk->opt;
    e=pbp_epoch_id(obs[0].time); if (e<0) return 0;
    ntrop=pbp_neq_ntrop(opt);
    Ri=mat(nv,nv); gmap=(int*)malloc(sizeof(int)*rtk->nx);
    if (!Ri||!gmap){free(Ri);free(gmap);return 0;}
    matcpy(Ri,R,nv,nv);
    if (matinv(Ri,nv)){free(Ri);free(gmap);return 0;}

    /* --- gmap: rtk->x index → global NEQ column ----------------------- */
    for (int i=0;i<rtk->nx;i++) gmap[i]=-1;
    gmap[0]=0; gmap[1]=1; gmap[2]=2;                      /* XYZ */
    int clk_col_e = g_pbp_neq.clock_col[e*g_pbp_neq.n_clk_sys+0];
    gmap[pbp_NP(opt)+0] = clk_col_e;                      /* GPS clock */
    if (pbp_NT(opt)>0) {
        int hb=pbp_ztd_id(obs[0].time);
        if (hb>=0) gmap[pbp_NP(opt)+pbp_NC(opt)] =
            g_pbp_neq.ztd_col[hb*PBP_MAX2(1,ntrop)];
    }
    /* Arc ambiguities – lazy build, no x[ib]==0 skip */
    for (int i=0;i<n&&i<MAXOBS;i++) {
        int sat=obs[i].sat, ib=pbp_IB(sat,0,opt);
        if (ib<0||ib>=rtk->nx) continue;
        int col=pbp_lazy_get_arc_col(sat,obs[i].time,pbp_has_slip(&obs[i]));
        if (col>=0) gmap[ib]=col;
    }

    /* --- Save GPS clock linearisation point (metres) for this epoch --- */
    /* The NEQ solves for CORRECTIONS dx to rtk->x (the EKF state).
     * To recover the absolute fixed clock we need:
     *   x_clk_abs[e] = clk_ref[e] + xhat[clk_col_e]
     * clk_ref is set on the FIRST call for epoch e (later calls are filtered
     * by the EKF anyway; overwriting with latest value is also fine since
     * the EKF GPS clock state evolves smoothly). */
    if (g_pbp_neq.clk_ref && clk_col_e>=0)
        g_pbp_neq.clk_ref[e] = rtk->x[pbp_NP(opt)+0];

    /* --- Accumulate H^T R^{-1} H and H^T R^{-1} v -------------------- */
    for (int a=0;a<rtk->nx;a++) {
        int ga=gmap[a]; if (ga<0) continue;
        double wa=0.0;
        for (int k=0;k<nv;k++) {
            double rv=0.0;
            for (int l=0;l<nv;l++) rv+=Ri[k+l*nv]*v[l];
            wa+=H[a+k*rtk->nx]*rv;
        }
        g_pbp_neq.w[ga]+=wa;
        for (int b=0;b<rtk->nx;b++) {
            int gb=gmap[b]; if (gb<0) continue;
            double Nab=0.0;
            for (int k=0;k<nv;k++) {
                double rh=0.0;
                for (int l=0;l<nv;l++) rh+=Ri[k+l*nv]*H[b+l*rtk->nx];
                Nab+=H[a+k*rtk->nx]*rh;
            }
            g_pbp_neq.N[ga+gb*g_pbp_neq.n_total]+=Nab;
        }
    }
    free(Ri); free(gmap);
    return 1;
}

/* pbp_add_one_dd_constraint: time-based arc lookup (safe for arc_id=-1) */
static int pbp_add_one_dd_constraint(int sat1, int sat2, double bc, double wt)
{
    int r0=-1,r1=-1,s0=-1,s1=-1;
    int c[4]; double d[4]={+1.0,-1.0,-1.0,+1.0};
    if(!select_best_arc_pair(satamb,sat1,&r0,&r1))return 0;
    if(!select_best_arc_pair(satamb,sat2,&s0,&s1))return 0;
    c[0]=pbp_find_arc_by_time(sat1,satamb[sat1-1].arc[r0].ts);
    c[1]=pbp_find_arc_by_time(sat1,satamb[sat1-1].arc[r1].ts);
    c[2]=pbp_find_arc_by_time(sat2,satamb[sat2-1].arc[s0].ts);
    c[3]=pbp_find_arc_by_time(sat2,satamb[sat2-1].arc[s1].ts);
    if(c[0]<0){gtime_t m=timeadd(satamb[sat1-1].arc[r0].ts,timediff(satamb[sat1-1].arc[r0].te,satamb[sat1-1].arc[r0].ts)*0.5);c[0]=pbp_find_arc_by_time(sat1,m);}
    if(c[1]<0){gtime_t m=timeadd(satamb[sat1-1].arc[r1].ts,timediff(satamb[sat1-1].arc[r1].te,satamb[sat1-1].arc[r1].ts)*0.5);c[1]=pbp_find_arc_by_time(sat1,m);}
    if(c[2]<0){gtime_t m=timeadd(satamb[sat2-1].arc[s0].ts,timediff(satamb[sat2-1].arc[s0].te,satamb[sat2-1].arc[s0].ts)*0.5);c[2]=pbp_find_arc_by_time(sat2,m);}
    if(c[3]<0){gtime_t m=timeadd(satamb[sat2-1].arc[s1].ts,timediff(satamb[sat2-1].arc[s1].te,satamb[sat2-1].arc[s1].ts)*0.5);c[3]=pbp_find_arc_by_time(sat2,m);}
    if(c[0]<0||c[1]<0||c[2]<0||c[3]<0)return 0;
    for(int a=0;a<4;a++){g_pbp_neq.w[c[a]]+=wt*d[a]*bc;for(int b=0;b<4;b++)g_pbp_neq.N[c[a]+c[b]*g_pbp_neq.n_total]+=wt*d[a]*d[b];}
    return 1;
}

extern int pbp_store_fixed_constraints(const ddamb_t *dd, int n_dd, double Pb)
{
    pbp_pb_weight = Pb>0.0 ? Pb : 1.0e10;
    if (!dd || n_dd<=0) return 0;
    if (!pbp_build_arc_columns()) return 0;
    g_pbp_neq.n_ddc = 0;
    for (int i=0;i<n_dd && g_pbp_neq.n_ddc<PBP_MAX_DD_CONSTR;i++) {
        if (!dd[i].fixed_WL || !dd[i].fixed_NL) continue;
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].sat1 = dd[i].sat1;
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].sat2 = dd[i].sat2;
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].bc   = dd[i].DD_IF_fix;
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].weight = pbp_pb_weight;
        g_pbp_neq.n_ddc++;
    }
    fprintf(stderr,"[PBP-NEQ] selected DD constraints=%d\n",g_pbp_neq.n_ddc);
    return g_pbp_neq.n_ddc;
}

/*===========================================================================
 * pbp_finalize_final_neq  –  Schur complement + Tikhonov solver
 *
 * Layout of compressed NEQ (ncomp columns):
 *   [xyz(3) | arc_used | clock_epoch | ztd]
 *   "slow"   = xyz + arc_used + ztd   (n_s ≈ 230)
 *   "clock"  = clock_epoch            (n_c ≈ 5760, white-noise, diagonal N_cc)
 *
 * Steps:
 *   1. Inject DD pseudo-obs into N, w
 *   2. First compression (remove unused arc slots)
 *   3. Schur: eliminate clock columns analytically (O(n_s² × n_c))
 *   4. Tikhonov: Ns += λ·I  (guarantees invertibility)
 *   5. Solve n_s×n_s system
 *   6. Back-substitute per-epoch clocks
 *   7. Absolute clock = clk_ref[e] + dx_clk[e]
 *      (clk_ref saved in pbp_neq_add_epoch = EKF state at that epoch)
 *=========================================================================*/
extern int pbp_finalize_final_neq(void)
{
    int nclkztd_start = g_pbp_neq.n_xyz + PBP_MAX_ARC_PARAM;
    if (!g_pbp_neq.N || !g_pbp_neq.w) {
        fprintf(stderr,"[PBP-NEQ] ERROR: NEQ not initialized\n"); return 0;
    }

    /* 1. DD constraints */
    int dd_ok=0;
    for (int i=0;i<g_pbp_neq.n_ddc;i++)
        if (pbp_add_one_dd_constraint(g_pbp_neq.ddc[i].sat1,g_pbp_neq.ddc[i].sat2,
                                      g_pbp_neq.ddc[i].bc,g_pbp_neq.ddc[i].weight)) dd_ok++;
    fprintf(stderr,"[PBP-NEQ] DD constraints: %d/%d OK\n",dd_ok,g_pbp_neq.n_ddc);

    /* 2. First compression */
    int *map=(int*)malloc(sizeof(int)*g_pbp_neq.n_total);
    int *invmap=(int*)malloc(sizeof(int)*g_pbp_neq.n_total);
    if (!map||!invmap){free(map);free(invmap);return 0;}
    for (int i=0;i<g_pbp_neq.n_total;i++) invmap[i]=-1;
    int ncomp=0;
    for (int i=0;i<g_pbp_neq.n_xyz;i++) {map[ncomp]=i;invmap[i]=ncomp;ncomp++;}
    for (int i=0;i<g_pbp_neq.n_arc_used;i++) {
        int oldc=g_pbp_neq.arc_cols[i].amb_col;
        if(oldc<0||oldc>=g_pbp_neq.n_total) continue;
        map[ncomp]=oldc; invmap[oldc]=ncomp;
        g_pbp_neq.arc_cols[i].amb_col=ncomp; ncomp++;
    }
    for (int oldc=nclkztd_start;oldc<g_pbp_neq.n_total;oldc++)
        {map[ncomp]=oldc;invmap[oldc]=ncomp;ncomp++;}

    double *Nc=zeros(ncomp,ncomp), *wc=zeros(ncomp,1);
    if (!Nc||!wc){
        fprintf(stderr,"[PBP-NEQ] zeros(Nc) failed (%.0f MB)\n",(double)ncomp*ncomp*8.0/1e6);
        free(map);free(invmap);free(Nc);free(wc);return 0;
    }
    for (int i=0;i<ncomp;i++){
        wc[i]=g_pbp_neq.w[map[i]];
        for (int j=0;j<ncomp;j++) Nc[i+j*ncomp]=g_pbp_neq.N[map[i]+map[j]*g_pbp_neq.n_total];
    }
    /* Remap clock_col/ztd_col before freeing original */
    for (int e=0;e<g_pbp_neq.n_clk_epoch;e++)
        g_pbp_neq.clock_col[e]=invmap[g_pbp_neq.clock_col[e]];
    {
        int nt=pbp_neq_ntrop(&(prcopt_t){0});
        for (int e=0;e<g_pbp_neq.n_ztd*PBP_MAX2(1,nt);e++)
            if (g_pbp_neq.ztd_col&&g_pbp_neq.ztd_col[e]>=0&&g_pbp_neq.ztd_col[e]<g_pbp_neq.n_total)
                g_pbp_neq.ztd_col[e]=invmap[g_pbp_neq.ztd_col[e]];
    }
    free(g_pbp_neq.N);g_pbp_neq.N=NULL;
    free(g_pbp_neq.w);g_pbp_neq.w=NULL;
    free(map);free(invmap);map=NULL;invmap=NULL;
    fprintf(stderr,"[PBP-NEQ] compression: %d->%d (%.0f MB)\n",
            g_pbp_neq.n_total,ncomp,(double)ncomp*ncomp*8.0/1e6);

    /* 3. Schur complement: eliminate per-epoch clock columns */
    int clk0=g_pbp_neq.n_xyz+g_pbp_neq.n_arc_used;
    int clk1=clk0+g_pbp_neq.n_clk_epoch;
    int n_s=ncomp-g_pbp_neq.n_clk_epoch;
    int *slow_map=(int*)malloc(sizeof(int)*ncomp);
    if (!slow_map){free(Nc);free(wc);return 0;}
    {int si=0; for(int i=0;i<ncomp;i++) slow_map[i]=(i>=clk0&&i<clk1)?-1:si++;}
    double *Ns=zeros(n_s,n_s), *ws=zeros(n_s,1);
    if (!Ns||!ws){free(Nc);free(wc);free(slow_map);free(Ns);free(ws);return 0;}
    for (int i=0;i<ncomp;i++){int si=slow_map[i];if(si<0)continue; ws[si]=wc[i];
        for(int j=0;j<ncomp;j++){int sj=slow_map[j];if(sj<0)continue;Ns[si+sj*n_s]=Nc[i+j*ncomp];}}
    int n_elim=0,n_skip=0;
    for (int e=0;e<g_pbp_neq.n_clk_epoch;e++){
        int c=clk0+e; double Ncc=Nc[c+c*ncomp];
        if(fabs(Ncc)<1e-20){n_skip++;continue;}
        double wce=wc[c],inv_Ncc=1.0/Ncc;
        for(int i=0;i<ncomp;i++){int si=slow_map[i];if(si<0)continue;
            double Nsc_i=Nc[i+c*ncomp];if(fabs(Nsc_i)<1e-30)continue;
            ws[si]-=Nsc_i*wce*inv_Ncc;
            for(int j=0;j<ncomp;j++){int sj=slow_map[j];if(sj<0)continue;
                double Nsc_j=Nc[j+c*ncomp];if(fabs(Nsc_j)<1e-30)continue;
                Ns[si+sj*n_s]-=Nsc_i*Nsc_j*inv_Ncc;}}
        n_elim++;
    }
    fprintf(stderr,"[PBP-NEQ] Schur: n_s=%d elim=%d skip=%d\n",n_s,n_elim,n_skip);

    /* 4. Direct solve of the reduced n_s×n_s system (no Tikhonov)
     * After fixing OMC-based w construction and compression to only
     * actual arc parameters, N_s should be full-rank and directly invertible.
     * (Paper Eq.6: [x;b;u] = Q * [wx;wb;wu]) */
    fprintf(stderr,"[PBP-NEQ] Direct solve: n_s=%d (no Tikhonov regularization)\n",n_s);

    /* 5. Solve n_s×n_s system: Q_s = N_s^{-1}, then x_s = Q_s * w_s */
    double *Qs=mat(n_s,n_s), *xs=zeros(n_s,1);
    if(!Qs||!xs){free(Nc);free(wc);free(slow_map);free(Ns);free(ws);free(Qs);free(xs);return 0;}

    /* Pre-inversion rank check: verify no zero diagonal entries */
    {
        int n_zero_diag=0;
        for(int i=0;i<n_s;i++){
            if(fabs(Ns[i+i*n_s])<1e-30) n_zero_diag++;
        }
        if(n_zero_diag>0){
            fprintf(stderr,"[PBP-NEQ] WARNING: %d/%d zero diagonal entries in Ns before inversion\n",
                    n_zero_diag,n_s);
            for(int i=0;i<n_s&&i<30;i++)
                fprintf(stderr,"  Ns_diag[%d]=%.6e\n",i,Ns[i+i*n_s]);
        }
    }

    matcpy(Qs,Ns,n_s,n_s); free(Ns);Ns=NULL;
    if(matinv(Qs,n_s)){
        fprintf(stderr,"[PBP-NEQ] ERROR: Ns singular (n_s=%d). "
                "Check that OMC (prefit) residuals are used and arc columns are correctly built.\n",n_s);
        free(Nc);free(wc);free(slow_map);free(ws);free(Qs);free(xs);return 0;}
    matmul("NN",n_s,1,n_s,Qs,ws,xs);
    free(Qs);Qs=NULL; free(ws);ws=NULL;

    /* 6. Back-substitute: dx_clk[e] = (w_c_e - N_cs_e^T * x_s) / N_cc_e */
    double *xfull=zeros(ncomp,1);
    if(!xfull){free(Nc);free(wc);free(slow_map);free(xs);return 0;}
    for(int i=0;i<ncomp;i++){int si=slow_map[i];if(si>=0) xfull[i]=xs[si];}
    free(xs);xs=NULL;
    int n_rec=0;
    for(int e=0;e<g_pbp_neq.n_clk_epoch;e++){
        int c=clk0+e; double Ncc=Nc[c+c*ncomp];
        if(fabs(Ncc)<1e-20){xfull[c]=0.0;continue;}
        double rhs=wc[c];
        for(int i=0;i<ncomp;i++){int si=slow_map[i];if(si<0)continue; rhs-=Nc[c+i*ncomp]*xfull[i];}
        xfull[c]=rhs/Ncc; n_rec++;
    }
    free(Nc);Nc=NULL; free(wc);wc=NULL; free(slow_map);slow_map=NULL;

    /* 7. Absolute clock = clk_ref[e] + dx_clk[e]
     *    clk_ref[e] = rtk->x[NP] at epoch e (EKF linearisation point, metres)
     *    xfull[clk_col] = dx_clk[e]  (correction, metres)
     *    Together they give the batch-LS absolute clock in metres.
     *    The .stat file outputs dtr[0]*1e9 ns where dtr[0]=rtk->x[NP]/CLIGHT,
     *    so our absolute value is directly comparable after *M2NS conversion.  */
    free(g_pbp_neq.xhat); g_pbp_neq.xhat=xfull; xfull=NULL;

    for(int e=0;e<g_pbp_neq.n_epoch;e++){
        for(int s=0;s<g_pbp_neq.n_clk_sys;s++){
            int col=g_pbp_neq.clock_col[e*g_pbp_neq.n_clk_sys+s];
            double dx = (col>=0&&col<ncomp) ? g_pbp_neq.xhat[col] : 0.0;
            double ref = (g_pbp_neq.clk_ref && e<g_pbp_neq.n_epoch) ?
                         g_pbp_neq.clk_ref[e*g_pbp_neq.n_clk_sys+s] : 0.0;
            /* absolute clock in metres */
            g_pbp_neq.fixed_clk[e*g_pbp_neq.n_clk_sys+s] = ref + dx;
        }
    }
    g_pbp_neq.day1_epoch_start=pbp_epoch_offset[1];
    g_pbp_neq.day1_epoch_end  =pbp_epoch_offset[1]+pbp_day_epoch_n[1]-1;
    if(g_pbp_neq.day1_epoch_start<0) g_pbp_neq.day1_epoch_start=0;
    if(g_pbp_neq.day1_epoch_end>=g_pbp_neq.n_epoch)
        g_pbp_neq.day1_epoch_end=g_pbp_neq.n_epoch-1;
    g_pbp_neq.day1_epoch_count=
        g_pbp_neq.day1_epoch_end-g_pbp_neq.day1_epoch_start+1;
    for(int e=0;e<g_pbp_neq.day1_epoch_count;e++){
        int ge=g_pbp_neq.day1_epoch_start+e;
        g_pbp_neq.day1_fixed_clock[e]=g_pbp_neq.fixed_clk[ge*g_pbp_neq.n_clk_sys+0];
    }
    g_pbp_neq.ready=1;
    fprintf(stderr,
        "[PBP-NEQ] solve OK: arcs=%d dd_ok=%d ncomp=%d n_s=%d clk_rec=%d day1=%d\n",
        g_pbp_neq.n_arc_used,dd_ok,ncomp,n_s,n_rec,g_pbp_neq.day1_epoch_count);
    printf("[PBP-NEQ] solve OK: arcs=%d dd=%d day1_epochs=%d\n",
           g_pbp_neq.n_arc_used,g_pbp_neq.n_ddc,g_pbp_neq.day1_epoch_count);
    return 1;
}

extern int pbp_write_day1_fixed_clock_file(const char *path)
{
    /* fixed_clk is in metres (c*dt); convert to ns: m * (1e9/CLIGHT) */
    const double M2NS = 1e9 / CLIGHT;
    FILE *fp; int y,m,d,hh,mm; double ss;
    if (!path||!*path||!g_pbp_neq.ready||g_pbp_neq.day1_epoch_count<=0) {
        fprintf(stderr,"[PBP-NEQ] ERROR: invalid write request path=%s ready=%d day1_epochs=%d\n",
                path?path:"(null)",g_pbp_neq.ready,g_pbp_neq.day1_epoch_count);
        printf("[PBP-NEQ] ERROR: invalid write request path=%s ready=%d day1_epochs=%d\n",
               path?path:"(null)",g_pbp_neq.ready,g_pbp_neq.day1_epoch_count);
        return 0;
    }
    createdir(path); fp=fopen(path,"w");
    if (!fp){fprintf(stderr,"[PBP-NEQ] ERROR: fopen failed %s\n",path);return 0;}
    fprintf(fp,"# PBP fixed receiver clock (GPS, day1)  unit: ns\n");
    fprintf(fp,"# YYYY/MM/DD HH:MM:SS  clock_ns\n");
    for (int e=0;e<g_pbp_neq.day1_epoch_count;e++) {
        gtime_t t=g_pbp_neq.epoch_time[g_pbp_neq.day1_epoch_start+e];
        double ep[6]; time2epoch(t,ep);
        y=(int)ep[0];m=(int)ep[1];d=(int)ep[2];
        hh=(int)ep[3];mm=(int)ep[4];ss=ep[5];
        fprintf(fp,"%04d/%02d/%02d %02d:%02d:%02d %.6f\n",
                y,m,d,hh,mm,(int)floor(ss+0.5),
                g_pbp_neq.day1_fixed_clock[e]*M2NS);
    }
    fclose(fp);
    fprintf(stderr,"[PBP-NEQ] wrote fixed clock (ns): %s (epochs=%d)\n",
            path,g_pbp_neq.day1_epoch_count);
    printf("[PBP-NEQ] wrote fixed clock (ns): %s (epochs=%d)\n",
           path,g_pbp_neq.day1_epoch_count);
    return 1;
}

extern int pbp_get_fixed_clock(gtime_t t, int sys_idx, double *clk)
{
    int e;
    if (!g_pbp_neq.ready || !clk) return 0;
    e = pbp_epoch_id(t);
    if (e<0) return 0;
    if (sys_idx<0 || sys_idx>=g_pbp_neq.n_clk_sys) return 0;
    *clk = g_pbp_neq.fixed_clk[e*g_pbp_neq.n_clk_sys+sys_idx];
    return 1;
}

extern int pbp_get_fixed_arc_bias(gtime_t t, int sat, double *bias, double *var)
{
    (void)t; (void)sat; (void)bias; (void)var;
    return 0; /* arc-level IF prior path removed */
}

/* ── Legacy link stubs ───────────────────────────────────────────────────── */
extern int pbp_apply_session_pseudoobs(rtk_t *rtk) { (void)rtk; return 0; }
extern int apply_ar_fixed(rtk_t *rtk, const ddamb_t *ddamb, int n_dd)
{
    (void)rtk; (void)ddamb; (void)n_dd; return 0;
}
extern int pbp_bds_is_sidereal(int sat)
{
    int prn = 0;
    if (satsys(sat, &prn) != SYS_CMP) return 0;
    if (prn >= 1  && prn <= 10) return 1;
    if (prn == 13 || prn == 16) return 1;
    if (prn >= 38 && prn <= 40) return 1;
    if (prn >= 59 && prn <= 62) return 1;
    return 0;
}
