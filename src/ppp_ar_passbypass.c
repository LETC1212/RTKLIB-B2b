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

            /* ── A. IF: last epoch value + conditional variance ─────────────────
         *
         * VALUE: always overwrite with current epoch (last epoch = best EKF
         * estimate incorporating all prior observations).
         *
         * VARIANCE: use the CONDITIONAL variance P[b|clk], NOT the marginal
         * variance P[b,b].
         *
         * Why: The marginal variance P[b,b] is dominated by the clock-ambiguity
         * correlation (receiver clock is ~0.3-1m white noise).  Typical values:
         *   P[b,b] ≈ (5-10mm)²     → DD std_NL ≈ 0.1-0.2 cycles
         *   ξ test fails at 99.9%  → ZERO ambiguities fixed.
         *
         * In the DD computation, receiver clocks cancel (same-station,
         * same-epoch between-satellite difference).  The physically relevant
         * uncertainty is the CONDITIONAL variance with clock removed:
         *   P[b|clk] = P[b,b] - P[b,clk]² / P[clk,clk]
         * Typical values:
         *   P[b|clk] ≈ (0.5-2mm)²  → DD std_NL ≈ 0.01-0.04 cycles
         *   ξ test passes easily for correct integers.
         *
         * The previous inverse-variance weighted mean accidentally achieved
         * small variance (by treating correlated epochs as independent), but
         * the resulting value was physically meaningless and too small,
         * causing wrong integers to be accepted.  The conditional variance
         * is both physically correct AND the right magnitude.
         *
         * Variance floor (3mm)² prevents degenerate cases.
         */
        if (var_IF > 0.0) {
            const double VAR_IF_FLOOR = 9.0e-6; /* (3 mm)^2 */

            /* Determine system-specific clock index for this satellite */
            int sys_s = satsys(sat, NULL), clk_k = 0;
            switch (sys_s) {
                case SYS_GPS: clk_k = 0; break;
                case SYS_GLO: clk_k = 1; break;
                case SYS_GAL: clk_k = 2; break;
                case SYS_CMP: clk_k = 3; break;
                case SYS_IRN: clk_k = 4; break;
                default:      clk_k = 0; break;
            }
#ifdef BDS2BDS3
            if (sys_s == SYS_CMP) {
                int prn_s = 0; satsys(sat, &prn_s);
                if (prn_s > 16) clk_k = NSYS;
            }
#endif
            int ic = pbp_NP(&rtk->opt) + clk_k;

            /* Conditional variance: remove clock contribution via Schur */
            double var_cond = var_IF;
            if (ic < rtk->nx && idx_IF < rtk->nx) {
                double P_bc = rtk->P[idx_IF + (size_t)ic * rtk->nx];
                double P_cc = rtk->P[ic + (size_t)ic * rtk->nx];
                if (P_cc > 1e-20) {
                    double vc = var_IF - P_bc * P_bc / P_cc;
                    if (vc > VAR_IF_FLOOR) var_cond = vc;
                    else                    var_cond = VAR_IF_FLOOR;
                }
            }

            a->N_IF   = N_IF;
            a->var_IF = (var_cond > VAR_IF_FLOOR) ? var_cond : VAR_IF_FLOOR;
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

#ifndef PBP_SIDEREAL_SHIFT_GPS
#define PBP_SIDEREAL_SHIFT_GPS 86164.09  /* GPS repeat ≈ 1 sidereal day (23h56m4s) */
#endif
#ifndef PBP_SIDEREAL_SHIFT_BDS
#define PBP_SIDEREAL_SHIFT_BDS 86170.00  /* BDS GEO/IGSO repeat period */
#endif

static double pbp_default_shift_sec(int sat)
{
    int sys = satsys(sat, NULL);
    if (sys == SYS_CMP) return PBP_SIDEREAL_SHIFT_BDS;
    return PBP_SIDEREAL_SHIFT_GPS;
}

static double pbp_overlap_sec(gtime_t a0, gtime_t a1, gtime_t b0, gtime_t b1)
{
    double A0 = 0.0, A1 = timediff(a1, a0);
    double B0 = -timediff(a0, b0), B1 = -timediff(a0, b1);
    return fmax(0.0, fmin(A1, B1) - fmax(A0, B0));
}

/* ── select_best_arc_pair ─────────────────────────────────────────────────
 * For a given satellite, find the day0 arc and day1 arc that correspond
 * to the same sidereal repeat pass.
 *
 * Method: shift day0 arc times by the sidereal day offset, then find
 * the day1 arc with maximum overlap. This is the ONLY selection criterion.
 * The old "frac-based fallback" (select by WL small fractional part)
 * had no theoretical basis and could match wrong arc pairs.
 * ──────────────────────────────────────────────────────────────────────── */
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
            if (ov < 3.0 * 30.0) continue; /* minimum 90s overlap */
            int minobs = a0->nobs < a1->nobs ? a0->nobs : a1->nobs;
            if (ov > best_overlap + 1E-6 ||
                (fabs(ov - best_overlap) <= 1E-6 && minobs > best_minobs)) {
                best_overlap = ov; best0 = i0; best1 = i1; best_minobs = minobs;
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

/* ══════════════════════════════════════════════════════════════════════════
 * Day1 NEQ v9: ALL parameters in ONE matrix, Cholesky solve
 *
 * Matrix: [xyz(3) | arcs(≤100) | ALL_clk(n_ep×n_aclk) | ztd(≤50)]
 * ALL clock systems (GPS+BDS2+BDS3) are in the global matrix.
 * NO Schur elimination. NO parameter separation. ZERO approximation.
 *
 * Cholesky solve: L*L^T decomposition + forward/back substitution
 *   - O(n^3/6) ≈ 57 seconds for n=8793
 *   - In-place: no extra n×n allocation (saves 619 MB)
 *   - Numerically stable for positive-definite matrices
 * ══════════════════════════════════════════════════════════════════════════ */

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
int                  pbp_epoch_collected = 0;
int                  pbp_current_day = -1;
gtime_t              pbp_day_start_win[2]={{0}};
gtime_t              pbp_day_end_win[2]={{0}};
int                  pbp_epoch_offset[2]={0};
int                  pbp_ztd_offset[2]={0};
int                  pbp_day_epoch_n[2]={0};
int                  pbp_day_ztd_n[2]={0};

#define PBP_MAX_ARC_DAY1   100
#define PBP_MAX_DD_CONSTR  256
#define PBP_MAX_CLK_SYS      8
#define PBP_MAX_ZTD_DAY1    50

typedef struct {
    int sat, day, arc_id, amb_col;
    gtime_t ts, te;
    double xlin_wsum, xlin_wdenom;
} pbp_arc_col_t;

typedef struct {
    int sat1, sat2;
    double bc, weight;
    int col_ref, col_sat;
    double DD_IF_fix;      /* integer-fixed DD IF value */
    double b_ref_d0, b_sat_d0; /* day0 float arc values (known constants) */
} pbp_ddcon_t;

typedef struct {
    int n_xyz, n_aclk, n_epoch, n_ztd, n_trop_per_block, n_total;
    int clk_map[PBP_MAX_CLK_SYS]; /* EKF k → active clock index */
    int n_arc_used;
    pbp_arc_col_t arc_cols[PBP_MAX_ARC_DAY1];
    pbp_ddcon_t ddc[PBP_MAX_DD_CONSTR];
    int n_ddc;

    double *N;        /* [n_total × n_total] */
    double *w;        /* [n_total]           */
    double *clk_lin;  /* [n_epoch * n_aclk]  */
    double *clk_float;/* [n_epoch * n_aclk]  */
    int *epoch_valid;   /* 1 = EKF produced a solution at this epoch */
    int *epoch_has_neq; /* 1 = NEQ data accumulated (converged epoch) */
    gtime_t *epoch_time;

    int day1_epoch_count;
    double *day1_fixed_clock;
    double *day1_neqfloat_clock; /* NEQ float (no DD constraints) */
    int ready;
    gtime_t t0, t1; double ti;
} pbp_neq_t;

static pbp_neq_t g_pbp_neq = {0};

/* ══════════════════════════════════════════════════════════════════════════
 * Diagnostic instrumentation (logging & counters only, NO math changes)
 * ══════════════════════════════════════════════════════════════════════════ */
#define PBP_DIAG_MAX_EPOCH 8000

typedef struct {
    /* Global counters */
    int add_call_total;
    int add_ok_total;
    int fail_bad_input;
    int fail_epoch_id;
    int fail_R_invert;
    int fail_no_active_param;
    int fail_vrms_too_large;   /* hard reject (nwrms >= T2) */
    int fail_alloc;
    int schur_fail_total;        /* always 0 in v9 (no Schur) */
    int epoch_collision_total;
    int backsub_ok_count;
    int backsub_fallback_count;
    /* Normalized-residual gate statistics */
    int gate_soft_weight_count;  /* T1 < nwrms < T2: soft downweighted */
    double nwrms_sum;
    int    nwrms_count;
    double nwrms_max;
    /* Per-epoch arrays */
    int    epoch_hit_count[PBP_DIAG_MAX_EPOCH];
    int    epoch_nv_last[PBP_DIAG_MAX_EPOCH];
    int    epoch_nact_last[PBP_DIAG_MAX_EPOCH];
    int    epoch_backsub_ok[PBP_DIAG_MAX_EPOCH];
} pbp_diag_t;

static pbp_diag_t g_pbp_diag = {0};

static void pbp_diag_reset(void)
{
    memset(&g_pbp_diag, 0, sizeof(g_pbp_diag));
}

/* ── Cholesky solve: L*L^T*x = b, in-place ─────────────────────────── */
static int cholesky_solve(double *A, double *b, int n)
{
    int i, j, k;
    double s;
    /* L*L^T decomposition (lower triangle stored in A) */
    for (j = 0; j < n; j++) {
        s = A[j + j*(size_t)n];
        for (k = 0; k < j; k++) s -= A[j + k*(size_t)n] * A[j + k*(size_t)n];
        if (s <= 1e-30) return -1;
        A[j + j*(size_t)n] = sqrt(s);
        for (i = j+1; i < n; i++) {
            s = A[i + j*(size_t)n];
            for (k = 0; k < j; k++) s -= A[i + k*(size_t)n] * A[j + k*(size_t)n];
            A[i + j*(size_t)n] = s / A[j + j*(size_t)n];
        }
    }
    /* Forward: L*y = b */
    for (i = 0; i < n; i++) {
        s = b[i];
        for (k = 0; k < i; k++) s -= A[i + k*(size_t)n] * b[k];
        b[i] = s / A[i + i*(size_t)n];
    }
    /* Backward: L^T*x = y */
    for (i = n-1; i >= 0; i--) {
        s = b[i];
        for (k = i+1; k < n; k++) s -= A[k + i*(size_t)n] * b[k];
        b[i] = s / A[i + i*(size_t)n];
    }
    return 0;
}

static void pbp_neq_free(void)
{
    free(g_pbp_neq.N); free(g_pbp_neq.w);
    free(g_pbp_neq.clk_lin); free(g_pbp_neq.clk_float);
    free(g_pbp_neq.epoch_valid); free(g_pbp_neq.epoch_has_neq);
    free(g_pbp_neq.epoch_time); free(g_pbp_neq.day1_fixed_clock);
    free(g_pbp_neq.day1_neqfloat_clock);
    memset(&g_pbp_neq,0,sizeof(g_pbp_neq));
}

extern void pbp_clear_fixed_constraints(void)
{
    pbp_pb_weight=1e10; pbp_neq_free(); pbp_current_day=-1;
    memset(pbp_day_start_win,0,sizeof(pbp_day_start_win));
    memset(pbp_day_end_win,0,sizeof(pbp_day_end_win));
    memset(pbp_epoch_offset,0,sizeof(pbp_epoch_offset));
    memset(pbp_ztd_offset,0,sizeof(pbp_ztd_offset));
    memset(pbp_day_epoch_n,0,sizeof(pbp_day_epoch_n));
    memset(pbp_day_ztd_n,0,sizeof(pbp_day_ztd_n));
}
extern int pbp_has_fixed_constraints(void){return g_pbp_neq.ready;}

extern void pbp_set_day_window(int day,gtime_t ts,gtime_t te,double ti)
{
    if(day<0||day>1)return;
    if(ti<=0){
        fprintf(stderr,"[PBP-DIAG] WARNING: ti=%.1f<=0 in set_day_window(day=%d), "
                "replacing with 30.0s\n",ti,day);
        ti=30.0;
    }
    pbp_current_day=day;
    pbp_day_start_win[day]=ts; pbp_day_end_win[day]=te;
    pbp_day_epoch_n[day]=(int)floor(timediff(te,ts)/ti+0.5)+1;
    pbp_day_ztd_n[day]=(int)floor(timediff(te,ts)/3600.0+1.0)+1;
    if(day==0){pbp_epoch_offset[0]=0;pbp_ztd_offset[0]=0;}
    else{pbp_epoch_offset[1]=pbp_day_epoch_n[0];pbp_ztd_offset[1]=pbp_day_ztd_n[0];}
    fprintf(stderr,"[PBP-DIAG] set_day_window: day=%d ti=%.1fs n_epoch=%d "
            "ep_offset=%d ztd_offset=%d\n",
            day,ti,pbp_day_epoch_n[day],pbp_epoch_offset[day],pbp_ztd_offset[day]);
}

extern int pbp_neq_init(gtime_t t0,gtime_t t1,double ti,const prcopt_t *opt)
{
    int nc,k;
    pbp_neq_free();
    pbp_diag_reset();
    if(!opt||timediff(t1,t0)<0)return 0;
    if(ti<=0){
        fprintf(stderr,"[PBP-DIAG] WARNING: ti=%.1f<=0 in pbp_neq_init, "
                "replacing with 30.0s\n",ti);
        ti=30.0;
    }
    g_pbp_neq.t0=t0; g_pbp_neq.t1=t1; g_pbp_neq.ti=ti;
    g_pbp_neq.n_xyz=3;
    g_pbp_neq.n_epoch=(int)floor(timediff(t1,t0)/ti+0.5)+1;
    g_pbp_neq.n_ztd=(int)floor(timediff(t1,t0)/3600.0+1.0)+1;
    if(g_pbp_neq.n_ztd>PBP_MAX_ZTD_DAY1)g_pbp_neq.n_ztd=PBP_MAX_ZTD_DAY1;

    /* Map ALL nc clock systems including BDS3 */
    nc=pbp_NC(opt);
    g_pbp_neq.n_aclk=0;
    for(k=0;k<PBP_MAX_CLK_SYS;k++)g_pbp_neq.clk_map[k]=-1;
    for(k=0;k<nc&&k<PBP_MAX_CLK_SYS;k++){
        int sys;
        switch(k){
            case 0:sys=SYS_GPS;break;case 1:sys=SYS_GLO;break;
            case 2:sys=SYS_GAL;break;case 3:sys=SYS_CMP;break;
            case 4:sys=SYS_IRN;break;default:sys=SYS_CMP;break;
        }
        if(opt->navsys&sys) g_pbp_neq.clk_map[k]=g_pbp_neq.n_aclk++;
    }
    if(g_pbp_neq.n_aclk<=0)return 0;

    /* Troposphere: 1 param (ZTD only) for EST, 3 params (ZTD+GN+GE) for ESTG */
    g_pbp_neq.n_trop_per_block = pbp_NT(opt); /* 0, 1, or 3 */
    if(g_pbp_neq.n_trop_per_block<0) g_pbp_neq.n_trop_per_block=0;

    int ne=g_pbp_neq.n_epoch, na=g_pbp_neq.n_aclk;
    int ntrop_total = g_pbp_neq.n_ztd * g_pbp_neq.n_trop_per_block;
    g_pbp_neq.n_total=3+PBP_MAX_ARC_DAY1+ne*na+ntrop_total;
    int nt=g_pbp_neq.n_total;

    g_pbp_neq.N        =(double*)calloc((size_t)nt*nt,sizeof(double));
    g_pbp_neq.w        =(double*)calloc(nt,sizeof(double));
    g_pbp_neq.clk_lin  =(double*)calloc((size_t)ne*na,sizeof(double));
    g_pbp_neq.clk_float=(double*)calloc((size_t)ne*na,sizeof(double));
    g_pbp_neq.epoch_valid=(int*)calloc(ne,sizeof(int));
    g_pbp_neq.epoch_has_neq=(int*)calloc(ne,sizeof(int));
    g_pbp_neq.epoch_time=(gtime_t*)calloc(ne,sizeof(gtime_t));
    g_pbp_neq.day1_fixed_clock=(double*)calloc(ne,sizeof(double));
    g_pbp_neq.day1_neqfloat_clock=(double*)calloc(ne,sizeof(double));
    if(!g_pbp_neq.N||!g_pbp_neq.w||!g_pbp_neq.clk_lin||!g_pbp_neq.clk_float||
       !g_pbp_neq.epoch_valid||!g_pbp_neq.epoch_has_neq||
       !g_pbp_neq.epoch_time||!g_pbp_neq.day1_fixed_clock||
       !g_pbp_neq.day1_neqfloat_clock){
        pbp_neq_free();return 0;
    }
    for(k=0;k<ne;k++) g_pbp_neq.epoch_time[k]=timeadd(t0,k*ti);
    fprintf(stderr,"[PBP-NEQ] init: ep=%d aclk=%d trop=%dx%d nt=%d (%.0f MB) Cholesky~%.0fs\n",
            ne,na,g_pbp_neq.n_ztd,g_pbp_neq.n_trop_per_block,
            nt,(double)nt*nt*8/1e6,(double)nt*nt*(double)nt/6.0/2e9);
    fprintf(stderr,"[PBP-DIAG] init: ti=%.1fs n_epoch=%d day_offsets=[%d,%d] "
            "ztd_offsets=[%d,%d]\n",
            ti,ne,pbp_epoch_offset[0],pbp_epoch_offset[1],
            pbp_ztd_offset[0],pbp_ztd_offset[1]);
    return 1;
}

/* ── Column helpers ────────────────────────────────────────────────────── */
static int pbp_epoch_id_day1(gtime_t t){
    int e=(int)floor(timediff(t,g_pbp_neq.t0)/g_pbp_neq.ti+0.5);
    if(e<0)e=0;if(e>=g_pbp_neq.n_epoch)e=g_pbp_neq.n_epoch-1;return e;
}
static int pbp_clk_col(int epoch,int ekf_k){
    int a=g_pbp_neq.clk_map[ekf_k]; if(a<0)return -1;
    return g_pbp_neq.n_xyz+PBP_MAX_ARC_DAY1+epoch*g_pbp_neq.n_aclk+a;
}
/* Troposphere column: trop_idx=0→ZTD, 1→GN, 2→GE within each hour block */
static int pbp_trop_col(gtime_t t, int trop_idx){
    int ntpb=g_pbp_neq.n_trop_per_block;
    if(ntpb<=0||trop_idx<0||trop_idx>=ntpb)return -1;
    int h=(int)floor(timediff(t,g_pbp_neq.t0)/3600.0);
    if(h<0)h=0;if(h>=g_pbp_neq.n_ztd)h=g_pbp_neq.n_ztd-1;
    return g_pbp_neq.n_xyz+PBP_MAX_ARC_DAY1+g_pbp_neq.n_epoch*g_pbp_neq.n_aclk
           +h*ntpb+trop_idx;
}

/* ── Arc column ────────────────────────────────────────────────────────── */
static int pbp_get_arc_col_day1(int sat,gtime_t t){
    int arc_id,k;
    if(sat<=0||sat>MAXSAT||satamb[sat-1].n<=0)return -1;
    arc_id=-1;
    for(k=0;k<satamb[sat-1].n;k++){
        const ambarc_t *a=&satamb[sat-1].arc[k];
        if(a->day!=1)continue;
        if(timediff(t,a->ts)<-DTTOL||timediff(t,a->te)>PBP_EPOCH_GAP_SEC+DTTOL)continue;
        arc_id=k;break;
    }
    if(arc_id<0)return -1;
    for(k=0;k<g_pbp_neq.n_arc_used;k++)
        if(g_pbp_neq.arc_cols[k].sat==sat&&g_pbp_neq.arc_cols[k].arc_id==arc_id)
            return g_pbp_neq.arc_cols[k].amb_col;
    if(g_pbp_neq.n_arc_used>=PBP_MAX_ARC_DAY1)return -1;
    k=g_pbp_neq.n_arc_used;
    g_pbp_neq.arc_cols[k].sat=sat;g_pbp_neq.arc_cols[k].day=1;
    g_pbp_neq.arc_cols[k].arc_id=arc_id;
    g_pbp_neq.arc_cols[k].ts=satamb[sat-1].arc[arc_id].ts;
    g_pbp_neq.arc_cols[k].te=satamb[sat-1].arc[arc_id].te;
    g_pbp_neq.arc_cols[k].amb_col=g_pbp_neq.n_xyz+k;
    g_pbp_neq.arc_cols[k].xlin_wsum=0;g_pbp_neq.arc_cols[k].xlin_wdenom=0;
    g_pbp_neq.n_arc_used++;
    return g_pbp_neq.arc_cols[k].amb_col;
}
static int pbp_find_arc_col_day1(int sat,int arc_id){
    int k;
    for(k=0;k<g_pbp_neq.n_arc_used;k++)
        if(g_pbp_neq.arc_cols[k].sat==sat&&g_pbp_neq.arc_cols[k].arc_id==arc_id)
            return g_pbp_neq.arc_cols[k].amb_col;
    return -1;
}
extern int pbp_build_arc_columns(void){
    int sat,j,k,added=0;
    for(sat=1;sat<=MAXSAT;sat++)
        for(j=0;j<satamb[sat-1].n;j++){
            const ambarc_t *a=&satamb[sat-1].arc[j];
            if(a->day!=1||a->nobs<10)continue;
            if(pbp_find_arc_col_day1(sat,j)>=0)continue;
            if(g_pbp_neq.n_arc_used>=PBP_MAX_ARC_DAY1)continue;
            k=g_pbp_neq.n_arc_used;
            g_pbp_neq.arc_cols[k].sat=sat;g_pbp_neq.arc_cols[k].day=1;
            g_pbp_neq.arc_cols[k].arc_id=j;g_pbp_neq.arc_cols[k].ts=a->ts;
            g_pbp_neq.arc_cols[k].te=a->te;g_pbp_neq.arc_cols[k].amb_col=g_pbp_neq.n_xyz+k;
            g_pbp_neq.arc_cols[k].xlin_wsum=0;g_pbp_neq.arc_cols[k].xlin_wdenom=0;
            g_pbp_neq.n_arc_used++;added++;
        }
    fprintf(stderr,"[PBP-NEQ] arcs: %d\n",g_pbp_neq.n_arc_used);
    return g_pbp_neq.n_arc_used>0;
}

static double pbp_htrh(const double *H,const double *Ri,int a,int b,int nx,int nv){
    double s=0;int k,l;
    for(k=0;k<nv;k++){double r=0;for(l=0;l<nv;l++)r+=Ri[k+l*nv]*H[b+l*nx];s+=H[a+k*nx]*r;}
    return s;
}

/*===========================================================================
 * pbp_neq_add_epoch — ALL params in one global matrix, no quality gate.
 *
 * CRITICAL ordering:
 *   1. Compute epoch index, save epoch_time
 *   2. Save clk_lin / clk_float / epoch_valid  (BEFORE gate!)
 *   3. Gate: soft downweight only (no hard reject in stage 1)
 *   4. Build gmap, accumulate N/w
 *=========================================================================*/
extern int pbp_neq_add_epoch(rtk_t *rtk,const obsd_t *obs,int n,
                              const double *v,const double *H,
                              const double *R,int nv,const double *x_lin)
{
    const prcopt_t *opt;
    double *Ri=NULL;
    int *gmap=NULL;
    int e,nc,nt,i,k;

    g_pbp_diag.add_call_total++;

    if(!rtk||!obs||n<=0||!v||!H||!R||nv<=0||!x_lin){
        g_pbp_diag.fail_bad_input++; return 0;
    }
    if(!pbp_neq_accum_flag||!g_pbp_neq.N){
        g_pbp_diag.fail_bad_input++; return 0;
    }
    opt=&rtk->opt; nc=pbp_NC(opt);
    nt=g_pbp_neq.n_total;

    /* ── Step 1: epoch index and epoch_time ─────────────────────────── */
    e=pbp_epoch_id_day1(obs[0].time);
    if(e<0||e>=g_pbp_neq.n_epoch){
        g_pbp_diag.fail_epoch_id++; return 0;
    }
    g_pbp_neq.epoch_time[e]=obs[0].time;

    /* ── Step 2: ALWAYS save clk_lin/clk_float/epoch_valid ──────────
     * These MUST be saved even if the gate later downweights heavily.
     * Without clk_lin[e], the backsub step has no linearization point
     * and falls back to zero → backsub_fallback_ratio ≈ 100%.
     * Without epoch_valid[e]=1, the epoch is excluded from output. */
    {
        int na=g_pbp_neq.n_aclk;
        for(k=0;k<nc&&k<PBP_MAX_CLK_SYS;k++){
            int a=g_pbp_neq.clk_map[k]; if(a<0)continue;
            int ic=pbp_NP(opt)+k;
            if(ic<rtk->nx){
                g_pbp_neq.clk_lin[e*na+a]=x_lin[ic];
                g_pbp_neq.clk_float[e*na+a]=rtk->x[ic];
            }
        }
    }
    g_pbp_neq.epoch_valid[e]=1;

    /* ── Step 3: Convergence check ────────────────────────────────────
     * During early EKF convergence (first ~10min), x_predicted has large
     * position errors (~km) → v_omc = y - h(x_predicted) has RMS >> 50m.
     * After convergence, position is correct; v_omc RMS is ~1-5m (dominated
     * by clock prediction error, which is normal for the NEQ).
     *
     * Unconverged epochs: skip NEQ accumulation, output clk_float directly.
     * This avoids injecting huge v_omc into the law equations while still
     * producing a complete clock series (float values for early epochs,
     * fixed values for converged epochs). */
    {
        double vrms=0;
        for(k=0;k<nv;k++) vrms+=v[k]*v[k];
        vrms=sqrt(vrms/(nv>0?nv:1));
        if(vrms>=50.0){
            /* Unconverged: epoch_valid=1 but epoch_has_neq=0 → output float */
            g_pbp_diag.fail_vrms_too_large++;
            return 0;
        }
    }
    g_pbp_neq.epoch_has_neq[e]=1;

    /* ── Step 4: Invert R and accumulate N/w ─────────────────────────── */
    Ri=mat(nv,nv);if(!Ri){g_pbp_diag.fail_alloc++;return 0;}
    matcpy(Ri,R,nv,nv);
    if(matinv(Ri,nv)){free(Ri);g_pbp_diag.fail_R_invert++;return 0;}

    /* ── Step 4: Build gmap and accumulate N/w ─────────────────────── */
    gmap=(int*)malloc(sizeof(int)*rtk->nx);
    if(!gmap){free(Ri);g_pbp_diag.fail_alloc++;return 0;}
    for(i=0;i<rtk->nx;i++) gmap[i]=-1;

    gmap[0]=0;gmap[1]=1;gmap[2]=2; /* XYZ */

    /* ALL clock systems → global columns */
    for(k=0;k<nc&&k<PBP_MAX_CLK_SYS;k++){
        if(g_pbp_neq.clk_map[k]<0)continue;
        int ic=pbp_NP(opt)+k;
        if(ic>=rtk->nx)continue;
        int has=0;
        for(i=0;i<nv;i++)if(fabs(H[ic+i*rtk->nx])>1e-30){has=1;break;}
        if(has) gmap[ic]=pbp_clk_col(e,k);
    }

    /* Troposphere: map all NT params (ZTD, and GN/GE if ESTG) */
    {
        int nt_trop=pbp_NT(opt);
        for(int t_idx=0;t_idx<nt_trop;t_idx++){
            int ekf_idx=pbp_NP(opt)+pbp_NC(opt)+t_idx;
            if(ekf_idx<rtk->nx){
                int col=pbp_trop_col(obs[0].time,t_idx);
                if(col>=0) gmap[ekf_idx]=col;
            }
        }
    }

    /* Day1 arcs */
    for(i=0;i<n&&i<MAXOBS;i++){
        int sat=obs[i].sat,ib=pbp_IB(sat,0,opt);
        if(ib<0||ib>=rtk->nx)continue;
        int ok=0;
        for(k=0;k<nv;k++)if(fabs(H[ib+k*rtk->nx])>1e-30){ok=1;break;}
        if(!ok)continue;
        int col=pbp_get_arc_col_day1(sat,obs[i].time);
        if(col>=0) gmap[ib]=col;
    }

    /* Collect active params */
    int nact=0;
    int act_x[600],act_g[600];
    for(i=0;i<rtk->nx;i++)
        if(gmap[i]>=0&&nact<600){act_x[nact]=i;act_g[nact]=gmap[i];nact++;}

    if(nact<=0){
        free(Ri);free(gmap);
        g_pbp_diag.fail_no_active_param++;return 0;
    }

    /* Ri*v */
    double *Riv=mat(nv,1);
    if(!Riv){free(Ri);free(gmap);g_pbp_diag.fail_alloc++;return 0;}
    for(k=0;k<nv;k++){double s=0;for(i=0;i<nv;i++)s+=Ri[k+i*nv]*v[i];Riv[k]=s;}

    /* Direct accumulation: N += H^T Ri H,  w += H^T Ri v */
    for(i=0;i<nact;i++){
        int a=act_x[i],gi=act_g[i];
        double wa=0;
        for(k=0;k<nv;k++) wa+=H[a+k*rtk->nx]*Riv[k];
        g_pbp_neq.w[gi]+=wa;
        for(int j=0;j<nact;j++){
            int b=act_x[j],gj=act_g[j];
            g_pbp_neq.N[gi+(size_t)gj*nt]+=pbp_htrh(H,Ri,a,b,rtk->nx,nv);
        }
        /* Track xlin for arcs */
        if(gi>=g_pbp_neq.n_xyz && gi<g_pbp_neq.n_xyz+PBP_MAX_ARC_DAY1){
            int ai=gi-g_pbp_neq.n_xyz;
            double Nab=pbp_htrh(H,Ri,a,a,rtk->nx,nv);
            g_pbp_neq.arc_cols[ai].xlin_wsum+=Nab*x_lin[a];
            g_pbp_neq.arc_cols[ai].xlin_wdenom+=Nab;
        }
    }

    free(Riv);free(Ri);free(gmap);

    /* Diag: record per-epoch stats */
    if(e>=0 && e<PBP_DIAG_MAX_EPOCH){
        g_pbp_diag.epoch_hit_count[e]++;
        if(g_pbp_diag.epoch_hit_count[e]>1)
            g_pbp_diag.epoch_collision_total++;
        g_pbp_diag.epoch_nv_last[e]=nv;
        g_pbp_diag.epoch_nact_last[e]=nact;
    }
    g_pbp_diag.add_ok_total++;
    return 1;
}

/*===========================================================================
 * pbp_store_fixed_constraints
 *=========================================================================*/
extern int pbp_store_fixed_constraints(const ddamb_t *dd,int n_dd,double Pb)
{
    int i,r0,r1,s0,s1,ref,sat,col_ref,col_sat;
    pbp_pb_weight=Pb>0?Pb:1e10;
    if(!dd||n_dd<=0)return 0;
    if(!pbp_build_arc_columns())return 0;
    g_pbp_neq.n_ddc=0;
    for(i=0;i<n_dd&&g_pbp_neq.n_ddc<PBP_MAX_DD_CONSTR;i++){
        if(!dd[i].fixed_WL||!dd[i].fixed_NL)continue;
        ref=dd[i].sat1;sat=dd[i].sat2;
        r0=-1;r1=-1;s0=-1;s1=-1;
        if(!select_best_arc_pair(satamb,ref,&r0,&r1))continue;
        if(!select_best_arc_pair(satamb,sat,&s0,&s1))continue;
        col_ref=pbp_find_arc_col_day1(ref,r1);
        col_sat=pbp_find_arc_col_day1(sat,s1);
        if(col_ref<0||col_sat<0)continue;
        {int ai_r=col_ref-g_pbp_neq.n_xyz,ai_s=col_sat-g_pbp_neq.n_xyz;
         if(ai_r<0||ai_s<0)continue;
         if(g_pbp_neq.arc_cols[ai_r].xlin_wdenom<=0||
            g_pbp_neq.arc_cols[ai_s].xlin_wdenom<=0)continue;}

        /* Store DD info — bc will be computed in finalize() after Pass A,
         * using xhat_float to ensure perfect consistency with the NEQ. */
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].sat1=ref;
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].sat2=sat;
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].bc=0.0; /* placeholder */
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].weight=pbp_pb_weight;
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].col_ref=col_ref;
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].col_sat=col_sat;
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].DD_IF_fix=dd[i].DD_IF_fix;
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].b_ref_d0=satamb[ref-1].arc[r0].N_IF;
        g_pbp_neq.ddc[g_pbp_neq.n_ddc].b_sat_d0=satamb[sat-1].arc[s0].N_IF;
        g_pbp_neq.n_ddc++;
    }
    fprintf(stderr,"[PBP-NEQ] DD stored: %d (bc deferred to finalize)\n",g_pbp_neq.n_ddc);
    return g_pbp_neq.n_ddc;
}

/*===========================================================================
 * pbp_finalize_final_neq — constrain, inject DD, compress, Cholesky solve
 *=========================================================================*/
extern int pbp_finalize_final_neq(void)
{
    int nt=g_pbp_neq.n_total,ne=g_pbp_neq.n_epoch,na=g_pbp_neq.n_aclk,i,j;
    if(!g_pbp_neq.N||!g_pbp_neq.w)return 0;

    /* 0. Constrain XYZ and ZTD — VERY tight to prevent DD→position→clock leak.
     * DD constraints modify ambiguities; through N cross-terms this pulls
     * position/ZTD; through geometry-dependent coupling this creates slow
     * clock variations → excess TDEV at long τ.
     * Making these essentially rigid forces DD effect into ambiguities+clocks
     * only, eliminating low-frequency contamination. */
    {
        double wt_xyz=1.0/(1e-4*1e-4); /* σ=0.1mm — nearly rigid */
        for(i=0;i<3;i++) g_pbp_neq.N[i+(size_t)i*nt]+=wt_xyz;
        /* Trop: ZTD σ=1mm, gradients σ=0.5mm */
        {
            int ntpb=g_pbp_neq.n_trop_per_block;
            int trop_base=g_pbp_neq.n_xyz+PBP_MAX_ARC_DAY1+ne*na;
            for(i=0;i<g_pbp_neq.n_ztd;i++){
                for(int t_idx=0;t_idx<ntpb;t_idx++){
                    int col=trop_base+i*ntpb+t_idx;
                    if(col<nt&&fabs(g_pbp_neq.N[col+(size_t)col*nt])>1e-25){
                        double sigma=(t_idx==0)?0.001:0.0005;
                        g_pbp_neq.N[col+(size_t)col*nt]+=1.0/(sigma*sigma);
                    }
                }
            }
        }
    }

    /* ════════════════════════════════════════════════════════════════════
     * Two-pass solve WITH compression:
     *   Pass A: NEQ float (no DD) — compress + Cholesky
     *   Pass B: NEQ fixed (with DD) — inject DD, compress + Cholesky
     * N/w preserved between passes, freed after Pass B.
     * ════════════════════════════════════════════════════════════════════ */

    double *xhat_float=NULL, *xhat_fixed=NULL;

    /* ── Pass A: NEQ float (no DD) ─────────────────────────────────── */
    {
        int *cmap=(int*)calloc(nt,sizeof(int));
        if(!cmap)return 0;
        for(i=0;i<nt;i++)cmap[i]=-1;
        int ncomp=0;
        for(i=0;i<nt;i++)if(fabs(g_pbp_neq.N[i+(size_t)i*nt])>1e-25)cmap[i]=ncomp++;
        if(ncomp<=0){free(cmap);return 0;}

        double *Nc=(double*)calloc((size_t)ncomp*ncomp,sizeof(double));
        double *wc=(double*)calloc(ncomp,sizeof(double));
        if(!Nc||!wc){free(cmap);free(Nc);free(wc);return 0;}
        for(i=0;i<nt;i++){int ci=cmap[i];if(ci<0)continue;
            wc[ci]=g_pbp_neq.w[i];
            for(j=0;j<nt;j++){int cj=cmap[j];if(cj<0)continue;
                Nc[ci+(size_t)cj*ncomp]=g_pbp_neq.N[i+(size_t)j*nt];}
        }
        fprintf(stderr,"[PBP-NEQ] Pass A (float): compress %d→%d, Cholesky...\n",nt,ncomp);
        if(cholesky_solve(Nc,wc,ncomp)){
            fprintf(stderr,"[PBP-NEQ] Pass A Cholesky failed\n");
            free(cmap);free(Nc);free(wc);return 0;
        }
        free(Nc);
        xhat_float=zeros(nt,1);
        if(!xhat_float){free(cmap);free(wc);return 0;}
        for(i=0;i<nt;i++){int ci=cmap[i];if(ci>=0) xhat_float[i]=wc[ci];}
        free(wc);free(cmap);
    }

    /* ── Compute bc, validate against float solution, then inject DD ──────
     *
     * Correct bc in the NEQ correction space (dx = x_true - x_lin_eff):
     *   dx_ref - dx_sat = -DD_IF_fix + (b_ref_d0 - b_sat_d0)
     *                     + (xlin_sat - xlin_ref)
     *
     * xlin_ref/sat = weighted-mean x_lin for each arc, tracked exactly
     * during NEQ accumulation (xlin_wsum/xlin_wdenom).
     *
     * VALIDATION (new): After Pass A, the float solution xhat_float gives
     * the best batch estimate.  For each DD constraint, check if the float
     * arc values are consistent with the integer constraint:
     *   residual = (dx_ref_float - dx_sat_float) - bc
     * In NL cycles: res_nl = residual / C_nl
     * If |res_nl| > 0.25, the integer is likely wrong → skip this DD.
     *
     * This catches wrong integers that passed the initial ξ test due to
     * variance underestimation or systematic biases. */

    /* NL wavelength C for validation threshold (GPS L1/L2) */
    const double f1_v=FREQ_GPS_L1, f2_v=FREQ_GPS_L2;
    const double C_nl = f1_v*f1_v/(f1_v*f1_v-f2_v*f2_v)*(CLIGHT/f1_v)
                      - f2_v*f2_v/(f1_v*f1_v-f2_v*f2_v)*(CLIGHT/f2_v);
    const double DD_VAL_THRESH_NL = 0.25; /* reject if |res| > 0.25 NL cycles */
    int dd_ok=0, dd_rejected=0;
    for(i=0;i<g_pbp_neq.n_ddc;i++){
        int cr=g_pbp_neq.ddc[i].col_ref, cs=g_pbp_neq.ddc[i].col_sat;
        double wt=g_pbp_neq.ddc[i].weight;
        if(cr<0||cr>=nt||cs<0||cs>=nt)continue;

        /* NEQ linearization points for day1 arcs */
        int ai_r=cr-g_pbp_neq.n_xyz, ai_s=cs-g_pbp_neq.n_xyz;
        double xlin_ref=0,xlin_sat=0;
        if(ai_r>=0&&g_pbp_neq.arc_cols[ai_r].xlin_wdenom>0)
            xlin_ref=g_pbp_neq.arc_cols[ai_r].xlin_wsum/g_pbp_neq.arc_cols[ai_r].xlin_wdenom;
        if(ai_s>=0&&g_pbp_neq.arc_cols[ai_s].xlin_wdenom>0)
            xlin_sat=g_pbp_neq.arc_cols[ai_s].xlin_wsum/g_pbp_neq.arc_cols[ai_s].xlin_wdenom;

        double bc = -g_pbp_neq.ddc[i].DD_IF_fix
                  + (g_pbp_neq.ddc[i].b_ref_d0 - g_pbp_neq.ddc[i].b_sat_d0)
                  + (xlin_sat - xlin_ref);
        g_pbp_neq.ddc[i].bc = bc;

        /* ── DD validation against Pass A float solution ──────────────
         * The float NEQ solution gives the optimal unconstrained arc values.
         * If the DD constraint disagrees with the float solution by more than
         * 0.25 NL cycles, the integer is almost certainly wrong.  Injecting
         * a wrong constraint would introduce a systematic bias (~10 cm in IF)
         * that propagates to clock through NEQ cross-correlations, degrading
         * MDEV at long τ far more than leaving the ambiguity float. */
        {
            double dx_ref_f = (cr>=0&&cr<nt) ? xhat_float[cr] : 0.0;
            double dx_sat_f = (cs>=0&&cs<nt) ? xhat_float[cs] : 0.0;
            double res = (dx_ref_f - dx_sat_f) - bc;
            double res_nl = (fabs(C_nl) > 1e-12) ? fabs(res) / C_nl : 9.9;
            if (res_nl > DD_VAL_THRESH_NL) {
                char sid1[8], sid2[8];
                satno2id(g_pbp_neq.ddc[i].sat1, sid1);
                satno2id(g_pbp_neq.ddc[i].sat2, sid2);
                fprintf(stderr, "[PBP-NEQ] DD REJECTED %s-%s: "
                        "res=%.4f m (%.3f NL cyc) > %.2f threshold\n",
                        sid1, sid2, res, res_nl, DD_VAL_THRESH_NL);
                dd_rejected++;
                continue; /* skip this DD constraint */
            }
        }

        g_pbp_neq.N[cr+(size_t)cr*nt]+=wt;g_pbp_neq.N[cs+(size_t)cs*nt]+=wt;
        g_pbp_neq.N[cr+(size_t)cs*nt]-=wt;g_pbp_neq.N[cs+(size_t)cr*nt]-=wt;
        g_pbp_neq.w[cr]+=wt*bc;g_pbp_neq.w[cs]-=wt*bc;
        dd_ok++;
    }
    fprintf(stderr,"[PBP-NEQ] DD injected: %d/%d rejected: %d (Pb=%.0e, bc from xlin)\n",
            dd_ok,g_pbp_neq.n_ddc,dd_rejected,pbp_pb_weight);

    /* ── Pass B: NEQ fixed (with DD) ───────────────────────────────── */
    {
        int *cmap=(int*)calloc(nt,sizeof(int));
        if(!cmap){free(xhat_float);return 0;}
        for(i=0;i<nt;i++)cmap[i]=-1;
        int ncomp=0;
        for(i=0;i<nt;i++)if(fabs(g_pbp_neq.N[i+(size_t)i*nt])>1e-25)cmap[i]=ncomp++;
        fprintf(stderr,"[PBP-NEQ] Pass B (fixed): compress %d→%d DD=%d, Cholesky...\n",nt,ncomp,dd_ok);
        if(ncomp<=0){free(cmap);free(xhat_float);return 0;}

        double *Nc=(double*)calloc((size_t)ncomp*ncomp,sizeof(double));
        double *wc=(double*)calloc(ncomp,sizeof(double));
        if(!Nc||!wc){free(cmap);free(Nc);free(wc);free(xhat_float);return 0;}
        for(i=0;i<nt;i++){int ci=cmap[i];if(ci<0)continue;
            wc[ci]=g_pbp_neq.w[i];
            for(j=0;j<nt;j++){int cj=cmap[j];if(cj<0)continue;
                Nc[ci+(size_t)cj*ncomp]=g_pbp_neq.N[i+(size_t)j*nt];}
        }
        free(g_pbp_neq.N);g_pbp_neq.N=NULL;
        free(g_pbp_neq.w);g_pbp_neq.w=NULL;

        if(cholesky_solve(Nc,wc,ncomp)){
            fprintf(stderr,"[PBP-NEQ] Pass B Cholesky failed\n");
            free(cmap);free(Nc);free(wc);free(xhat_float);return 0;
        }
        free(Nc);
        xhat_fixed=zeros(nt,1);
        if(!xhat_fixed){free(cmap);free(wc);free(xhat_float);return 0;}
        for(i=0;i<nt;i++){int ci=cmap[i];if(ci>=0) xhat_fixed[i]=wc[ci];}
        free(wc);free(cmap);
    }

    /* ── Extract clocks from both solutions ────────────────────────── */
    int gps_a=g_pbp_neq.clk_map[0];
    g_pbp_neq.day1_epoch_count=ne;
    {
        int n_fixed=0, n_float_fb=0, n_none=0;
        double sum_off_fix=0, sum_off_flt=0;
        int cnt_off_fix=0, cnt_off_flt=0;

        for(i=0;i<ne;i++){
            int col=g_pbp_neq.n_xyz+PBP_MAX_ARC_DAY1+i*na+gps_a;
            double ekf_flt=g_pbp_neq.clk_float[i*na+gps_a];
            double lin=g_pbp_neq.clk_lin[i*na+gps_a];

            if(g_pbp_neq.epoch_has_neq && g_pbp_neq.epoch_has_neq[i]){
                double dx_fix=(col>=0&&col<nt)?xhat_fixed[col]:0.0;
                g_pbp_neq.day1_fixed_clock[i]=lin+dx_fix;
                sum_off_fix+=(lin+dx_fix-ekf_flt); cnt_off_fix++;

                double dx_flt=(col>=0&&col<nt)?xhat_float[col]:0.0;
                g_pbp_neq.day1_neqfloat_clock[i]=lin+dx_flt;
                sum_off_flt+=(lin+dx_flt-ekf_flt); cnt_off_flt++;

                n_fixed++;
                g_pbp_diag.backsub_ok_count++;
                if(i<PBP_DIAG_MAX_EPOCH) g_pbp_diag.epoch_backsub_ok[i]=1;
            } else if(g_pbp_neq.epoch_valid[i]){
                g_pbp_neq.day1_fixed_clock[i]=ekf_flt;
                g_pbp_neq.day1_neqfloat_clock[i]=ekf_flt;
                n_float_fb++;
                g_pbp_diag.backsub_fallback_count++;
                if(i<PBP_DIAG_MAX_EPOCH) g_pbp_diag.epoch_backsub_ok[i]=0;
            } else {
                g_pbp_neq.day1_fixed_clock[i]=0.0;
                g_pbp_neq.day1_neqfloat_clock[i]=0.0;
                n_none++;
            }
        }

        if(cnt_off_fix>0){
            double mean_off=sum_off_fix/cnt_off_fix;
            for(i=0;i<ne;i++)
                if(g_pbp_neq.epoch_has_neq && g_pbp_neq.epoch_has_neq[i])
                    g_pbp_neq.day1_fixed_clock[i]-=mean_off;
            fprintf(stderr,"[PBP-NEQ] datum(fixed): %.3f ns from %d epochs\n",
                    mean_off*1e9/CLIGHT,cnt_off_fix);
        }
        if(cnt_off_flt>0){
            double mean_off=sum_off_flt/cnt_off_flt;
            for(i=0;i<ne;i++)
                if(g_pbp_neq.epoch_has_neq && g_pbp_neq.epoch_has_neq[i])
                    g_pbp_neq.day1_neqfloat_clock[i]-=mean_off;
            fprintf(stderr,"[PBP-NEQ] datum(neqfloat): %.3f ns from %d epochs\n",
                    mean_off*1e9/CLIGHT,cnt_off_flt);
        }
        fprintf(stderr,"[PBP-NEQ] clk: fixed=%d float_fb=%d none=%d / %d\n",
                n_fixed,n_float_fb,n_none,ne);
    }
    free(xhat_float); free(xhat_fixed);

    g_pbp_neq.ready=1;

    /* ── Diagnostic summary ────────────────────────────────────────────── */
    {
        int total_bs=g_pbp_diag.backsub_ok_count+g_pbp_diag.backsub_fallback_count;
        double fb_ratio = total_bs>0
            ? (double)g_pbp_diag.backsub_fallback_count/total_bs*100.0 : 0.0;
        fprintf(stderr,
            "\n[PBP-DIAG] ══════════ NEQ SUMMARY ══════════\n"
            "[PBP-DIAG] add_epoch: ok=%d / total=%d (%.1f%%)\n"
            "[PBP-DIAG]   fail_bad_input      = %d\n"
            "[PBP-DIAG]   fail_epoch_id        = %d\n"
            "[PBP-DIAG]   fail_unconverged     = %d (vrms>=50m, output float)\n"
            "[PBP-DIAG]   fail_R_invert        = %d\n"
            "[PBP-DIAG]   fail_no_active_param = %d\n"
            "[PBP-DIAG]   fail_alloc           = %d\n"
            "[PBP-DIAG] epoch_collision_total  = %d\n"
            "[PBP-DIAG] output: fixed=%d  float_fallback=%d  ratio=%.1f%%\n"
            "[PBP-DIAG] ══════════════════════════════════\n",
            g_pbp_diag.add_ok_total, g_pbp_diag.add_call_total,
            g_pbp_diag.add_call_total>0
                ? 100.0*g_pbp_diag.add_ok_total/g_pbp_diag.add_call_total : 0.0,
            g_pbp_diag.fail_bad_input,
            g_pbp_diag.fail_epoch_id,
            g_pbp_diag.fail_vrms_too_large,
            g_pbp_diag.fail_R_invert,
            g_pbp_diag.fail_no_active_param,
            g_pbp_diag.fail_alloc,
            g_pbp_diag.epoch_collision_total,
            g_pbp_diag.backsub_ok_count,g_pbp_diag.backsub_fallback_count,fb_ratio);
    }
    return 1;
}

extern int pbp_write_day1_fixed_clock_file(const char *path)
{
    const double M2NS=1e9/CLIGHT;
    int e,n_w=0;
    if(!path||!*path||!g_pbp_neq.ready||g_pbp_neq.day1_epoch_count<=0)return 0;
    createdir(path);
    FILE *fp=fopen(path,"w");if(!fp)return 0;
    fprintf(fp,"# PBP fixed GPS clock (day1, all-in-one Cholesky)  unit: ns\n");
    for(e=0;e<g_pbp_neq.day1_epoch_count;e++){
        if(!g_pbp_neq.epoch_valid||!g_pbp_neq.epoch_valid[e])continue;
        gtime_t t=g_pbp_neq.epoch_time[e];double ep[6];time2epoch(t,ep);
        fprintf(fp,"%04d/%02d/%02d %02d:%02d:%02d %.3f\n",
                (int)ep[0],(int)ep[1],(int)ep[2],
                (int)ep[3],(int)ep[4],(int)floor(ep[5]+0.5),
                g_pbp_neq.day1_fixed_clock[e]*M2NS);
        n_w++;
    }
    fclose(fp);
    fprintf(stderr,"[PBP-NEQ] wrote: %s (%d ep)\n",path,n_w);
    return 1;
}

extern int pbp_write_day1_neqfloat_clock_file(const char *path)
{
    const double M2NS=1e9/CLIGHT;
    int e,n_w=0;
    if(!path||!*path||!g_pbp_neq.ready||g_pbp_neq.day1_epoch_count<=0)return 0;
    if(!g_pbp_neq.day1_neqfloat_clock)return 0;
    createdir(path);
    FILE *fp=fopen(path,"w");if(!fp)return 0;
    fprintf(fp,"# PBP NEQ float GPS clock (day1, no DD constraints)  unit: ns\n");
    for(e=0;e<g_pbp_neq.day1_epoch_count;e++){
        if(!g_pbp_neq.epoch_valid||!g_pbp_neq.epoch_valid[e])continue;
        gtime_t t=g_pbp_neq.epoch_time[e];double ep[6];time2epoch(t,ep);
        fprintf(fp,"%04d/%02d/%02d %02d:%02d:%02d %.3f\n",
                (int)ep[0],(int)ep[1],(int)ep[2],
                (int)ep[3],(int)ep[4],(int)floor(ep[5]+0.5),
                g_pbp_neq.day1_neqfloat_clock[e]*M2NS);
        n_w++;
    }
    fclose(fp);
    fprintf(stderr,"[PBP-NEQ] wrote neqfloat: %s (%d ep)\n",path,n_w);
    return 1;
}

extern int pbp_get_fixed_clock(gtime_t t,int sys_idx,double *clk)
{
    if(!g_pbp_neq.ready||!clk||sys_idx!=0)return 0;
    int e=pbp_epoch_id_day1(t);
    if(e<0||e>=g_pbp_neq.n_epoch||!g_pbp_neq.epoch_valid[e])return 0;
    *clk=g_pbp_neq.day1_fixed_clock[e];return 1;
}

extern int pbp_get_fixed_arc_bias(gtime_t t,int sat,double *bias,double *var)
{(void)t;(void)sat;(void)bias;(void)var;return 0;}

/* ── Diagnostic CSV output (requirement 6) ─────────────────────────────── */
extern int pbp_write_diag_csv(const char *path)
{
    int e,ne;
    if(!path||!*path)return 0;
    ne=g_pbp_neq.n_epoch;
    if(ne<=0)return 0;
    createdir(path);
    FILE *fp=fopen(path,"w");
    if(!fp){fprintf(stderr,"[PBP-DIAG] cannot open %s\n",path);return 0;}
    fprintf(fp,"epoch_id,time,hit_count,nv_last,nact_last,backsub_ok\n");
    for(e=0;e<ne&&e<PBP_DIAG_MAX_EPOCH;e++){
        gtime_t t=g_pbp_neq.epoch_time[e];
        double ep[6]; time2epoch(t,ep);
        fprintf(fp,"%d,%04d/%02d/%02d %02d:%02d:%02d,%d,%d,%d,%d\n",
                e,(int)ep[0],(int)ep[1],(int)ep[2],
                (int)ep[3],(int)ep[4],(int)floor(ep[5]+0.5),
                g_pbp_diag.epoch_hit_count[e],
                g_pbp_diag.epoch_nv_last[e],
                g_pbp_diag.epoch_nact_last[e],
                g_pbp_diag.epoch_backsub_ok[e]);
    }
    fclose(fp);
    fprintf(stderr,"[PBP-DIAG] wrote: %s (%d rows)\n",path,ne<PBP_DIAG_MAX_EPOCH?ne:PBP_DIAG_MAX_EPOCH);
    return 1;
}

/* ── Legacy stubs ─────────────────────────────────────────────────────── */
extern int pbp_apply_session_pseudoobs(rtk_t *rtk){(void)rtk;return 0;}
extern int apply_ar_fixed(rtk_t *rtk,const ddamb_t *ddamb,int n_dd)
{(void)rtk;(void)ddamb;(void)n_dd;return 0;}
extern int pbp_bds_is_sidereal(int sat)
{
    int prn=0;
    if(satsys(sat,&prn)!=SYS_CMP)return 0;
    if(prn>=1&&prn<=10)return 1;
    if(prn==13||prn==16)return 1;
    if(prn>=38&&prn<=40)return 1;
    if(prn>=59&&prn<=62)return 1;
    return 0;
}
