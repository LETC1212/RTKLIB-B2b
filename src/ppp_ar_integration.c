/*------------------------------------------------------------------------------
* ppp_ar_integration.c : Integration wrapper for Pass-by-Pass AR
*
* v3.0 changes (pbp_write_dd_wlnl):
*   [ADDED] 5-bin fraction distribution for DD WL and DD NL
*   [ADDED] GPS / BDS per-system split statistics
*   [ADDED] Both "all independent DD" and "fixed only" sets
*   [KEPT]  Original |frac|<0.2 counting intact
*
*   Bin definitions (pbp_frac output = signed distance to nearest integer):
*     bin 0: [-0.2, 0.2]    near integer
*     bin 1: [-0.4,-0.2)    left mid-range
*     bin 2: ( 0.2, 0.4]    right mid-range
*     bin 3: (-inf,-0.4)    far left
*     bin 4: ( 0.4,+inf)    far right
*
*   BDS count = entries already in ddamb[] (MEO sats failing arc-pairing excluded)
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include "ppp_ar_passbypass.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>

/* signed fractional part to nearest integer  (-0.5, 0.5] */
static double pbp_frac(double x)
{
    double r = x - floor(x + 0.5);
    if (r <= -0.5) r += 1.0;
    if (r >   0.5) r -= 1.0;
    return r;
}

/* classify frac into bin 0-4 */
static int frac_bin(double frac)
{
    if (frac >= -0.2 && frac <=  0.2) return 0;
    if (frac >= -0.4 && frac <  -0.2) return 1;
    if (frac >   0.2 && frac <=  0.4) return 2;
    if (frac <  -0.4)                  return 3;
    return 4;
}

/* write one 5-bin percentage line */
static void print_frac_bins(FILE *fp, const char *prefix,
                             const int cnt[5], int total)
{
    if (!fp || total <= 0) return;
    double t = (double)total;
    fprintf(fp,
        "%s: [-0.2,0.2]=%5.1f%%  [-0.4,-0.2)=%5.1f%%  (0.2,0.4]=%5.1f%%"
        "  <-0.4=%5.1f%%  >0.4=%5.1f%%  (n=%d)\n",
        prefix,
        100.0*cnt[0]/t, 100.0*cnt[1]/t, 100.0*cnt[2]/t,
        100.0*cnt[3]/t, 100.0*cnt[4]/t, total);
}

/*------------------------------------------------------------------------------
* pbp_write_dd_wlnl : write DD WL/NL table + extended statistics
*
*   [NEW §2] 5-bin fraction distribution for all/fixed, total + GPS + BDS
*   [KEPT §1] Original |frac|<0.2 summary block unchanged
*-----------------------------------------------------------------------------*/
static int pbp_write_dd_wlnl(FILE *fp, const ddamb_t *dd, int ndd, const char *tag)
{
    if (!fp || !dd || ndd <= 0) return 0;

    /* legacy counters */
    int n_wl_all=0, n_wl_02=0, n_wl_fix=0, n_wl_fix_02=0;
    int n_nl_all=0, n_nl_02=0, n_nl_fix=0, n_nl_fix_02=0;

    /* 5-bin counters: [sys_idx][subset][bin]
     *   sys_idx : 0=all  1=GPS  2=BDS
     *   subset  : 0=all independent DD  1=fixed only             */
    int wl_bins[3][2][5], nl_bins[3][2][5];
    int wl_cnt[3][2],     nl_cnt[3][2];
    memset(wl_bins, 0, sizeof(wl_bins));
    memset(nl_bins, 0, sizeof(nl_bins));
    memset(wl_cnt,  0, sizeof(wl_cnt));
    memset(nl_cnt,  0, sizeof(nl_cnt));

    fprintf(fp, "\n# ===== %s =====\n", tag ? tag : "DD_WL/DD_NL");
    fprintf(fp, "# sat_ref sat_sat arc_ref arc_sat "
                "DD_WL frac_WL fixed_WL DD_NL frac_NL fixed_NL DD_IF_fix(m)\n");

    for (int i = 0; i < ndd; i++) {
        if (dd[i].sat1 <= 0 || dd[i].sat2 <= 0) continue;

        char sref[16]="", ssat[16]="";
        satno2id(dd[i].sat1, sref);
        satno2id(dd[i].sat2, ssat);

        const double wl  = dd[i].DD_WL;
        const double nl  = dd[i].DD_NL;
        const double fwl = pbp_frac(wl);
        const double fnl = pbp_frac(nl);

        fprintf(fp, "%s %s %d %d %.6f %.6f %d %.6f %.6f %d %.6f\n",
                sref, ssat, dd[i].arc1, dd[i].arc2,
                wl, fwl, dd[i].fixed_WL ? 1 : 0,
                nl, fnl, dd[i].fixed_NL ? 1 : 0,
                dd[i].DD_IF_fix);

        /* legacy |frac|<0.2 */
        const int wl_ok = (fabs(fwl) < 0.2);
        const int nl_ok = (fabs(fnl) < 0.2);
        n_wl_all++; if (wl_ok) n_wl_02++;
        n_nl_all++; if (nl_ok) n_nl_02++;
        if (dd[i].fixed_WL) { n_wl_fix++; if (wl_ok) n_wl_fix_02++; }
        if (dd[i].fixed_NL) { n_nl_fix++; if (nl_ok) n_nl_fix_02++; }

        /* 5-bin accounting */
        int sys_sat = satsys(dd[i].sat1, NULL);
        int si = (sys_sat==SYS_GPS) ? 1 : (sys_sat==SYS_CMP) ? 2 : 0;
        int bwl = frac_bin(fwl);
        int bnl = frac_bin(fnl);

        /* WL: all independent DD */
        wl_bins[0][0][bwl]++; wl_bins[si][0][bwl]++;
        wl_cnt[0][0]++;       wl_cnt[si][0]++;
        /* WL: fixed */
        if (dd[i].fixed_WL) {
            wl_bins[0][1][bwl]++; wl_bins[si][1][bwl]++;
            wl_cnt[0][1]++;       wl_cnt[si][1]++;
        }
        /* NL: count when WL was fixed (DD_NL is then meaningful) */
        if (dd[i].fixed_WL) {
            nl_bins[0][0][bnl]++; nl_bins[si][0][bnl]++;
            nl_cnt[0][0]++;       nl_cnt[si][0]++;
        }
        /* NL: fully fixed */
        if (dd[i].fixed_NL) {
            nl_bins[0][1][bnl]++; nl_bins[si][1][bnl]++;
            nl_cnt[0][1]++;       nl_cnt[si][1]++;
        }
    }

    /* ── §1  Legacy |frac|<0.2 summary (unchanged) ─────────────────────── */
    fprintf(fp, "# --- Summary (|frac|<0.2) ---\n");
    if (n_wl_all > 0)
        fprintf(fp, "# WL: all=%d, in0.2=%d (%.2f%%)\n",
                n_wl_all, n_wl_02, 100.0*(double)n_wl_02/(double)n_wl_all);
    if (n_wl_fix > 0)
        fprintf(fp, "# WL: fixed=%d, in0.2=%d (%.2f%%)\n",
                n_wl_fix, n_wl_fix_02, 100.0*(double)n_wl_fix_02/(double)n_wl_fix);
    if (n_nl_all > 0)
        fprintf(fp, "# NL: all=%d, in0.2=%d (%.2f%%)\n",
                n_nl_all, n_nl_02, 100.0*(double)n_nl_02/(double)n_nl_all);
    if (n_nl_fix > 0)
        fprintf(fp, "# NL: fixed=%d, in0.2=%d (%.2f%%)\n",
                n_nl_fix, n_nl_fix_02, 100.0*(double)n_nl_fix_02/(double)n_nl_fix);

    /* ── §2  5-bin: all systems ─────────────────────────────────────────── */
    fprintf(fp, "# --- Fraction bins (all independent DD) ---\n");
    print_frac_bins(fp, "# WL bins", wl_bins[0][0], wl_cnt[0][0]);
    print_frac_bins(fp, "# NL bins", nl_bins[0][0], nl_cnt[0][0]);
    fprintf(fp, "# --- Fraction bins (fixed only) ---\n");
    print_frac_bins(fp, "# WL bins", wl_bins[0][1], wl_cnt[0][1]);
    print_frac_bins(fp, "# NL bins", nl_bins[0][1], nl_cnt[0][1]);

    /* ── §3  GPS ─────────────────────────────────────────────────────────── */
    if (wl_cnt[1][0] > 0 || wl_cnt[1][1] > 0) {
        fprintf(fp, "# --- GPS Fraction bins (all independent DD) ---\n");
        print_frac_bins(fp, "# GPS WL bins", wl_bins[1][0], wl_cnt[1][0]);
        print_frac_bins(fp, "# GPS NL bins", nl_bins[1][0], nl_cnt[1][0]);
        fprintf(fp, "# --- GPS Fraction bins (fixed only) ---\n");
        print_frac_bins(fp, "# GPS WL bins", wl_bins[1][1], wl_cnt[1][1]);
        print_frac_bins(fp, "# GPS NL bins", nl_bins[1][1], nl_cnt[1][1]);
    }

    /* ── §4  BDS (~1-day repeat sats only, others excluded by arc-pairing) ── */
    if (wl_cnt[2][0] > 0 || wl_cnt[2][1] > 0) {
        fprintf(fp,
            "# --- BDS Fraction bins (all, ~1-day sats after arc-pairing filter) ---\n");
        print_frac_bins(fp, "# BDS WL bins", wl_bins[2][0], wl_cnt[2][0]);
        print_frac_bins(fp, "# BDS NL bins", nl_bins[2][0], nl_cnt[2][0]);
        fprintf(fp, "# --- BDS Fraction bins (fixed only) ---\n");
        print_frac_bins(fp, "# BDS WL bins", wl_bins[2][1], wl_cnt[2][1]);
        print_frac_bins(fp, "# BDS NL bins", nl_bins[2][1], nl_cnt[2][1]);
    }

    return 1;
}

/* ── Driver flags ─────────────────────────────────────────────────────────── */
int pbp_day_tag      = -1;
int pbp_collect_flag =  0;
int pbp_apply_flag   =  0;
extern int pbp_resolve_flag;

extern satamb_t satamb[MAXSAT];
extern int      n_ddamb;
extern ddamb_t  ddamb[MAXSAT*MAXSAT];
extern int      refsat;

extern void print_arc_summary(void);
extern int  compute_dd_ambiguities(const satamb_t *satamb, int refsat, ddamb_t *ddamb, int *n_dd);
extern int  fix_wl_nl_ambiguities(ddamb_t *ddamb, int n_dd);
extern void pbp_clear_fixed_constraints(void);
extern int  pbp_store_fixed_constraints(const ddamb_t *ddamb, int n_dd, double Pb);
extern int  collect_ambiguities(const rtk_t *rtk, const obsd_t *obs, int n, int day, satamb_t *satamb);
extern int pbp_bds_is_sidereal(int sat);

/* ── Fixed DD IF day-2 sorted output ─────────────────────────────────────────
 * This diagnostic writes ALL successfully fixed DD WL/NL constraints.  The
 * previous epoch-expanded version required the day-1 reference arc and the
 * day-1 target arc to overlap at the same epoch; that condition is too strict
 * for pass-by-pass arc-level DD constraints and can make the output contain
 * only one DD per system although many DDs were actually fixed.
 *
 * Here each fixed DD is written once and sorted by its second-day target arc
 * start time.  Both the reference satellite and target satellite day-2 arc
 * windows are printed, so the user can see where in day 2 this fixed DD comes
 * from without discarding non-overlapping arc pairs.
 *-----------------------------------------------------------------------------*/
#ifndef PBP_SIDEREAL_SHIFT_SEC
#define PBP_SIDEREAL_SHIFT_SEC 85920.0
#endif

typedef struct {
    gtime_t key_time;
    int sat1, sat2;
    int arc1_d0, arc1_d1;
    int arc2_d0, arc2_d1;
    int fixed_wl, fixed_nl;
    double dd_if, dd_if_fix;
    double dd_wl, dd_wl_fix;
    double dd_nl, dd_nl_fix;
    double var_dd_if, var_dd_wl, var_dd_nl;
} pbp_fixed_ddif_rec_t;

static int pbp_time_is_valid(gtime_t t)
{
    return t.time != 0 || t.sec != 0.0;
}

static double pbp_overlap_sec(gtime_t a0, gtime_t a1, gtime_t b0, gtime_t b1)
{
    gtime_t s = timediff(a0,b0) >= 0.0 ? a0 : b0;
    gtime_t e = timediff(a1,b1) <= 0.0 ? a1 : b1;
    double dt = timediff(e,s);
    return dt > 0.0 ? dt : 0.0;
}

static int pbp_find_day1_arc_from_day0(int sat, int arc0)
{
    if (sat <= 0 || sat > MAXSAT) return -1;
    if (arc0 < 0 || arc0 >= satamb[sat-1].n) return -1;

    const ambarc_t *a0 = &satamb[sat-1].arc[arc0];
    if (a0->day != 0) return -1;

    gtime_t sh_s = timeadd(a0->ts, PBP_SIDEREAL_SHIFT_SEC);
    gtime_t sh_e = timeadd(a0->te, PBP_SIDEREAL_SHIFT_SEC);

    int best = -1;
    double best_ov = -1.0;
    for (int i = 0; i < satamb[sat-1].n; i++) {
        const ambarc_t *a1 = &satamb[sat-1].arc[i];
        if (a1->day != 1 || a1->nobs < 10) continue;
        double ov = pbp_overlap_sec(sh_s, sh_e, a1->ts, a1->te);
        if (ov > best_ov) {
            best_ov = ov;
            best = i;
        }
    }
    if (best >= 0 && best_ov > 0.0) return best;

    /* Fallback: choose the closest day-1 arc by midpoint.  This is only for
     * diagnostic output; the fixed DD itself has already been accepted by
     * compute_dd_ambiguities() and fix_wl_nl_ambiguities(). */
    double mid0 = time2gpst(sh_s, NULL) + 0.5 * timediff(sh_e, sh_s);
    double best_abs = 1.0e100;
    best = -1;
    for (int i = 0; i < satamb[sat-1].n; i++) {
        const ambarc_t *a1 = &satamb[sat-1].arc[i];
        if (a1->day != 1 || a1->nobs < 10) continue;
        double mid1 = time2gpst(a1->ts, NULL) + 0.5 * timediff(a1->te, a1->ts);
        double d = fabs(mid1 - mid0);
        if (d < best_abs) {
            best_abs = d;
            best = i;
        }
    }
    return best;
}

static int pbp_cmp_fixed_ddif_rec(const void *pa, const void *pb)
{
    const pbp_fixed_ddif_rec_t *a = (const pbp_fixed_ddif_rec_t *)pa;
    const pbp_fixed_ddif_rec_t *b = (const pbp_fixed_ddif_rec_t *)pb;
    double dt = timediff(a->key_time, b->key_time);
    if (dt < -1E-9) return -1;
    if (dt >  1E-9) return  1;
    if (a->sat1 != b->sat1) return a->sat1 - b->sat1;
    if (a->sat2 != b->sat2) return a->sat2 - b->sat2;
    return a->arc2_d0 - b->arc2_d0;
}

static const char *pbp_sys_name_for_sat(int sat)
{
    int sys = satsys(sat, NULL);
    if (sys == SYS_GPS) return "GPS";
    if (sys == SYS_CMP) return "BDS";
    if (sys == SYS_GLO) return "GLO";
    if (sys == SYS_GAL) return "GAL";
    if (sys == SYS_QZS) return "QZS";
    return "UNK";
}

static void pbp_time_to_str_or_dash(gtime_t t, char *s, int n)
{
    if (!s || n <= 0) return;
    if (!pbp_time_is_valid(t)) {
        snprintf(s, n, "-");
        return;
    }
    time2str(t, s, 1);
}

static int pbp_write_fixed_ddif_day2(const char *path, const ddamb_t *dd, int ndd)
{
    FILE *fp = NULL;
    pbp_fixed_ddif_rec_t *rec = NULL;
    int nrec = 0, nfixed = 0;
    char cwd[PATH_MAX] = "";

    if (!path || !*path || !dd || ndd <= 0) return 0;

    rec = (pbp_fixed_ddif_rec_t *)calloc((size_t)ndd, sizeof(pbp_fixed_ddif_rec_t));
    if (!rec) return 0;

    for (int i = 0; i < ndd; i++) {
        if (!dd[i].fixed_WL || !dd[i].fixed_NL) continue;
        if (dd[i].sat1 <= 0 || dd[i].sat2 <= 0) continue;
        nfixed++;

        memset(&rec[nrec], 0, sizeof(rec[nrec]));
        rec[nrec].sat1 = dd[i].sat1;
        rec[nrec].sat2 = dd[i].sat2;
        rec[nrec].arc1_d0 = dd[i].arc1;
        rec[nrec].arc2_d0 = dd[i].arc2;
        rec[nrec].arc1_d1 = pbp_find_day1_arc_from_day0(dd[i].sat1, dd[i].arc1);
        rec[nrec].arc2_d1 = pbp_find_day1_arc_from_day0(dd[i].sat2, dd[i].arc2);
        rec[nrec].fixed_wl = dd[i].fixed_WL ? 1 : 0;
        rec[nrec].fixed_nl = dd[i].fixed_NL ? 1 : 0;
        rec[nrec].dd_if = dd[i].DD_IF;
        rec[nrec].dd_if_fix = dd[i].DD_IF_fix;
        rec[nrec].dd_wl = dd[i].DD_WL;
        rec[nrec].dd_wl_fix = dd[i].DD_WL_fix;
        rec[nrec].dd_nl = dd[i].DD_NL;
        rec[nrec].dd_nl_fix = dd[i].DD_NL_fix;
        rec[nrec].var_dd_if = dd[i].var_DD_IF;
        rec[nrec].var_dd_wl = dd[i].var_DD_WL;
        rec[nrec].var_dd_nl = dd[i].var_DD_NL;

        if (rec[nrec].arc2_d1 >= 0) {
            rec[nrec].key_time = satamb[dd[i].sat2-1].arc[rec[nrec].arc2_d1].ts;
        }
        else if (rec[nrec].arc1_d1 >= 0) {
            rec[nrec].key_time = satamb[dd[i].sat1-1].arc[rec[nrec].arc1_d1].ts;
        }
        else {
            rec[nrec].key_time = (gtime_t){0};
        }
        nrec++;
    }

    qsort(rec, nrec, sizeof(rec[0]), pbp_cmp_fixed_ddif_rec);

    fp = fopen(path, "w");
    if (!fp) { free(rec); return 0; }
    if (!getcwd(cwd, sizeof(cwd))) strncpy(cwd, ".", sizeof(cwd)-1);

    fprintf(fp, "%% PBP fixed DD IF table sorted by day-2 arc time (GPST)\n");
    fprintf(fp, "%% CWD: %s\n", cwd);
    fprintf(fp, "%% IMPORTANT: all fixed DD constraints are written once; no epoch-overlap filter is applied.\n");
    fprintf(fp, "%% DD_IF_FIX_M is the fixed DD IF pseudo-observation in metres; DD_WL_INT and DD_NL_INT are the fixed integers.\n");
    fprintf(fp, "%% FIX_TIME_KEY REF_D2_START REF_D2_END SAT_D2_START SAT_D2_END SYS REF SAT ");
    fprintf(fp, "ARC_REF_D0 ARC_REF_D1 ARC_SAT_D0 ARC_SAT_D1 DD_IF_FLOAT_M DD_IF_FIX_M DD_IF_RES_M DD_WL_FLOAT DD_WL_INT DD_NL_FLOAT DD_NL_INT SIG_DD_IF_M SIG_DD_WL_CYC SIG_DD_NL_CYC\n");

    for (int i = 0; i < nrec; i++) {
        char key[64], rts[64], rte[64], sts[64], ste[64];
        char sref[16] = "", ssat[16] = "";
        gtime_t rt0 = {0}, rt1 = {0}, st0 = {0}, st1 = {0};
        if (rec[i].arc1_d1 >= 0) {
            rt0 = satamb[rec[i].sat1-1].arc[rec[i].arc1_d1].ts;
            rt1 = satamb[rec[i].sat1-1].arc[rec[i].arc1_d1].te;
        }
        if (rec[i].arc2_d1 >= 0) {
            st0 = satamb[rec[i].sat2-1].arc[rec[i].arc2_d1].ts;
            st1 = satamb[rec[i].sat2-1].arc[rec[i].arc2_d1].te;
        }
        pbp_time_to_str_or_dash(rec[i].key_time, key, sizeof(key));
        pbp_time_to_str_or_dash(rt0, rts, sizeof(rts));
        pbp_time_to_str_or_dash(rt1, rte, sizeof(rte));
        pbp_time_to_str_or_dash(st0, sts, sizeof(sts));
        pbp_time_to_str_or_dash(st1, ste, sizeof(ste));
        satno2id(rec[i].sat1, sref);
        satno2id(rec[i].sat2, ssat);

        fprintf(fp,
            "%s %s %s %s %s %3s %3s %3s %10d %10d %10d %10d "
            "%14.6f %14.6f %12.6f %12.6f %10.0f %12.6f %10.0f %12.6f %12.6f %12.6f\n",
            key, rts, rte, sts, ste,
            pbp_sys_name_for_sat(rec[i].sat1), sref, ssat,
            rec[i].arc1_d0, rec[i].arc1_d1, rec[i].arc2_d0, rec[i].arc2_d1,
            rec[i].dd_if, rec[i].dd_if_fix, rec[i].dd_if_fix - rec[i].dd_if,
            rec[i].dd_wl, rec[i].dd_wl_fix,
            rec[i].dd_nl, rec[i].dd_nl_fix,
            rec[i].var_dd_if > 0.0 ? sqrt(rec[i].var_dd_if) : 0.0,
            rec[i].var_dd_wl > 0.0 ? sqrt(rec[i].var_dd_wl) : 0.0,
            rec[i].var_dd_nl > 0.0 ? sqrt(rec[i].var_dd_nl) : 0.0);
    }

    fprintf(fp, "%% summary: input_dd=%d fixed_dd_seen=%d output_fixed_dd=%d\n", ndd, nfixed, nrec);
    fclose(fp);
    free(rec);
    return nrec;
}

static int select_refsat_auto_sys(const satamb_t *sa, int sys)
{
    int ref = 0, max_obs = 0;
    for (int i = 0; i < MAXSAT; i++) {
        if (sa[i].n <= 0) continue;
        int sat = i + 1;
        if (satsys(sat, NULL) != sys) continue;
        int has0=0, has1=0, tot=0;
        for (int j = 0; j < sa[i].n; j++) {
            if(!pbp_bds_is_sidereal(sat)&&satsys(sat,NULL)==SYS_CMP)  continue;
            if (sa[i].arc[j].day==0 && sa[i].arc[j].nobs>=10) { has0=1; tot+=sa[i].arc[j].nobs; }
            if (sa[i].arc[j].day==1 && sa[i].arc[j].nobs>=10) { has1=1; tot+=sa[i].arc[j].nobs; }
        }
        if (has0 && has1 && tot > max_obs) { max_obs=tot; ref=sat; }
    }
    return ref;
}

/*------------------------------------------------------------------------------
* ppp_ar_48h
*-----------------------------------------------------------------------------*/
extern int ppp_ar_48h(const prcopt_t *popt, rtk_t *rtk, const obs_t *obs)
{
    int n_fixed = 0;
    (void)obs; (void)rtk;
    if (!popt) return 0;
    if (popt->armode_pbp == 0) {
        trace(2, "ppp_ar_48h: disabled\n"); return 0;
    }

    trace(1, "\n========== Pass-by-Pass Ambiguity Resolution ==========\n");
    printf("\n========== Pass-by-Pass Ambiguity Resolution ==========\n");
    print_arc_summary();

    int ref_gps=0, ref_cmp=0;
    if (popt->pbp_refsat != 0) {
        refsat = popt->pbp_refsat;
        if (refsat<=0||refsat>MAXSAT) {
            printf("Error: invalid pbp_refsat=%d\n", popt->pbp_refsat); return 0;
        }
        ref_gps = (satsys(refsat,NULL)==SYS_GPS) ? refsat : 0;
        ref_cmp = (satsys(refsat,NULL)==SYS_CMP) ? refsat : 0;
    }
    else {
        ref_gps = select_refsat_auto_sys(satamb, SYS_GPS);
        ref_cmp = select_refsat_auto_sys(satamb, SYS_CMP);
        refsat  = ref_gps ? ref_gps : ref_cmp;
    }

    if (refsat == 0) {
        printf("Error: no satellite has both day0 and day1 data (>=10 obs)\n");
        return 0;
    }
    if (ref_gps) { char s[16]=""; satno2id(ref_gps,s); printf("Reference (GPS): %s\n",s); }
    if (ref_cmp) { char s[16]=""; satno2id(ref_cmp,s); printf("Reference (BDS): %s\n",s); }

    /* BDS arc diagnostic: print every BDS satellite's arc info to stderr */
    if (ref_cmp == 0) {
        fprintf(stderr, "[PBP] WARNING: No BDS reference satellite found!\n");
        fprintf(stderr, "[PBP] BDS arc summary:\n");
        for (int i = 0; i < MAXSAT; i++) {
            if (satamb[i].n == 0) continue;
            int sat = i+1;
            if (satsys(sat, NULL) != SYS_CMP) continue;
            char sid[8]=""; satno2id(sat, sid);
            int has0=0, has1=0;
            for (int j = 0; j < satamb[i].n; j++) {
                if (satamb[i].arc[j].day==0 && satamb[i].arc[j].nobs>=10) has0=1;
                if (satamb[i].arc[j].day==1 && satamb[i].arc[j].nobs>=10) has1=1;
            }
            fprintf(stderr, "[PBP]   %s: %d arcs  has_day0=%d has_day1=%d\n",
                    sid, satamb[i].n, has0, has1);
            for (int j = 0; j < satamb[i].n; j++) {
                fprintf(stderr, "[PBP]     arc%d: day=%d nobs=%d "
                        "N_WL=%.3f var_WL=%.6f N_IF=%.3f var_IF=%.2e\n",
                        j, satamb[i].arc[j].day, satamb[i].arc[j].nobs,
                        satamb[i].arc[j].N_WL, satamb[i].arc[j].var_WL,
                        satamb[i].arc[j].N_IF, satamb[i].arc[j].var_IF);
            }
        }
    }

    /* diagnostic file */
    char cwd[PATH_MAX]="", ddpath[PATH_MAX]="pbp_dd_wlnl.txt";
    if (!getcwd(cwd,sizeof(cwd))) strncpy(cwd,".",sizeof(cwd)-1);
    FILE *fpdd = fopen(ddpath,"w");
    if (!fpdd) {
        snprintf(ddpath,sizeof(ddpath),"/tmp/pbp_dd_wlnl_%d.txt",(int)getpid());
        fpdd = fopen(ddpath,"w");
    }
    if (fpdd) {
        fprintf(fpdd,"# PBP arc-level DD WL/NL (day0-day1)\n");
        fprintf(fpdd,"# CWD: %s\n",cwd);
        fprintf(fpdd,"# BDS = sats passing arc-pairing (non-1-day repeat excluded)\n");
        fflush(fpdd);   /* ensure header is written even if program exits early */
    }

    n_ddamb = 0;
    int n_fixed_total=0, n_dd=0;

    if (ref_gps) {
        n_dd=0;
        if (compute_dd_ambiguities(satamb,ref_gps,ddamb+n_ddamb,&n_dd)&&n_dd>0) {
            int nfx=fix_wl_nl_ambiguities(ddamb+n_ddamb,n_dd);
            printf("DD/fixed GPS: %d/%d\n",nfx,n_dd);
            if (fpdd) pbp_write_dd_wlnl(fpdd,ddamb+n_ddamb,n_dd,"GPS");
            n_fixed_total+=nfx; n_ddamb+=n_dd;
        }
    }
    if (ref_cmp) {
        n_dd=0;
        if (compute_dd_ambiguities(satamb,ref_cmp,ddamb+n_ddamb,&n_dd)&&n_dd>0) {
            int nfx=fix_wl_nl_ambiguities(ddamb+n_ddamb,n_dd);
            printf("DD/fixed BDS: %d/%d\n",nfx,n_dd);
            if (fpdd) pbp_write_dd_wlnl(fpdd,ddamb+n_ddamb,n_dd,"BDS");
            n_fixed_total+=nfx; n_ddamb+=n_dd;
        }
    }

    printf("DD/fixed TOTAL: %d/%d\n",n_fixed_total,n_ddamb);
    n_fixed=n_fixed_total;

    if (n_fixed<=0||n_ddamb<=0) {
        printf("Warning: no ambiguities fixed\n");
        if (fpdd) { fprintf(fpdd,"\n# no ambiguities fixed\n"); fclose(fpdd); }
        return 0;
    }

    /* store for NEQ re-solve */
    {
        int n_ind=pbp_store_fixed_constraints(ddamb,n_ddamb,1.0e10);
        printf("Stored independent fixed DD constraints: %d\n",n_ind);
    }

    {
        char ddifpath[PATH_MAX] = "pbp_fixed_ddif_day2.txt";
        int nw = pbp_write_fixed_ddif_day2(ddifpath, ddamb, n_ddamb);
        if (nw > 0) {
            printf("PBP fixed DD IF day2 sorted table written: %s  (rows=%d)\n", ddifpath, nw);
        }
        else {
            printf("[PBP] Warning: no fixed DD IF row written to %s\n", ddifpath);
        }
    }

    if (fpdd) { fclose(fpdd); printf("PBP DD WL/NL written: %s\n",ddpath); }
    printf("=======================================================\n\n");
    return n_fixed;
}

/*------------------------------------------------------------------------------
* collect_ambiguities_epoch
*-----------------------------------------------------------------------------*/
extern int collect_ambiguities_epoch(const rtk_t *rtk, const obsd_t *obs,
                                      int n, int day)
{
    return collect_ambiguities(rtk, obs, n, day, satamb);
}
