/*------------------------------------------------------------------------------
* postpos.c : post-processing positioning
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/05/08  1.0  new
*           2008/06/16  1.1  support binary inputs
*           2009/01/02  1.2  support new rtk positioning api
*           2009/09/03  1.3  fix bug on combined mode of moving-baseline
*           2009/12/04  1.4  fix bug on obs data buffer overflow
*           2010/07/26  1.5  support ppp-kinematic and ppp-static
*                            support multiple sessions
*                            support sbas positioning
*                            changed api:
*                                postpos()
*                            deleted api:
*                                postposopt()
*           2010/08/16  1.6  fix bug sbas message synchronization (2.4.0_p4)
*           2010/12/09  1.7  support qzss lex and ssr corrections
*           2011/02/07  1.8  fix bug on sbas navigation data conflict
*           2011/03/22  1.9  add function reading g_tec file
*           2011/08/20  1.10 fix bug on freez if solstatic=single and combined
*           2011/09/15  1.11 add function reading stec file
*           2012/02/01  1.12 support keyword expansion of rtcm ssr corrections
*           2013/03/11  1.13 add function reading otl and erp data
*           2014/06/29  1.14 fix problem on overflow of # of satellites
*           2015/03/23  1.15 fix bug on ant type replacement by rinex header
*                            fix bug on combined filter for moving-base mode
*           2015/04/29  1.16 fix bug on reading rtcm ssr corrections
*                            add function to read satellite fcb
*                            add function to read stec and troposphere file
*                            add keyword replacement in dcb, erp and ionos file
*           2015/11/13  1.17 add support of L5 antenna phase center parameters
*                            add *.stec and *.trp file for ppp correction
*           2015/11/26  1.18 support opt->freqopt(disable L2)
*           2016/01/12  1.19 add carrier-phase bias correction by ssr
*           2016/07/31  1.20 fix error message problem in rnx2rtkp
*           2016/08/29  1.21 suppress warnings
*           2016/10/10  1.22 fix bug on identification of file fopt->blq
*           2017/06/13  1.23 add smoother of velocity solution
*           2020/11/30  1.24 use API sat2freq() to get carrier frequency
*                            fix bug on select best solution in static mode
*                            delete function to use L2 instead of L5 PCV
*                            writing solution file in binary mode
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#include "B2b.h"
#include "ppp_ar_passbypass.h"

#define MIN(x,y)    ((x)<(y)?(x):(y))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

#define MAXPRCDAYS  100          /* max days of continuous processing */
#define MAXINFILE   1000         /* max number of input files */
#define MAXINVALIDTM 100         /* max number of invalid time marks */

#define EPSILON 1e-3

/* constants/global variables ------------------------------------------------*/

static pcvs_t pcvss={0};        /* satellite antenna parameters */
static pcvs_t pcvsr={0};        /* receiver antenna parameters */
static obs_t obss={0};          /* observation data */
static nav_t navs={0};          /* navigation data */
static sbs_t sbss={0};          /* sbas messages */
static sta_t stas[MAXRCV];      /* station information */
static int nepoch=0;            /* number of observation epochs */
static int nitm  =0;            /* number of invalid time marks */
static int iobsu =0;            /* current rover observation data index */
static int iobsr =0;            /* current reference observation data index */
static int isbs  =0;            /* current sbas message index */
static int iitm  =0;            /* current invalid time mark index */
static int reverse=0;           /* analysis direction (0:forward,1:backward) */
static int aborts=0;            /* abort status */
static sol_t *solf;             /* forward solutions */
static sol_t *solb;             /* backward solutions */
static double *rbf;             /* forward base positions */
static double *rbb;             /* backward base positions */
static int isolf=0;             /* current forward solutions index */
static int isolb=0;             /* current backward solutions index */
static char proc_rov [64]="";   /* rover for current processing */
static char proc_base[64]="";   /* base station for current processing */
static char rtcm_file[1024]=""; /* rtcm data file */
static char rtcm_path[1024]=""; /* rtcm data path */
static gtime_t invalidtm[MAXINVALIDTM]={{0}};/* invalid time marks */
static rtcm_t rtcm;             /* rtcm control struct */
static FILE *fp_rtcm=NULL;      /* rtcm data file pointer */

static char B2b_files[MAXINFILE][1024]; /* B2b data files array */
static int n_B2b_files = 0;     /* number of B2b files */
static int cur_B2b_idx = 0;     /* current B2b file index */
static char B2b_path[1024]=""; /* B2b data path */
// static B2b_t B2b;              /* B2b control struct */
static raw_t B2braw;              /* raw control struct */

static FILE *fp_B2b=NULL;      /* B2b data file pointer */

/* Reset B2b reader state (used for Pass-by-Pass re-processing).
 * This prevents Pass2 from inheriting EOF state from Pass1 and
 * accidentally reusing stale SSR corrections (which can make solutions
 * appear only near the last SSR epoch).
 */
static void reset_B2b_reader(int day_idx)
{
    if (n_B2b_files <= 0) return;

    if (day_idx < 0) day_idx = 0;
    if (day_idx >= n_B2b_files) day_idx = n_B2b_files - 1;

    /* close any previously-open B2b file so the next update will reopen */
    if (fp_B2b) {
        fclose(fp_B2b);
        fp_B2b = NULL;
    }
    B2b_path[0] = '\0';
    cur_B2b_idx = day_idx;

    trace(3, "[B2B] reset reader to index %d/%d\n", cur_B2b_idx + 1, n_B2b_files);
}

static int processed_days = 0; /* counter for processed days */
static FILE *fp_URA=NULL;      /* B2b URA data file pointer */

/* ambiguity resolution data - defined globally for AR module access ---------*/
satamb_t satamb[MAXSAT] = {{0}}; /* ambiguity arcs for all satellites */
int n_ddamb = 0;                 /* number of DD ambiguities */
ddamb_t ddamb[MAXSAT*MAXSAT] = {{0}}; /* DD ambiguity data */
int refsat = 0;                  /* reference satellite for DD */
/* show message and check break ----------------------------------------------*/
static int checkbrk(const char *format, ...)
{
    va_list arg;
    char buff[1024],*p=buff;
    if (!*format) return showmsg("");
    va_start(arg,format);
    p+=vsprintf(p,format,arg);
    va_end(arg);
    if (*proc_rov&&*proc_base) sprintf(p," (%s-%s)",proc_rov,proc_base);
    else if (*proc_rov ) sprintf(p," (%s)",proc_rov );
    else if (*proc_base) sprintf(p," (%s)",proc_base);
    return showmsg(buff);
}
/* Solution option to field separator ----------------------------------------*/
/* Repeated from solution.c */
static const char *opt2sep(const solopt_t *opt)
{
    if (!*opt->sep) return " ";
    else if (!strcmp(opt->sep,"\\t")) return "\t";
    return opt->sep;
}
/* output reference position -------------------------------------------------*/
static void outrpos(FILE *fp, const double *r, const solopt_t *opt)
{
    double pos[3],dms1[3],dms2[3];

    trace(3,"outrpos :\n");

    const char *sep = opt2sep(opt);
    if (opt->posf==SOLF_LLH||opt->posf==SOLF_ENU) {
        ecef2pos(r,pos);
        if (opt->degf) {
            deg2dms(pos[0]*R2D,dms1,5);
            deg2dms(pos[1]*R2D,dms2,5);
            fprintf(fp,"%3.0f%s%02.0f%s%08.5f%s%4.0f%s%02.0f%s%08.5f%s%10.4f",
                    dms1[0],sep,dms1[1],sep,dms1[2],sep,dms2[0],sep,dms2[1],
                    sep,dms2[2],sep,pos[2]);
        }
        else {
            fprintf(fp,"%13.9f%s%14.9f%s%10.4f",pos[0]*R2D,sep,pos[1]*R2D,
                    sep,pos[2]);
        }
    }
    else if (opt->posf==SOLF_XYZ) {
        fprintf(fp,"%14.4f%s%14.4f%s%14.4f",r[0],sep,r[1],sep,r[2]);
    }
}
/* output header -------------------------------------------------------------*/
static void outheader(FILE *fp, const char **file, int n, const prcopt_t *popt,
                      const solopt_t *sopt)
{
    const char *s1[]={"GPST","UTC","JST"};
    gtime_t ts,te;
    double t1,t2;
    int i,j,w1,w2;
    char s2[32],s3[32];

    trace(3,"outheader: n=%d\n",n);

    if (sopt->posf==SOLF_NMEA||sopt->posf==SOLF_STAT) {
        return;
    }
    if (sopt->outhead) {
        if (!*sopt->prog) {
            fprintf(fp,"%s program   : RTKLIB ver.%s %s\n",COMMENTH,VER_RTKLIB,PATCH_LEVEL);
        }
        else {
            fprintf(fp,"%s program   : %s\n",COMMENTH,sopt->prog);
        }
        for (i=0;i<n;i++) {
            fprintf(fp,"%s inp file  : %s\n",COMMENTH,file[i]);
        }
        for (i=0;i<obss.n;i++)    if (obss.data[i].rcv==1) break;
        for (j=obss.n-1;j>=0;j--) if (obss.data[j].rcv==1) break;
        if (j<i) {fprintf(fp,"\n%s no rover obs data\n",COMMENTH); return;}
        ts=obss.data[i].time;
        te=obss.data[j].time;
        t1=time2gpst(ts,&w1);
        t2=time2gpst(te,&w2);
        if (sopt->times>=1) {
            ts=gpst2utc(ts);
            te=gpst2utc(te);
        }
        if (sopt->times==2) {
            ts=timeadd(ts,9*3600.0);
            te=timeadd(te,9*3600.0);
        }
        time2str(ts,s2,1);
        time2str(te,s3,1);
        fprintf(fp,"%s obs start : %s %s (week%04d %8.1fs)\n",COMMENTH,s2,s1[sopt->times],w1,t1);
        fprintf(fp,"%s obs end   : %s %s (week%04d %8.1fs)\n",COMMENTH,s3,s1[sopt->times],w2,t2);
    }
    if (sopt->outopt) {
        outprcopt(fp,popt);
    }
    if (PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED&&popt->mode!=PMODE_MOVEB) {
        fprintf(fp,"%s ref pos   :",COMMENTH);
        outrpos(fp,popt->rb,sopt);
        fprintf(fp,"\n");
    }
    if (sopt->outhead||sopt->outopt) fprintf(fp,"%s\n",COMMENTH);

    outsolhead(fp,sopt);
}
/* search next observation data index ----------------------------------------*/
static int nextobsf(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;

    for (;*i<obs->n;(*i)++) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);
        if (obs->data[*i+n].rcv!=rcv||tt>DTTOL) break;
    }
    return n;
}
static int nextobsb(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;

    for (;*i>=0;(*i)--) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i-n>=0;n++) {
        tt=timediff(obs->data[*i-n].time,obs->data[*i].time);
        if (obs->data[*i-n].rcv!=rcv||tt<-DTTOL) break;
    }
    return n;
}
/* update rtcm ssr correction ------------------------------------------------*/
static void update_rtcm_ssr(gtime_t time)
{
    char path[1024];
    int i;

    /* open or swap rtcm file */
    reppath(rtcm_file,path,time,"","");

    if (strcmp(path,rtcm_path)) {
        strcpy(rtcm_path,path);

        if (fp_rtcm) fclose(fp_rtcm);
        fp_rtcm=fopen(path,"rb");
        if (fp_rtcm) {
            rtcm.time=time;
            input_rtcm3f(&rtcm,fp_rtcm);
            trace(2,"rtcm file open: %s\n",path);
        }
    }
    if (!fp_rtcm) return;

    /* read rtcm file until current time */
    while (timediff(rtcm.time,time)<1E-3) {
        if (input_rtcm3f(&rtcm,fp_rtcm)<-1) break;

        /* update ssr corrections */
        for (i=0;i<MAXSAT;i++) {
            if (!rtcm.ssr[i].update||
                rtcm.ssr[i].iod[0]!=rtcm.ssr[i].iod[1]||
                timediff(time,rtcm.ssr[i].t0[0])<-1E-3) continue;
            navs.ssr[i]=rtcm.ssr[i];
            rtcm.ssr[i].update=0;
        }
    }
}

/* update B2b ssr correction ------------------------------------------------*/
static void update_B2b_ssr(gtime_t time, int format)
{
    char path[1024];
    char obstime_str[128] = {0};
    char B2btime_str[128] = {0};
    char t0_str[6][64];
    char satid[8];
    int i, j;
    int found = 0;

    if (format != STRFMT_SINO && format != STRFMT_UNICORE) {
        printf("Error Format: %d \n", format);
    }

    time2str(time, obstime_str, 3);
    trace(22, "update_B2b_ssr  : obstime=%s \n", obstime_str);
    trace(22, "---------------------------------------------------------\n");

    /* If current file is still open, continue using it - nothing to do here */
    if (fp_B2b != NULL) {
        /* File already open, skip to data reading section */
    } else {
        /* Need to open a new file - start from cur_B2b_idx (next file after EOF) */
        found = 0;
        for (j = cur_B2b_idx; j < n_B2b_files; j++) {
            reppath(B2b_files[j], path, time, "", "");

            /* Try to open this file */
            FILE *fp_test = fopen(path, "rb");
            if (fp_test) {
                fclose(fp_test);
                found = 1;
                cur_B2b_idx = j;
                trace(2, "Found B2b file [%d/%d]: %s\n", j+1, n_B2b_files, path);
                printf("Found B2b file [%d/%d]: %s\n", j+1, n_B2b_files, path);
                break;
            }
        }

        /* If no file found, return */
        if (!found) {
            trace(3, "No B2b file available for time: %s\n", obstime_str);
            return;
        }

        /* Open the found file */
        fp_B2b = fopen(path, "rb");
        if (!fp_B2b) {
            trace(1, "Failed to open B2b file: %s\n", path);
            printf("Error: Failed to open B2b file: %s\n", path);
            return;
        }

        /* Initialize for new file */
        strcpy(B2b_path, path);
        init_raw(&B2braw, format);
        B2braw.time = time;

        /* Read first data from file and check validity */
        int first_ret = input_rawf(&B2braw, format, fp_B2b);
        if (first_ret == -2 || first_ret == -1) {
            /* File is EOF or has error on first read - skip this file */
            trace(1, "ERROR: B2b file empty or immediate error: %s (ret=%d)\n", path, first_ret);
            printf("ERROR: B2b file empty or immediate error: %s (ret=%d)\n", path, first_ret);
            printf("  Skipping this file, will try next file on next call\n");
            fclose(fp_B2b);
            fp_B2b = NULL;
            B2b_path[0] = '\0';
            cur_B2b_idx++;  /* Skip to next file */
            return;
        }
        /* If first_ret == 0 (no complete message yet), continue - this is normal for binary format syncing */

        trace(2, "Opened B2b file [%d/%d]: %s, first B2b time=%s\n",
              cur_B2b_idx+1, n_B2b_files, path, time_str(B2braw.time,3));
        printf("Opened B2b file [%d/%d]: %s\n", cur_B2b_idx+1, n_B2b_files, path);
    }

    if (!fp_B2b) {
        trace(3, "No B2b file available for time: %s\n", obstime_str);
        return;
    }
    
    /* read B2b data until current time (assuming B2b data is directly available) */
    int max_iterations = 10000;  /* Reduced from 100000 to prevent long hangs */
    int iter_count = 0;
    gtime_t last_time = B2braw.time;
    int no_progress_count = 0;

    trace(3, "Starting B2b read loop: B2braw.time=%s, obs_time=%s\n",
          time_str(B2braw.time,3), obstime_str);

    while (timediff(B2braw.time, time) < 1E-3) {
        int ret;

        /* Check iteration limit to prevent infinite loop */
        if (++iter_count > max_iterations) {
            trace(1, "ERROR: B2b read loop exceeded %d iterations at %s\n",
                  max_iterations, obstime_str);
            printf("ERROR: B2b read loop exceeded %d iterations, stopping\n", max_iterations);
            printf("  B2braw.time=%s, obs_time=%s, file=%s\n",
                   time_str(B2braw.time,3), obstime_str, B2b_path);
            break;
        }

        /* Check if B2braw.time is progressing */
        if (timediff(B2braw.time, last_time) < 1E-6) {
            if (++no_progress_count > 1000) {  /* Reduced from 50000 */
                trace(1, "ERROR: B2b time not progressing after 1000 reads at %s\n", obstime_str);
                printf("ERROR: B2b time not progressing after 1000 reads\n");
                printf("  B2braw.time=%s stuck, obs_time=%s, file=%s\n",
                       time_str(B2braw.time,3), obstime_str, B2b_path);
                break;
            }
        } else {
            no_progress_count = 0;
            last_time = B2braw.time;
        }

        ret = input_rawf(&B2braw, format, fp_B2b);

        /* Handle EOF, error, or no data (ret <= 0) */
        if (ret <= 0) {
            if (ret == -2) {
                /* EOF - close file and move to next file index */
                trace(2, "B2b file EOF reached: %s (cur_B2b_idx=%d/%d)\n",
                      B2b_path, cur_B2b_idx, n_B2b_files);
                printf("B2b file EOF reached: %s\n", B2b_path);
                cur_B2b_idx++;
                trace(2, "Advanced to next B2b file index: %d/%d\n", cur_B2b_idx, n_B2b_files);
                printf("  Advanced to next file index [%d/%d] for next day\n",
                       cur_B2b_idx+1, n_B2b_files);
                if (fp_B2b) {
                    fclose(fp_B2b);
                    fp_B2b = NULL;
                }
                B2b_path[0] = '\0';
            } else if (ret == -1) {
                /* Error - close file and skip to next file */
                trace(1, "B2b file read error: %s\n", B2b_path);
                printf("B2b file read error: %s\n", B2b_path);
                cur_B2b_idx++;  /* Skip problematic file */
                if (fp_B2b) {
                    fclose(fp_B2b);
                    fp_B2b = NULL;
                }
                B2b_path[0] = '\0';
            } else if (ret == 0) {
                /* No complete message yet - keep file open for next call */
                trace(3, "B2b no complete message read, will continue next call\n");
                /* Do NOT close file or increment index - continue from here next time */
            }
            break;
        }
        time2str(B2braw.time, B2btime_str, 3);

        /* Liu@APM: This output may log a message where timediff(time, B2b.B2bssr[i].t0[0]) > -1E-3, so be cautious */
        if (B2braw.num_PPPB2BINF01 != 0) trace(22, "Message 1(%d) Detected at %s GeoPRN is %d \n", B2braw.num_PPPB2BINF01, B2btime_str, B2braw.geoprn);
        if (B2braw.num_PPPB2BINF02 != 0) trace(22, "Message 2(%d) Detected at %s GeoPRN is %d \n", B2braw.num_PPPB2BINF02, B2btime_str, B2braw.geoprn);
        if (B2braw.num_PPPB2BINF03 != 0) trace(22, "Message 3(%d) Detected at %s GeoPRN is %d \n", B2braw.num_PPPB2BINF03, B2btime_str, B2braw.geoprn);
        if (B2braw.num_PPPB2BINF04 != 0) trace(22, "Message 4(%d) Detected at %s GeoPRN is %d \n", B2braw.num_PPPB2BINF04, B2btime_str, B2braw.geoprn);
        trace(22, "---------------------------------------------------------\n");

        /* update ssr corrections */
        for (i = 0; i < MAXSAT; i++) {
            if (i==0) continue;
            satno2id(i, satid);
            if (!B2braw.nav.B2bssr[i].update) continue;
            if (B2braw.num_PPPB2BINF02 != 0) {
                if (timediff(time, B2braw.time) < -1E-3) continue;
                if (abs(timediff(B2braw.nav.B2bssr[i].t0[0], navs.B2bssr[i].t0[0])) < 1E-3) {
                    if (!checkout_B2beph(&B2braw.nav.B2bssr[i], &navs.B2bssr[i])) {
                        printf("error: Prn(%d) B2b_Eph is different from before one!", B2braw.geoprn);
                        trace(22, "error: Prn(%d) B2b_Eph is different from before one!", B2braw.geoprn);
                    }
                    B2braw.nav.B2bssr[i].update = 0; /* Consider identical data as already updated */
                    continue; /* Do not update identical data */
                }
                else B2braw.nav.B2bssr[i].udi[0] = timediff(B2braw.nav.B2bssr[i].t0[0], navs.B2bssr[i].t0[0]);

                if (B2braw.nav.B2bssr[i].iodcorr[0] != B2braw.nav.B2bssr[i].iodcorr[1]) {
                    satno2id(i, satid);
                    printf("warning: %s eph(%d) and clk(%d) iodcorr is different! \n", satid, B2braw.nav.B2bssr[i].iodcorr[0], B2braw.nav.B2bssr[i].iodcorr[1]);
                    trace(22, "warning: %s eph(%d) and clk(%d) iodcorr is different!\n", satid, B2braw.nav.B2bssr[i].iodcorr[0], B2braw.nav.B2bssr[i].iodcorr[1]);
                    B2braw.nav.B2bssr[i].update = 0;
                    continue;
                }
            } 

            if (B2braw.num_PPPB2BINF03 != 0) {
                if (timediff(time, B2braw.time) < -1E-3) continue;

                if (abs(timediff(B2braw.nav.B2bssr[i].t0[1], navs.B2bssr[i].t0[1])) < 1E-3) {
                    if (!checkout_B2bcbia(&B2braw.nav.B2bssr[i], &navs.B2bssr[i])) {
                        printf("error: Prn(%d) B2b_Cbia is different from before one!", B2braw.geoprn);
                        trace(22, "error: Prn(%d) B2b_Cbia is different from before one!", B2braw.geoprn);
                    }
                    B2braw.nav.B2bssr[i].update = 0; /* Consider identical data as already updated; note that if the IOD for type 3 remains the same, udi will stay 0 since cbias values rarely change */
                    continue; /* Do not update identical data */
                }
                else B2braw.nav.B2bssr[i].udi[1] = timediff(B2braw.nav.B2bssr[i].t0[1], navs.B2bssr[i].t0[1]);
                if (B2braw.nav.B2bssr[i].iodcorr[0] != B2braw.nav.B2bssr[i].iodcorr[1]) {
                    satno2id(i, satid);
                    printf("warning: %s eph(%d) and clk(%d) iodcorr is different! \n", satid, B2braw.nav.B2bssr[i].iodcorr[0], B2braw.nav.B2bssr[i].iodcorr[1]);
                    trace(22, "warning: %s eph(%d) and clk(%d) iodcorr is different!\n", satid, B2braw.nav.B2bssr[i].iodcorr[0], B2braw.nav.B2bssr[i].iodcorr[1]);
                    B2braw.nav.B2bssr[i].update = 0;
                    continue;
                    /* It might have been carried over by other correction messages */
                }
            }

            if (B2braw.num_PPPB2BINF04 != 0) {
                if (timediff(time, B2braw.time) < -1E-3) continue;

                if (abs(timediff(B2braw.nav.B2bssr[i].t0[2], navs.B2bssr[i].t0[2])) < 1E-3) {
                    if (!checkout_B2bclk(&B2braw.nav.B2bssr[i], &navs.B2bssr[i])) {
                        printf("error: Prn(%d) B2b_Clk is different from before one!", B2braw.geoprn);
                        trace(22, "error: Prn(%d) B2b_Clk is different from before one!", B2braw.geoprn);
                    }
                    B2braw.nav.B2bssr[i].update = 0; /* Consider identical data as already updated */
                    continue; /* Do not update identical data */
                }
                else B2braw.nav.B2bssr[i].udi[2] = timediff(B2braw.nav.B2bssr[i].t0[2], navs.B2bssr[i].t0[2]);
                if (B2braw.nav.B2bssr[i].iodcorr[0] != B2braw.nav.B2bssr[i].iodcorr[1]) {
                    satno2id(i, satid);
                    printf("warning: %s eph(%d) and clk(%d) iodcorr is different! \n", satid, B2braw.nav.B2bssr[i].iodcorr[0], B2braw.nav.B2bssr[i].iodcorr[1]);
                    trace(22, "warning: %s eph(%d) and clk(%d) iodcorr is different!\n", satid, B2braw.nav.B2bssr[i].iodcorr[0], B2braw.nav.B2bssr[i].iodcorr[1]);
                    B2braw.nav.B2bssr[i].update = 0;
                    continue;
                    /* It might have been carried over by other correction messages */
                }
            }
            
            navs.B2bssr[i] = B2braw.nav.B2bssr[i];

            /* Iterate through the t0 array, convert all its elements to strings, and output them */
            for (int j = 0; j < 6; j++) {
                time2str(navs.B2bssr[i].t0[j], t0_str[j], 3);
            }
            satno2id(i, satid);
            trace(22, "navs.B2bssr[%d] update: %d\n", i, navs.B2bssr[i].update);
            trace(22, "B2bssr[%d]: %s tod = %d, verify_sow = %d ,t0 = [%s, %s, %s, %s, %s, %s], udi = [%f, %f, %f, %f, %f, %f], "
                      "iodssr = [%d, %d, %d, %d, %d, %d], iodp = [%d, %d], iodn = %d, iodcorr = [%d, %d, %d, %d], "
                      "deph = [%f, %f, %f], ddeph = [%f, %f, %f], "
                      "ura = %d, cbias = [%f, %f, %f], dclk = [%f, %f, %f], update = %d\n",
                      i, satid, navs.B2bssr[i].sow, navs.B2bssr[i].verify_sow,
                      t0_str[0], t0_str[1], t0_str[2], t0_str[3], t0_str[4], t0_str[5],  
                      navs.B2bssr[i].udi[0], navs.B2bssr[i].udi[1], navs.B2bssr[i].udi[2], navs.B2bssr[i].udi[3],
                      navs.B2bssr[i].udi[4], navs.B2bssr[i].udi[5],
                      navs.B2bssr[i].iodssr[0], navs.B2bssr[i].iodssr[1], navs.B2bssr[i].iodssr[2],
                      navs.B2bssr[i].iodssr[3], navs.B2bssr[i].iodssr[4], navs.B2bssr[i].iodssr[5],
                      navs.B2bssr[i].iodp[0], navs.B2bssr[i].iodp[1], navs.B2bssr[i].iodn,
                      navs.B2bssr[i].iodcorr[0], navs.B2bssr[i].iodcorr[1], navs.B2bssr[i].iodcorr[2], navs.B2bssr[i].iodcorr[3],
                      navs.B2bssr[i].deph[0], navs.B2bssr[i].deph[1], navs.B2bssr[i].deph[2],
                      navs.B2bssr[i].ddeph[0], navs.B2bssr[i].ddeph[1], navs.B2bssr[i].ddeph[2],
                      navs.B2bssr[i].ura,
                      navs.B2bssr[i].cbias[0], navs.B2bssr[i].cbias[1], navs.B2bssr[i].cbias[2],
                      navs.B2bssr[i].dclk[0], navs.B2bssr[i].dclk[1], navs.B2bssr[i].dclk[2],
                      navs.B2bssr[i].update);

            trace(22, "---------------------------------------------------------\n");

            B2braw.nav.B2bssr[i].update = 0;
        }
        B2braw.num_PPPB2BINF01 = 0;
        B2braw.num_PPPB2BINF02 = 0;
        B2braw.num_PPPB2BINF03 = 0;
        B2braw.num_PPPB2BINF04 = 0;
    }
}



/* input obs data, navigation messages and sbas correction -------------------*/
static int inputobs(obsd_t *obs, int solq, const prcopt_t *popt)
{
    gtime_t time={0};
    int i,nu,nr,n=0;
    double dt,dt_next;
    int format;

    trace(3,"\ninfunc  : dir=%d iobsu=%d iobsr=%d isbs=%d\n",reverse,iobsu,iobsr,isbs);
    // traceobs_impl(22,obss.data,obss.n);

    if (0<=iobsu&&iobsu<obss.n) {
        settime((time=obss.data[iobsu].time));
        if (checkbrk("processing : %s Q=%d",time_str(time,0),solq)) {
            aborts=1; showmsg("aborted"); return -1;
        }
    }
    if (!reverse) { /* input forward data */
        if ((nu=nextobsf(&obss,&iobsu,1))<=0) return -1;
        if (popt->intpref) {
            /* for interpolation, find first base timestamp after rover timestamp */
            for (;(nr=nextobsf(&obss,&iobsr,2))>0;iobsr+=nr)
                if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)>-DTTOL) break;
        }
        else {
            /* if not interpolating, find closest timestamp */
            dt=timediff(obss.data[iobsr].time,obss.data[iobsu].time);
            for (i=iobsr;(nr=nextobsf(&obss,&i,2))>0;iobsr=i,i+=nr) {
                dt_next=timediff(obss.data[i].time,obss.data[iobsu].time);
                if (fabs(dt_next)>fabs(dt)) break;
                dt=dt_next;
            }
        }
        nr=nextobsf(&obss,&iobsr,2);
        if (nr<=0) {
            nr=nextobsf(&obss,&iobsr,2);
        }
        for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu+i];
        for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr+i];
        iobsu+=nu;

        /* update sbas corrections */
        while (isbs<sbss.n) {
            time=gpst2time(sbss.msgs[isbs].week,sbss.msgs[isbs].tow);

            if (getbitu(sbss.msgs[isbs].msg,8,6)!=9) { /* except for geo nav */
                sbsupdatecorr(sbss.msgs+isbs,&navs);
            }
            if (timediff(time,obs[0].time)>-1.0-DTTOL) break;
            isbs++;
        }
        /* update rtcm ssr corrections */
        if (*rtcm_file) {
            update_rtcm_ssr(obs[0].time);
        }

         if (n_B2b_files > 0) {
            // int format = STRFMT_SINAN;
            // strcpy(rtcm_file, "../../../../GNSS_DataSet/20241002APMALIC/sinan/Cor_202410020000.log");
            format = popt->B2b_format;
            update_B2b_ssr(obs[0].time,format);
        }

    }
    else { /* input backward data */
        if ((nu=nextobsb(&obss,&iobsu,1))<=0) return -1;
        if (popt->intpref) {
            /* for interpolation, find first base timestamp before rover timestamp */
            for (;(nr=nextobsb(&obss,&iobsr,2))>0;iobsr-=nr)
                if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)<DTTOL) break;
        }
        else {
            /* if not interpolating, find closest timestamp */
            dt=iobsr>=0?timediff(obss.data[iobsr].time,obss.data[iobsu].time):0;
            for (i=iobsr;(nr=nextobsb(&obss,&i,2))>0;iobsr=i,i-=nr) {
                dt_next=timediff(obss.data[i].time,obss.data[iobsu].time);
                if (fabs(dt_next)>fabs(dt)) break;
                dt=dt_next;
            }
        }
        nr=nextobsb(&obss,&iobsr,2);
        for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu-nu+1+i];
        for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr-nr+1+i];
        iobsu-=nu;

        /* update sbas corrections */
        while (isbs>=0) {
            time=gpst2time(sbss.msgs[isbs].week,sbss.msgs[isbs].tow);

            if (getbitu(sbss.msgs[isbs].msg,8,6)!=9) { /* except for geo nav */
                sbsupdatecorr(sbss.msgs+isbs,&navs);
            }
            if (timediff(time,obs[0].time)<1.0+DTTOL) break;
            isbs--;
        }
    }
    return n;
}
/* output to file message of invalid time mark -------------------------------*/
static void outinvalidtm(FILE *fptm, const solopt_t *opt, const gtime_t tm)
{
    gtime_t time = tm;
    double gpst;
    const double secondsInAWeek = 604800;
    int week,timeu;
    char s[100];

    timeu=opt->timeu<0?0:(opt->timeu>20?20:opt->timeu);

    if (opt->times>=TIMES_UTC) time=gpst2utc(time);
    if (opt->times==TIMES_JST) time=timeadd(time,9*3600.0);

    if (opt->timef) time2str(time,s,timeu);
    else {
        gpst=time2gpst(time,&week);
        if (secondsInAWeek-gpst < 0.5/pow(10.0,timeu)) {
            week++;
            gpst=0.0;
        }
        sprintf(s,"%4d   %*.*f",week,6+(timeu<=0?0:timeu+1),timeu,gpst);
    }
    strcat(s, "   Q=0, Time mark is not valid\n");

    fwrite(s,strlen(s),1,fptm);
}
/* fill structure sol_t for time mark ----------------------------------------*/
static sol_t fillsoltm(const sol_t solold, const sol_t solnew, const gtime_t tm)
{
    gtime_t t1={0},t2={0};
    sol_t sol=solold;
    int i=0;

    if (solold.stat == 0 || solnew.stat == 0) {
        sol.stat = 0;
    } else {
        sol.stat = (solold.stat > solnew.stat) ? solold.stat : solnew.stat;
    }
    sol.ns = (solold.ns < solnew.ns) ? solold.ns : solnew.ns;
    sol.ratio = (solold.ratio < solnew.ratio) ? solold.ratio : solnew.ratio;

    /* interpolation position and speed of time mark */
    t1 = solold.time;
    t2 = solnew.time;
    sol.time = tm;

    for (i=0;i<6;i++)
    {
        sol.rr[i] = solold.rr[i] + timediff(tm,t1) / timediff(t2,t1) * (solnew.rr[i] - solold.rr[i]);
    }

    return sol;
}

/* carrier-phase bias correction by ssr --------------------------------------*/
static void corr_phase_bias_ssr(obsd_t *obs, int n, const nav_t *nav)
{
    double freq;
    uint8_t code;
    int i,j;

    for (i=0;i<n;i++) for (j=0;j<NFREQ;j++) {
        code=obs[i].code[j];

        if ((freq=sat2freq(obs[i].sat,code,nav))==0.0) continue;

        /* correct phase bias (cyc) */
        obs[i].L[j]-=nav->ssr[obs[i].sat-1].pbias[code-1]*freq/CLIGHT;
    }
}
/* process positioning -------------------------------------------------------*/
static void procpos(FILE *fp, FILE *fptm, const prcopt_t *popt, const solopt_t *sopt,
                    rtk_t *rtk, int mode)
{
    gtime_t time={0};
    sol_t sol={{0}},oldsol={{0}},newsol={{0}};
    obsd_t *obs_ptr = (obsd_t *)malloc(sizeof(obsd_t)*MAXOBS*2); /* for rover and base */
    double rb[3]={0};
    int i,nobs,n,solstatic,num=0,pri[]={6,1,2,3,4,5,1,6};
    static gtime_t day1_start = {0};  /* first day start time for AR */
    static int current_day = -1;      /* current processing day for AR */

    trace(3,"procpos : mode=%d\n",mode); /* 0=single dir, 1=combined */
    char ura_filename[256];
    sprintf(ura_filename, "./b2b_ura_data_%s.txt", popt->sationname);

    // if (!fp_URA) {
    //     fp_URA = fopen(ura_filename, "w+");
    //     if (!fp_URA) {
    //         fprintf(stderr, "Error: Failed to open B2b URA data file!\n");
    //         return;
    //     }
    //     trace(2, "B2b URA data file opened successfully.\n");
    // }

    solstatic=sopt->solstatic&&
              (popt->mode==PMODE_STATIC||popt->mode==PMODE_STATIC_START||popt->mode==PMODE_PPP_STATIC);
    
    rtcm_path[0]='\0';

    while ((nobs=inputobs(obs_ptr,rtk->sol.stat,popt))>=0) {

        /* exclude satellites */
        for (i=n=0;i<nobs;i++) {
            if ((satsys(obs_ptr[i].sat,NULL)&popt->navsys)&&
                popt->exsats[obs_ptr[i].sat-1]!=1) obs_ptr[n++]= obs_ptr[i];
        }
        if (n<=0) continue;

        /* carrier-phase bias correction */
        // if (!strstr(popt->pppopt,"-ENA_FCB")) {
        //     corr_phase_bias_ssr(obs_ptr,n,&navs);
        // }
        if (!rtkpos(rtk, obs_ptr,n,&navs)) {
            if (rtk->sol.eventime.time != 0) {
                if (mode == SOLMODE_SINGLE_DIR) {
                    outinvalidtm(fptm, sopt, rtk->sol.eventime);
                } else if (!reverse&&nitm<MAXINVALIDTM) {
                    invalidtm[nitm++] = rtk->sol.eventime;
                }
            }
            continue;
        }

         
/* Pass-by-Pass AR: collection only during Pass-1.
 * OLD CODE REMOVED HERE:
 *   - pbp_apply_flag driven epoch-wise apply_ar_fixed() update
 *   - all dynamic Kalman constraint injection inside procpos()
 */
{
    extern int pbp_day_tag, pbp_collect_flag;

    if (pbp_collect_flag) {
        int cur_day = pbp_day_tag;

        if (cur_day < 0) {
            extern int pbp_base_day_id; /* defined in ppp_ar_passbypass.c */
            int week = 0;
            double sow = time2gpst(obs_ptr[0].time, &week);
            int day_id = week * 7 + (int)floor(sow / 86400.0);
            if (pbp_base_day_id < 0) pbp_base_day_id = day_id;
            cur_day = day_id - pbp_base_day_id;
        }

        if ((cur_day == 0 || cur_day == 1) && rtk->sol.stat == SOLQ_PPP) {
            (void)collect_ambiguities_epoch(rtk, obs_ptr, n, cur_day);
        }
    }
}

if (mode==SOLMODE_SINGLE_DIR) {    /* forward or backward */
            if (!solstatic) {
                outsol(fp,&rtk->sol,rtk->rb,sopt);
            }
            else if (time.time==0||pri[rtk->sol.stat]<=pri[sol.stat]) {
                sol=rtk->sol;
                for (i=0;i<3;i++) rb[i]=rtk->rb[i];
                if (time.time==0||timediff(rtk->sol.time,time)<0.0) {
                    time=rtk->sol.time;
                }
            }
            /* check time mark */
            if (rtk->sol.eventime.time != 0)
            {
                newsol = fillsoltm(oldsol,rtk->sol,rtk->sol.eventime);
                num++;
                if (!solstatic&&mode==SOLMODE_SINGLE_DIR) {
                    outsol(fptm,&newsol,rb,sopt);
                }
            }
            oldsol = rtk->sol;
        }
        else if (!reverse) { /* combined-forward */
            if (isolf >= nepoch) {
                free(obs_ptr);
                return;
            }
            solf[isolf]=rtk->sol;
            for (i=0;i<3;i++) rbf[i+isolf*3]=rtk->rb[i];
            isolf++;
        }
        else { /* combined-backward */
            if (isolb>=nepoch) {
                free(obs_ptr);
                return;
            }
            solb[isolb]=rtk->sol;
            for (i=0;i<3;i++) rbb[i+isolb*3]=rtk->rb[i];
            isolb++;
        }
    }
    if (mode==SOLMODE_SINGLE_DIR && solstatic&&time.time!=0.0) {
        sol.time=time;
        outsol(fp,&sol,rb,sopt);
    }

    free(obs_ptr); /* moved from stack to heap to kill a stack overflow warning */
}
/* validation of combined solutions ------------------------------------------*/
static int valcomb(const sol_t *solf, const sol_t *solb, double *rbf,
        double *rbb, const prcopt_t *popt)
{
    double dr[3],var[3];
    int i;
    char tstr[32];

    trace(4,"valcomb :\n");

    /* compare forward and backward solution */
    for (i=0;i<3;i++) {
        dr[i]=solf->rr[i]-solb->rr[i];
        if (popt->mode==PMODE_MOVEB) dr[i]-=(rbf[i]-rbb[i]);
        var[i]=(double)solf->qr[i] + (double)solb->qr[i];
    }
    for (i=0;i<3;i++) {
        if (dr[i]*dr[i]<=16.0*var[i]) continue; /* ok if in 4-sigma */

        time2str(solf->time,tstr,2);
        trace(2,"degrade fix to float: %s dr=%.3f %.3f %.3f std=%.3f %.3f %.3f\n",
              tstr+11,dr[0],dr[1],dr[2],SQRT(var[0]),SQRT(var[1]),SQRT(var[2]));
        return 0;
    }
    return 1;
}
/* combine forward/backward solutions and save results ---------------------*/
static void combres(FILE *fp, FILE *fptm, const prcopt_t *popt, const solopt_t *sopt)
{
    gtime_t time={0};
    sol_t sols={{0}},sol={{0}},oldsol={{0}},newsol={{0}};
    double tt,Qf[9],Qb[9],Qs[9],rbs[3]={0},rb[3]={0},rr_f[3],rr_b[3],rr_s[3];
    int i,j,k,solstatic,num=0,pri[]={7,1,2,3,4,5,1,6};

    trace(3,"combres : isolf=%d isolb=%d\n",isolf,isolb);

    solstatic=sopt->solstatic&&
              (popt->mode==PMODE_STATIC||popt->mode==PMODE_STATIC_START||popt->mode==PMODE_PPP_STATIC);

    for (i=0,j=isolb-1;i<isolf&&j>=0;i++,j--) {
        if ((tt=timediff(solf[i].time,solb[j].time))<-DTTOL) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
            j++;
        }
        else if (tt>DTTOL) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
            i--;
        }
        else if (pri[solf[i].stat]<pri[solb[j].stat]) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
        }
        else if (pri[solf[i].stat]>pri[solb[j].stat]) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
        }
        else {
            sols=solf[i];
            sols.time=timeadd(sols.time,-tt/2.0);

            if ((popt->mode==PMODE_KINEMA||popt->mode==PMODE_MOVEB)&&
                sols.stat==SOLQ_FIX) {

                /* degrade fix to float if validation failed */
                if (!valcomb(solf+i,solb+j,rbf+i*3,rbb+j*3,popt)) sols.stat=SOLQ_FLOAT;
            }
            for (k=0;k<3;k++) {
                Qf[k+k*3]=solf[i].qr[k];
                Qb[k+k*3]=solb[j].qr[k];
            }
            Qf[1]=Qf[3]=solf[i].qr[3];
            Qf[5]=Qf[7]=solf[i].qr[4];
            Qf[2]=Qf[6]=solf[i].qr[5];
            Qb[1]=Qb[3]=solb[j].qr[3];
            Qb[5]=Qb[7]=solb[j].qr[4];
            Qb[2]=Qb[6]=solb[j].qr[5];

            if (popt->mode==PMODE_MOVEB) {
                for (k=0;k<3;k++) rr_f[k]=solf[i].rr[k]-rbf[k+i*3];
                for (k=0;k<3;k++) rr_b[k]=solb[j].rr[k]-rbb[k+j*3];
                if (smoother(rr_f,Qf,rr_b,Qb,3,rr_s,Qs)) continue;
                for (k=0;k<3;k++) sols.rr[k]=rbs[k]+rr_s[k];
            }
            else {
                if (smoother(solf[i].rr,Qf,solb[j].rr,Qb,3,sols.rr,Qs)) continue;
            }
            sols.qr[0]=(float)Qs[0];
            sols.qr[1]=(float)Qs[4];
            sols.qr[2]=(float)Qs[8];
            sols.qr[3]=(float)Qs[1];
            sols.qr[4]=(float)Qs[5];
            sols.qr[5]=(float)Qs[2];

            /* smoother for velocity solution */
            if (popt->dynamics) {
                for (k=0;k<3;k++) {
                    Qf[k+k*3]=solf[i].qv[k];
                    Qb[k+k*3]=solb[j].qv[k];
                }
                Qf[1]=Qf[3]=solf[i].qv[3];
                Qf[5]=Qf[7]=solf[i].qv[4];
                Qf[2]=Qf[6]=solf[i].qv[5];
                Qb[1]=Qb[3]=solb[j].qv[3];
                Qb[5]=Qb[7]=solb[j].qv[4];
                Qb[2]=Qb[6]=solb[j].qv[5];
                if (smoother(solf[i].rr+3,Qf,solb[j].rr+3,Qb,3,sols.rr+3,Qs)) continue;
                sols.qv[0]=(float)Qs[0];
                sols.qv[1]=(float)Qs[4];
                sols.qv[2]=(float)Qs[8];
                sols.qv[3]=(float)Qs[1];
                sols.qv[4]=(float)Qs[5];
                sols.qv[5]=(float)Qs[2];
            }
        }
        if (!solstatic) {
            outsol(fp,&sols,rbs,sopt);
        }
        else if (time.time==0||pri[sols.stat]<=pri[sol.stat]) {
            sol=sols;
            for (k=0;k<3;k++) rb[k]=rbs[k];
            if (time.time==0||timediff(sols.time,time)<0.0) {
                time=sols.time;
            }
        }
        if (iitm < nitm && timediff(invalidtm[iitm],sols.time)<0.0)
        {
            outinvalidtm(fptm,sopt,invalidtm[iitm]);
            iitm++;
        }
        if (sols.eventime.time != 0)
        {
            newsol = fillsoltm(oldsol,sols,sols.eventime);
            num++;
            if (!solstatic) {
                outsol(fptm,&newsol,rb,sopt);
            }
        }
        oldsol = sols;
    }
    if (solstatic&&time.time!=0.0) {
        sol.time=time;
        outsol(fp,&sol,rb,sopt);
    }
}
/* read prec ephemeris, sbas data, tec grid and open rtcm --------------------*/
static void readpreceph(const char **infile, int n, const prcopt_t *prcopt,
                        nav_t *nav, sbs_t *sbs)
{
    seph_t seph0={0};
    int i;
    const char *ext;

    trace(2,"readpreceph: n=%d\n",n);

    nav->ne=nav->nemax=0;
    nav->nc=nav->ncmax=0;
    sbs->n =sbs->nmax =0;

    /* read precise ephemeris files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        readsp3(infile[i],nav,0);
    }
    /* read precise clock files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        readrnxc(infile[i],nav);
    }
    /* read sbas message files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        sbsreadmsg(infile[i],prcopt->sbassatsel,sbs);
    }
    /* allocate sbas ephemeris */
    nav->ns=nav->nsmax=NSATSBS*2;
    if (!(nav->seph=(seph_t *)malloc(sizeof(seph_t)*nav->ns))) {
         showmsg("error : sbas ephem memory allocation");
         trace(1,"error : sbas ephem memory allocation");
         return;
    }
    for (i=0;i<nav->ns;i++) nav->seph[i]=seph0;

    /* set rtcm file and initialize rtcm struct */
    rtcm_file[0]=rtcm_path[0]='\0'; fp_rtcm=NULL;

    for (i=0;i<n;i++) {
        if ((ext=strrchr(infile[i],'.'))&&
            (!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) {
            strcpy(rtcm_file,infile[i]);
            init_rtcm(&rtcm);
            break;
        }
    }

       /* set B2b files and initialize B2b struct */
    /* Only initialize B2b files on first call (when n_B2b_files is 0) */
    /* This preserves cur_B2b_idx across multiple days of processing */
    if (n_B2b_files == 0) {
        cur_B2b_idx = 0;
        B2b_path[0]='\0';
        fp_B2b=NULL;

        /* Store all B2b files found in input */
        for (i=0;i<n && n_B2b_files < MAXINFILE;i++) {
            if ((ext=strrchr(infile[i],'.'))&&
                (!strcmp(ext,".b2b")||!strcmp(ext,".B2b"))) {
                strcpy(B2b_files[n_B2b_files],infile[i]);
                n_B2b_files++;
                printf("B2b file %d: %s\n", n_B2b_files, infile[i]);
            }
        }

        /* Initialize B2b raw struct if we have at least one B2b file */
        if (n_B2b_files > 0) {
            init_raw(&B2braw,prcopt->B2b_format);
            trace(2, "Initialized B2braw, Total B2b files found: %d\n", n_B2b_files);
            printf("Total B2b files found: %d\n", n_B2b_files);
        }
    } else {
        /* B2b files already initialized, keep existing cur_B2b_idx */
        /* BUT: need to re-initialize B2braw after free_raw() in freepreceph() */
        if (n_B2b_files > 0) {
            init_raw(&B2braw,prcopt->B2b_format);
            trace(2, "Re-initialized B2braw for new day processing\n");
        }
        trace(2, "B2b files already initialized, cur_B2b_idx=%d/%d\n", cur_B2b_idx, n_B2b_files);
        printf("Day %d: Continuing with B2b file index [%d/%d] (next file after previous EOF)\n",
               processed_days+1, cur_B2b_idx+1, n_B2b_files);
    }

}
/* free prec ephemeris and sbas data -----------------------------------------*/
static void freepreceph(nav_t *nav, sbs_t *sbs)
{
    int i;

    trace(3,"freepreceph:\n");

    free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
    free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
    free(sbs->msgs); sbs->msgs=NULL; sbs->n =sbs->nmax =0;
    for (i=0;i<nav->nt;i++) {
        free(nav->tec[i].data);
        free(nav->tec[i].rms );
    }
    free(nav->tec ); nav->tec =NULL; nav->nt=nav->ntmax=0;

    if (fp_rtcm) fclose(fp_rtcm);
    free_rtcm(&rtcm);
    /* close B2b file pointer */
    if (fp_B2b) {
        fclose(fp_B2b);
        fp_B2b = NULL;
    }
    B2b_path[0] = '\0';
    free_raw(&B2braw);
}
/* read obs and nav data -----------------------------------------------------*/
static int readobsnav(gtime_t ts, gtime_t te, double ti, const char **infile,
                      const int *index, int n, const prcopt_t *prcopt,
                      obs_t *obs, nav_t *nav, sta_t *sta)
{
    int i,j,ind=0,nobs=0,rcv=1;

    trace(3,"readobsnav: ts=%s n=%d\n",time_str(ts,0),n);

    obs->data=NULL; obs->n =obs->nmax =0;
    nav->eph =NULL; nav->n =nav->nmax =0;
    nav->geph=NULL; nav->ng=nav->ngmax=0;
    /* free(nav->seph); */ /* is this needed to avoid memory leak??? */
    nav->seph=NULL; nav->ns=nav->nsmax=0;
    nepoch=0;

    for (i=0;i<n;i++) {
        if (checkbrk("")) return 0;

        if (index[i]!=ind) {
            if (obs->n>nobs) rcv++;
            ind=index[i]; nobs=obs->n;
        }
        /* read rinex obs and nav file */
        if (readrnxt(infile[i],rcv,ts,te,ti,prcopt->rnxopt[rcv<=1?0:1],obs,nav,
                     rcv<=2?sta+rcv-1:NULL)<0) {
            checkbrk("error : insufficient memory");
            trace(1,"insufficient memory\n");
            return 0;
        }
    }
    if (obs->n<=0) {
        checkbrk("error : no obs data");
        trace(1,"\n");
        return 0;
    }
    if (nav->n<=0&&nav->ng<=0&&nav->ns<=0) {
        checkbrk("error : no nav data");
        trace(1,"\n");
        return 0;
    }
    /* sort observation data */
    nepoch=sortobs(obs);

    /* delete duplicated ephemeris */
    uniqnav(nav);

    /* set time span for progress display */
    if (ts.time==0||te.time==0) {
        for (i=0;   i<obs->n;i++) if (obs->data[i].rcv==1) break;
        for (j=obs->n-1;j>=0;j--) if (obs->data[j].rcv==1) break;
        if (i<j) {
            if (ts.time==0) ts=obs->data[i].time;
            if (te.time==0) te=obs->data[j].time;
            settspan(ts,te);
        }
    }
    return 1;
}
/* free obs and nav data -----------------------------------------------------*/
static void freeobsnav(obs_t *obs, nav_t *nav)
{
    trace(3,"freeobsnav:\n");

    free(obs->data); obs->data=NULL; obs->n =obs->nmax =0;
    free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;
    free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
}
/* average of single position ------------------------------------------------*/
static int avepos(double *ra, int rcv, const obs_t *obs, const nav_t *nav,
                  const prcopt_t *opt)
{
    obsd_t data[MAXOBS];
    gtime_t ts={0};
    sol_t sol={{0}};
    int i,j,n=0,m,iobs;
    char msg[128];

    trace(3,"avepos: rcv=%d obs.n=%d\n",rcv,obs->n);

    for (i=0;i<3;i++) ra[i]=0.0;

    for (iobs=0;(m=nextobsf(obs,&iobs,rcv))>0;iobs+=m) {

        for (i=j=0;i<m&&i<MAXOBS;i++) {
            data[j]=obs->data[iobs+i];
            if ((satsys(data[j].sat,NULL)&opt->navsys)&&
                opt->exsats[data[j].sat-1]!=1) j++;
        }
        if (j<=0||!screent(data[0].time,ts,ts,1.0)) continue; /* only 1 hz */

        if (!pntpos(data,j,nav,opt,&sol,NULL,NULL,msg)) continue;

        for (i=0;i<3;i++) ra[i]+=sol.rr[i];
        n++;
    }
    if (n<=0) {
        trace(1,"no average of base station position\n");
        return 0;
    }
    for (i=0;i<3;i++) ra[i]/=n;
    return 1;
}
/* station position from file ------------------------------------------------*/
static int getstapos(const char *file, const char *name, double *r)
{
    FILE *fp;
    char buff[256],sname[256],*p;
    const char *q;
    double pos[3];

    trace(3,"getstapos: file=%s name=%s\n",file,name);

    if (!(fp=fopen(file,"r"))) {
        trace(1,"station position file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if ((p=strchr(buff,'%'))) *p='\0';
        
        if (sscanf(buff,"%lf %lf %lf %255s",pos,pos+1,pos+2,sname)<4) continue;
        
        for (p=sname,q=name;*p&&*q;p++,q++) {
            if (toupper((int)*p)!=toupper((int)*q)) break;
        }
        if (!*p) {
            pos[0]*=D2R;
            pos[1]*=D2R;
            pos2ecef(pos,r);
            fclose(fp);
            return 1;
        }
    }
    fclose(fp);
    trace(1,"no station position: %s %s\n",name,file);
    return 0;
}
/* antenna phase center position ---------------------------------------------*/
static int antpos(prcopt_t *opt, int rcvno, const obs_t *obs, const nav_t *nav,
                  const sta_t *sta, const char *posfile)
{
    double *rr=rcvno==1?opt->ru:opt->rb,del[3],pos[3],dr[3]={0};
    int i,postype=rcvno==1?opt->rovpos:opt->refpos;
    char *name;

    trace(3,"antpos  : rcvno=%d\n",rcvno);

    if (postype==POSOPT_SINGLE) { /* average of single position */
        if (!avepos(rr,rcvno,obs,nav,opt)) {
            showmsg("error : station pos computation");
            return 0;
        }
    }
    else if (postype==POSOPT_FILE) { /* read from position file */
        name=stas[rcvno==1?0:1].name;
        if (!getstapos(posfile,name,rr)) {
            showmsg("error : no position of %s in %s",name,posfile);
            return 0;
        }
    }
    else if (postype==POSOPT_RINEX) { /* get from rinex header */
        if (norm(stas[rcvno==1?0:1].pos,3)<=0.0) {
            showmsg("error : no position in rinex header");
            trace(1,"no position in rinex header\n");
            return 0;
        }
        /* add antenna delta unless already done in antpcv() */
        if (!strcmp(opt->anttype[rcvno],"*")) {
            if (stas[rcvno==1?0:1].deltype==0) { /* enu */
                for (i=0;i<3;i++) del[i]=stas[rcvno==1?0:1].del[i];
                del[2]+=stas[rcvno==1?0:1].hgt;
                ecef2pos(stas[rcvno==1?0:1].pos,pos);
                enu2ecef(pos,del,dr);
            }  else { /* xyz */
                for (i=0;i<3;i++) dr[i]=stas[rcvno==1?0:1].del[i];
            }
        }
        for (i=0;i<3;i++) rr[i]=stas[rcvno==1?0:1].pos[i]+dr[i];
    }
    return 1;
}
/* open processing session ----------------------------------------------------*/
static int openses(const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    trace(3,"openses :\n");

    /* read satellite antenna parameters */
    if (*fopt->satantp&&!(readpcv(fopt->satantp,pcvs))) {
        showmsg("error : no sat ant pcv in %s",fopt->satantp);
        trace(1,"sat antenna pcv read error: %s\n",fopt->satantp);
        return 0;
    }
    /* read receiver antenna parameters */
    if (*fopt->rcvantp&&!(readpcv(fopt->rcvantp,pcvr))) {
        showmsg("error : no rec ant pcv in %s",fopt->rcvantp);
        trace(1,"rec antenna pcv read error: %s\n",fopt->rcvantp);
        return 0;
    }
    /* open geoid data */
    if (sopt->geoid>0&&*fopt->geoid) {
        if (!opengeoid(sopt->geoid,fopt->geoid)) {
            showmsg("error : no geoid data %s",fopt->geoid);
            trace(2,"no geoid data %s\n",fopt->geoid);
        }
    }
    return 1;
}
/* close processing session ---------------------------------------------------*/
static void closeses(nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    trace(3,"closeses:\n");

    /* free antenna parameters */
    free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
    free(pcvr->pcv); pcvr->pcv=NULL; pcvr->n=pcvr->nmax=0;

    /* close geoid data */
    closegeoid();

    /* free erp data */
    free(nav->erp.data); nav->erp.data=NULL; nav->erp.n=nav->erp.nmax=0;

    /* close solution statistics and debug trace */
    rtkclosestat();
    traceclose();
    B2b_traceclose();
}
/* set antenna parameters ----------------------------------------------------*/
static void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs,
                   const pcvs_t *pcvr, const sta_t *sta)
{
    pcv_t *pcv,pcv0={0};
    double pos[3],del[3];
    int i,j,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    char id[64];

    /* set satellite antenna parameters */
    for (i=0;i<MAXSAT;i++) {
        nav->pcvs[i]=pcv0;
        if (!(satsys(i+1,NULL)&popt->navsys)) continue;
        if (!(pcv=searchpcv(i+1,"",time,pcvs))) {
            satno2id(i+1,id);
            trace(4,"no satellite antenna pcv: %s\n",id);
            continue;
        }
        nav->pcvs[i]=*pcv;
    }
    for (i=0;i<(mode?2:1);i++) {
        popt->pcvr[i]=pcv0;
        if (!strcmp(popt->anttype[i],"*")) { /* set by station parameters */
            strcpy(popt->anttype[i],sta[i].antdes);
            if (sta[i].deltype==1) { /* xyz */
                if (norm(sta[i].pos,3)>0.0) {
                    ecef2pos(sta[i].pos,pos);
                    ecef2enu(pos,sta[i].del,del);
                    for (j=0;j<3;j++) popt->antdel[i][j]=del[j];
                }
            }
            else { /* enu */
                for (j=0;j<3;j++) popt->antdel[i][j]=stas[i].del[j];
            }
        }
        if (!(pcv=searchpcv(0,popt->anttype[i],time,pcvr))) {
            trace(2,"no receiver antenna pcv: %s\n",popt->anttype[i]);
            *popt->anttype[i]='\0';
            continue;
        }
        strcpy(popt->anttype[i],pcv->type);
        popt->pcvr[i]=*pcv;
    }
}
/* read ocean tide loading parameters ----------------------------------------*/
static void readotl(prcopt_t *popt, const char *file, const sta_t *sta)
{
    int i,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;

    for (i=0;i<(mode?2:1);i++) {
        readblq(file,sta[i].name,popt->odisp[i]);
    }
}
/* write header to output file -----------------------------------------------*/
static int outhead(const char *outfile, const char **infile, int n,
                   const prcopt_t *popt, const solopt_t *sopt)
{
    FILE *fp=stdout;
    int need_header=1;

    trace(3,"outhead: outfile=%s n=%d\n",outfile,n);

    if (*outfile) {
        createdir(outfile);

        /*
         * NOTE (patched): some drivers call outhead() multiple times for the same
         * output file when splitting processing by time unit. The original code
         * opened the file with "wb" unconditionally, which truncates the file and
         * can leave only the last segment in .pos/.stat outputs.
         *
         * Here we only write the header if the file is empty; otherwise we append.
         */
        {
            FILE *ft=fopen(outfile,"rb");
            if (ft) {
                fseek(ft,0,SEEK_END);
                if (ftell(ft)>0) need_header=0;
                fclose(ft);
            }
        }

        if (!(fp=fopen(outfile,need_header?"wb":"ab"))) {
            showmsg("error : open output file %s",outfile);
            return 0;
        }
    }

    /* output header (only when creating a new/empty file) */
    if (need_header) outheader(fp,infile,n,popt,sopt);

    if (*outfile) fclose(fp);

    return 1;
}
/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *outfile)
{
    trace(3,"openfile: outfile=%s\n",outfile);

    return !*outfile?stdout:fopen(outfile,"ab");
}
/* Name time marks file ------------------------------------------------------*/
static void namefiletm(char *outfiletm, const char *outfile)
{
    int i;

    for (i=(int)strlen(outfile);i>0;i--) {
        if (outfile[i] == '.') {
            break;
        }
    }
    /* if no file extension, then name time marks file as name of outfile + _events.pos */
    if (i == 0) {
        i = (int)strlen(outfile);
    }
    strncpy(outfiletm, outfile, i);
    strcat(outfiletm, "_events.pos");
}
/* execute processing session ------------------------------------------------*/
static int execses(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                   const solopt_t *sopt, const filopt_t *fopt, int flag,
                   const char **infile, const int *index, int n, const char *outfile)
{
    rtk_t *rtk_ptr = (rtk_t *)malloc(sizeof(rtk_t)); /* moved from stack to heap to avoid stack overflow warning */
    prcopt_t popt_=*popt;
    char B2btracefile[1024], tracefile[1024],statfile[1024],path[1024],outfiletm[1024]={0};
    const char *ext;
    int i,j,k,dcb_ok;

    trace(3,"execses : n=%d outfile=%s\n",n,outfile);

    /* open debug trace */
    if (flag&&sopt->trace>0) {
        if (*outfile) {
            strcpy(tracefile,outfile);
            strcat(tracefile,".trace");
        }
        else {
            strcpy(tracefile,fopt->trace);
        }
        traceclose();
        traceopen(tracefile);
        tracelevel(sopt->trace);
    }

    /* open B2b SSR/OSB trace (independent of main trace level) */
    if (flag && outfile && *outfile) {
        int b2b_trlev = sopt->trace;
        if (b2b_trlev < 22) b2b_trlev = 22;
        B2b_traceclose();
        strcpy(B2btracefile,outfile);
        strcat(B2btracefile,".B2bssr");
        B2b_tracelevel(b2b_trlev);
        B2b_traceopen(B2btracefile);
    }

/* read ionosphere data file */
    if (*fopt->iono&&(ext=strrchr(fopt->iono,'.'))) {
        if (strlen(ext)==4&&(ext[3]=='i'||ext[3]=='I'||
                             strcmp(ext,".INX")==0||strcmp(ext,".inx")==0)) {
            reppath(fopt->iono,path,ts,"","");
            readtec(path,&navs,1);
        }
    }
    /* read erp data */
    if (*fopt->eop) {
        free(navs.erp.data); navs.erp.data=NULL; navs.erp.n=navs.erp.nmax=0;
        reppath(fopt->eop,path,ts,"","");
        if (!readerp(path,&navs.erp)) {
            showmsg("error : no erp data %s",path);
            trace(2,"no erp data %s\n",path);
        }
    }
    /* read obs and nav data */
    if (!readobsnav(ts,te,ti,infile,index,n,&popt_,&obss,&navs,stas)) {
        /* free obs and nav data */
        freeobsnav(&obss, &navs);
        free(rtk_ptr);
        return 0;
    }

    /* read dcb parameters from DCB, BIA, BSX files */
    dcb_ok = 0;
    for (i=0;i<MAX_CODE_BIASES;i++) for (k=0;k<MAX_CODE_BIAS_FREQS;k++) {
        /* FIXME: cbias later initialized with 0 in readdcb()!  */
        for (j=0;j<MAXSAT;j++) navs.cbias[j][k][i]=-1;
        for (j=0;j<MAXRCV;j++) navs.rbias[j][k][i]=0;
        }
    for (i=0;i<n;i++) {  /* first check infiles for .BIA or .BSX files */
        if ((dcb_ok=readdcb(infile[i],&navs,stas))) break;
    }
    if (!dcb_ok&&*fopt->dcb) {  /* then check if DCB file specified */
        reppath(fopt->dcb,path,ts,"","");
        dcb_ok=readdcb(path,&navs,stas);
    }
    if (!dcb_ok) {

    }
    /* set antenna parameters */
    if (popt_.mode!=PMODE_SINGLE) {
        setpcv(obss.n>0?obss.data[0].time:timeget(),&popt_,&navs,&pcvss,&pcvsr,
               stas);
    }
    /* read ocean tide loading parameters */
    if (popt_.mode>PMODE_SINGLE&&*fopt->blq) {
        readotl(&popt_,fopt->blq,stas);
    }
    /* rover/reference fixed position */
    if (popt_.mode==PMODE_FIXED) {
        if (!antpos(&popt_,1,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
            free(rtk_ptr);
            return 0;
        }
        if (!antpos(&popt_,2,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
            free(rtk_ptr);
            return 0;
        }
    }
    else if (PMODE_DGPS<=popt_.mode&&popt_.mode<=PMODE_STATIC_START) {
        if (!antpos(&popt_,2,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
            free(rtk_ptr);
            return 0;
        }
    }
    /* open solution statistics */
    if (flag&&sopt->sstat>0) {
        strcpy(statfile,outfile);
        strcat(statfile,".stat");
        rtkclosestat();
        rtkopenstat(statfile,sopt->sstat);
    }
    /* write header to output file */
    if (flag&&!outhead(outfile,infile,n,&popt_,sopt)) {
        freeobsnav(&obss,&navs);
        free(rtk_ptr);
        return 0;
    }
    /* name time events file */
    namefiletm(outfiletm,outfile);
    /* write header to file with time marks */
    outhead(outfiletm,infile,n,&popt_,sopt);

    iobsu=iobsr=isbs=reverse=aborts=0;

    if (popt_.mode==PMODE_SINGLE||popt_.soltype==SOLTYPE_FORWARD) {
        FILE *fp=openfile(outfile);
        if (fp) {
            FILE *fptm=openfile(outfiletm);
            if (fptm) {
                rtkinit(rtk_ptr,&popt_);
                procpos(fp,fptm,&popt_,sopt,rtk_ptr,SOLMODE_SINGLE_DIR);
                 /* Pass-by-Pass AR: process 48h ambiguity resolution if enabled */
                if (popt_.armode_pbp > 0) {
                    int n_fixed = ppp_ar_48h(&popt_, rtk_ptr, &obss);
                    if (n_fixed > 0) {
                        trace(1, "Pass-by-Pass AR: %d ambiguities fixed\n", n_fixed);
                    }
                }
                rtkfree(rtk_ptr);
                fclose(fptm);
            }
            fclose(fp);
        }
    }
    else if (popt_.soltype==SOLTYPE_BACKWARD) {
        FILE *fp=openfile(outfile);
        if (fp) {
            FILE *fptm=openfile(outfiletm);
            if (fptm) {
                reverse=1; iobsu=iobsr=obss.n-1; isbs=sbss.n-1;
                rtkinit(rtk_ptr,&popt_);
                procpos(fp,fptm,&popt_,sopt,rtk_ptr,SOLMODE_SINGLE_DIR);
                rtkfree(rtk_ptr);
                fclose(fptm);
            }
            fclose(fp);
        }
    }
    else { /* combined or combined with no phase reset */
        solf=(sol_t *)malloc(sizeof(sol_t)*nepoch);
        solb=(sol_t *)malloc(sizeof(sol_t)*nepoch);
        rbf=(double *)malloc(sizeof(double)*nepoch*3);
        rbb=(double *)malloc(sizeof(double)*nepoch*3);

        if (solf&&solb) {
            isolf=isolb=0;
            rtkinit(rtk_ptr,&popt_);
            procpos(NULL,NULL,&popt_,sopt,rtk_ptr,SOLMODE_COMBINED); /* forward */
            reverse=1; iobsu=iobsr=obss.n-1; isbs=sbss.n-1;
            if (popt_.soltype!=SOLTYPE_COMBINED_NORESET) {
                /* Reset */
                rtkfree(rtk_ptr);
                rtkinit(rtk_ptr,&popt_);
            }
            procpos(NULL,NULL,&popt_,sopt,rtk_ptr,SOLMODE_COMBINED); /* backward */
             /* Pass-by-Pass AR: process 48h ambiguity resolution if enabled */
            if (popt_.armode_pbp > 0) {
                int n_fixed = ppp_ar_48h(&popt_, rtk_ptr, &obss);
                if (n_fixed > 0) {
                    trace(1, "Pass-by-Pass AR: %d ambiguities fixed\n", n_fixed);
                }
            }
            rtkfree(rtk_ptr);

            /* combine forward/backward solutions */
            if (!aborts) {
                FILE *fp=openfile(outfile);
                if (fp) {
                    FILE *fptm=openfile(outfiletm);
                    if (fptm) {
                        combres(fp,fptm,&popt_,sopt);
                        fclose(fptm);
                    }
                    fclose(fp);
                }
            }
        }
        else showmsg("error : memory allocation");
        free(solf);
        free(solb);
        free(rbf);
        free(rbb);
    }
    /* free rtk, obs and nav data */
    free(rtk_ptr);
    freeobsnav(&obss,&navs);

    return aborts?1:0;
}
/* execute processing session for each rover ---------------------------------*/
static int execses_r(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                     const solopt_t *sopt, const filopt_t *fopt, int flag,
                     const char **infile, const int *index, int n, const char *outfile,
                     const char *rov)
{
    gtime_t t0={0};
    int i,stat=0;
    char *ifile[MAXINFILE],ofile[1024],*rov_,*q,s[40]="";
    const char *p;

    trace(3,"execses_r: n=%d outfile=%s\n",n,outfile);

    for (i=0;i<n;i++) if (strstr(infile[i],"%r")) break;

    if (i<n) { /* include rover keywords */
        if (!(rov_=(char *)malloc(strlen(rov)+1))) return 0;
        strcpy(rov_,rov);

        for (i=0;i<n;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                free(rov_); for (;i>=0;i--) free(ifile[i]);
                return 0;
            }
        }
        for (p=rov_;;p=q+1) { /* for each rover */
            if ((q=strchr(p,' '))) *q='\0';

            if (*p) {
                strcpy(proc_rov,p);
                if (ts.time) time2str(ts,s,0); else *s='\0';
                if (checkbrk("reading    : %s",s)) {
                    stat=1;
                    break;
                }
                for (i=0;i<n;i++) reppath(infile[i],ifile[i],t0,p,"");
                reppath(outfile,ofile,t0,p,"");

                /* execute processing session */
                stat=execses(ts,te,ti,popt,sopt,fopt,flag,(const char **)ifile,index,n,ofile);
            }
            if (stat==1||!q) break;
        }
        free(rov_); for (i=0;i<n;i++) free(ifile[i]);
    }
    else {
        /* execute processing session */
        stat=execses(ts,te,ti,popt,sopt,fopt,flag,infile,index,n,outfile);
    }
    return stat;
}
/* execute processing session for each base station --------------------------*/
static int execses_b(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                     const solopt_t *sopt, const filopt_t *fopt, int flag,
                     const char **infile, const int *index, int n, const char *outfile,
                     const char *rov, const char *base)
{
    gtime_t t0={0};
    int i,stat=0;
    char *ifile[MAXINFILE],ofile[1024],*base_,*q,s[40];
    const char *p;

    trace(3,"execses_b: n=%d outfile=%s\n",n,outfile);

    /* read prec ephemeris and sbas data */
    readpreceph(infile,n,popt,&navs,&sbss);

    for (i=0;i<n;i++) if (strstr(infile[i],"%b")) break;

    if (i<n) { /* include base station keywords */
        if (!(base_=(char *)malloc(strlen(base)+1))) {
            freepreceph(&navs,&sbss);
            return 0;
        }
        strcpy(base_,base);

        for (i=0;i<n;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                free(base_); for (;i>=0;i--) free(ifile[i]);
                freepreceph(&navs,&sbss);
                return 0;
            }
        }
        for (p=base_;;p=q+1) { /* for each base station */
            if ((q=strchr(p,' '))) *q='\0';

            if (*p) {
                strcpy(proc_base,p);
                if (ts.time) time2str(ts,s,0); else *s='\0';
                if (checkbrk("reading    : %s",s)) {
                    stat=1;
                    break;
                }
                for (i=0;i<n;i++) reppath(infile[i],ifile[i],t0,"",p);
                reppath(outfile,ofile,t0,"",p);

                stat=execses_r(ts,te,ti,popt,sopt,fopt,flag,(const char **)ifile,index,n,(const char *)ofile,rov);
            }
            if (stat==1||!q) break;
        }
        free(base_); for (i=0;i<n;i++) free(ifile[i]);
    }
    else {
        stat=execses_r(ts,te,ti,popt,sopt,fopt,flag,infile,index,n,outfile,rov);
    }
    /* free prec ephemeris and sbas data */
    freepreceph(&navs,&sbss);

    return stat;
}
/* post-processing positioning -------------------------------------------------
* post-processing positioning
* args   : gtime_t ts       I   processing start time (ts.time==0: no limit)
*        : gtime_t te       I   processing end time   (te.time==0: no limit)
*          double ti        I   processing interval  (s) (0:all)
*          double tu        I   processing unit time (s) (0:all)
*          prcopt_t *popt   I   processing options
*          solopt_t *sopt   I   solution options
*          filopt_t *fopt   I   file options
*          char   **infile  I   input files (see below)
*          int    n         I   number of input files
*          char   *outfile  I   output file ("":stdout, see below)
*          char   *rov      I   rover id list        (separated by " ")
*          char   *base     I   base station id list (separated by " ")
* return : status (0:ok,0>:error,1:aborted)
* notes  : input files should contain observation data, navigation data, precise
*          ephemeris/clock (optional), sbas log file (optional), ssr message
*          log file (optional) and tec grid file (optional). only the first
*          observation data file in the input files is recognized as the rover
*          data.
*
*          the type of an input file is recognized by the file extension as ]
*          follows:
*              .sp3,.SP3,.eph*,.EPH*: precise ephemeris (sp3c)
*              .sbs,.SBS,.ems,.EMS  : sbas message log files (rtklib or ems)
*              .rtcm3,.RTCM3        : ssr message log files (rtcm3)
*              .*i,.*I              : tec grid files (ionex)
*              others               : rinex obs, nav, gnav, hnav, qnav or clock
*
*          inputs files can include wild-cards (*). if an file includes
*          wild-cards, the wild-card expanded multiple files are used.
*
*          inputs files can include keywords. if an file includes keywords,
*          the keywords are replaced by date, time, rover id and base station
*          id and multiple session analyses run. refer reppath() for the
*          keywords.
*
*          the output file can also include keywords. if the output file does
*          not include keywords. the results of all multiple session analyses
*          are output to a single output file.
*
*          ssr corrections are valid only for forward estimation.
*-----------------------------------------------------------------------------*/
extern int postpos(gtime_t ts, gtime_t te, double ti, double tu,
                   const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, const char **infile, int n,
                   const char *outfile, const char *rov, const char *base)
{
    gtime_t tts,tte,ttte;
    double tunit,tss;
    int i,j,k,nf,stat=0,week,flag=1,index[MAXINFILE]={0};
    char *ifile[MAXINFILE],ofile[1024];
    const char *ext;

    trace(3,"postpos : ti=%.0f tu=%.0f n=%d outfile=%s\n",ti,tu,n,outfile);

    /* open processing session */
    if (!openses(popt,sopt,fopt,&navs,&pcvss,&pcvsr)) return -1;

    /* --------------------------------------------------------------------- */
    /* case 1) time range specified (ts & te valid) */
    /* --------------------------------------------------------------------- */
    if (ts.time!=0 && te.time!=0 && tu>=0.0) {

        if (timediff(te,ts)<0.0) {
            showmsg("error : no period");
            closeses(&navs,&pcvss,&pcvsr);
            return 0;
        }

        for (i=0;i<MAXINFILE;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                for (;i>=0;i--) free(ifile[i]);
                closeses(&navs,&pcvss,&pcvsr);
                return -1;
            }
        }

        if (tu==0.0||tu>86400.0*MAXPRCDAYS) tu=86400.0*MAXPRCDAYS;

        /*
         * PATCH: If the output file name has no time keywords ('%'), do not split
         * processing into many small time-windows (tu). Otherwise each window may
         * re-initialize output/stat files in some variants and leave only the last
         * window in the final .pos/.stat files.
         */
        if (tu>0.0 && outfile && *outfile && !strchr(outfile,'%')) {
            tu=86400.0*MAXPRCDAYS;
        }
        settspan(ts,te);
        tunit=tu<86400.0?tu:86400.0;
        tss=tunit*(int)floor(time2gpst(ts,&week)/tunit);

        for (i=0;;i++) { /* for each periods */
            tts=gpst2time(week,tss+i*tu);
            tte=timeadd(tts,tu-DTTOL);
            if (timediff(tts,te)>0.0) break;
            if (timediff(tts,ts)<0.0) tts=ts;
            if (timediff(tte,te)>0.0) tte=te;

            strcpy(proc_rov ,"");
            strcpy(proc_base,"");
            if (checkbrk("reading    : %s",time_str(tts,0))) {
                stat=1;
                break;
            }

            for (j=k=nf=0;j<n;j++) {
                ext=strrchr(infile[j],'.');

                if (ext&&(!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) {
                    strcpy(ifile[nf++],infile[j]);
                }
                else {
                    ttte=tte;
                    if (ext&&(!strcmp(ext,".sp3")||!strcmp(ext,".SP3")||
                              !strcmp(ext,".eph")||!strcmp(ext,".EPH"))) {
                        ttte=timeadd(ttte,3600.0);
                    }
                    else if (strstr(infile[j],"brdc")) {
                        ttte=timeadd(ttte,7200.0);
                    }
                    nf+=reppaths(infile[j],ifile+nf,MAXINFILE-nf,tts,ttte,"","");
                }
                while (k<nf) index[k++]=j;

                if (nf>=MAXINFILE) {
                    trace(2,"too many input files. trancated\n");
                    break;
                }
            }

            if (!reppath(outfile,ofile,tts,"","") && i>0) flag=0;

            stat=execses_b(tts,tte,ti,popt,sopt,fopt,flag,
                           (const char **)ifile,index,nf,(const char *)ofile,rov,base);

            if (stat==0) {
                processed_days++;
                printf("Day %d processing completed\n", processed_days);
            }
            if (stat==1) break;
        }

        for (i=0;i<MAXINFILE;i++) free(ifile[i]);
    }

    /* --------------------------------------------------------------------- */
    /* case 2) only ts specified */
    /* --------------------------------------------------------------------- */
    else if (ts.time!=0) {

        for (i=0;i<n&&i<MAXINFILE;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                for (;i>=0;i--) free(ifile[i]);
                closeses(&navs,&pcvss,&pcvsr);
                return -1;
            }
            reppath(infile[i],ifile[i],ts,"","");
            index[i]=i;
        }
        reppath(outfile,ofile,ts,"","");

        stat=execses_b(ts,te,ti,popt,sopt,fopt,1,
                       (const char **)ifile,index,n,ofile,rov,base);

        if (stat==0) {
            processed_days++;
            printf("Day %d processing completed\n", processed_days);
        }

        for (i=0;i<n&&i<MAXINFILE;i++) free(ifile[i]);
    }

    /* --------------------------------------------------------------------- */
    /* case 3) ts/te not specified: auto-detect obs range and process by natural day */
    /* --------------------------------------------------------------------- */
    else {

        obs_t obs_temp={0};
        nav_t nav_temp={0};
        sta_t sta_temp[MAXRCV]={{0}};
        gtime_t ts_auto={0}, te_auto={0};

        for (i=0;i<n;i++) index[i]=i;

        for (i=0;i<MAXINFILE;i++) {
            if (!(ifile[i]=(char *)malloc(1024))) {
                for (;i>=0;i--) free(ifile[i]);
                closeses(&navs,&pcvss,&pcvsr);
                return -1;
            }
        }

        trace(3,"Reading obs files to determine time range\n");

        if (!readobsnav(ts,te,ti,infile,index,n,popt,&obs_temp,&nav_temp,sta_temp)) {
            if (obs_temp.data) free(obs_temp.data);
            for (i=0;i<MAXINFILE;i++) free(ifile[i]);
            closeses(&navs,&pcvss,&pcvsr);
            return -1;
        }

        if (obs_temp.n > 0) {

            /* robust scan min/max time */
            ts_auto = obs_temp.data[0].time;
            te_auto = obs_temp.data[0].time;
            for (i=1;i<obs_temp.n;i++) {
                if (timediff(obs_temp.data[i].time, ts_auto) < 0.0) ts_auto = obs_temp.data[i].time;
                if (timediff(obs_temp.data[i].time, te_auto) > 0.0) te_auto = obs_temp.data[i].time;
            }
            free(obs_temp.data);
            obs_temp.data = NULL;

            /* align to natural day boundary */
            {
                int week_start, week_end;
                double sow_start = time2gpst(ts_auto, &week_start);
                double sow_end   = time2gpst(te_auto, &week_end);

                int doy_start = (int)(sow_start/86400.0);
                int doy_end   = (int)(sow_end  /86400.0);

                gtime_t day_aligned_start = gpst2time(week_start, doy_start*86400.0);
                gtime_t day_aligned_end   = gpst2time(week_end, (doy_end+1)*86400.0 - DTTOL);

                int num_days = doy_end - doy_start + 1;
                if (week_end != week_start) {
                    num_days = (int)((timediff(day_aligned_end, day_aligned_start)+DTTOL)/86400.0) + 1;
                }

                printf("Auto-detected obs range: %s to %s\n",
                       time_str(ts_auto,0), time_str(te_auto,0));
                printf("Day-aligned processing: %s to %s (%d natural days)\n",
                       time_str(day_aligned_start,0), time_str(day_aligned_end,0), num_days);

                /* ============================================================= */
                /* Pass-by-Pass driver (only when >=2 days) */
                /* ============================================================= */
                if (popt->armode_pbp >= 1 && num_days >= 2) {

                    /* these are provided by your PBP modules */
                    extern void init_arc_data(void);
                    extern int  ppp_ar_48h(const prcopt_t *popt, rtk_t *rtk, const obs_t *obs);
                    extern int  pbp_day_tag, pbp_collect_flag, pbp_apply_flag, pbp_resolve_flag, pbp_neq_accum_flag;
                    extern int  pbp_neq_init(gtime_t t0, gtime_t t1, double ti, const prcopt_t *opt);
                    extern int  pbp_finalize_final_neq(void);
                    extern int  n_ddamb;
                    extern int  refsat;

                    int day_ok0 = 0, day_ok1 = 0;
                    prcopt_t popt_run = *popt;

                    popt_run.armode_pbp = 0; /* prevent nested ppp_ar_48h() call inside execses */

                    /* IMPORTANT: avoid ppp_ar_48h calling apply_ar_fixed on NULL rtk */
                    prcopt_t popt_ar  = *popt;
                    if (popt_ar.armode_pbp >= 3) popt_ar.armode_pbp = 2;

                    init_arc_data();
                    n_ddamb = 0;
                    refsat  = 0;

                    /* ---------- PASS1: day0/day1 float + collect ---------- */
                    pbp_collect_flag = 1;
                    pbp_apply_flag   = 0;

                    for (i=0;i<2;i++) {

                        gtime_t day_start, day_end;
                        int current_week;
                        double current_sow;

                        current_sow  = (doy_start + i) * 86400.0;
                        current_week = week_start;
                        while (current_sow >= 604800.0) {
                            current_sow -= 604800.0;
                            current_week++;
                        }
                        day_start = gpst2time(current_week, current_sow);
                        day_end   = timeadd(day_start, 86400.0 - DTTOL);

                        if (timediff(day_start, ts_auto) < 0.0) day_start = ts_auto;
                        if (timediff(day_end,   te_auto) > 0.0) day_end   = te_auto;

                        printf("\n[PBP] Pass1 (float+collect) day %d: %s to %s\n", i+1,
                               time_str(day_start,0), time_str(day_end,0));

                        /* build file list */
                        for (j=k=nf=0;j<n;j++) {

                            ttte = day_end;

                            /* IMPORTANT: B2b may need next day (you can adjust here if you want) */
                            if (strstr(infile[j],"b2b") || strstr(infile[j],"B2b")) {
                                ttte = timeadd(ttte, 86400.0);
                            }
                            else if (strstr(infile[j],"brdc")) {
                                ttte = timeadd(ttte, 7200.0);
                            }

                            nf += reppaths(infile[j], ifile+nf, MAXINFILE-nf, day_start, ttte, "", "");
                            while (k<nf) index[k++]=j;

                            if (nf>=MAXINFILE) {
                                trace(2,"too many input files. truncated\n");
                                break;
                            }
                        }

                        /* output float file (avoid overwrite) */
                        if (outfile && *outfile) {
                            char tmp[1024];
                            reppath(outfile, tmp, day_start, "", "");
                            snprintf(ofile, sizeof(ofile), "%s.pbpfloatD%d", tmp, i+1);
                        }
                        else {
                            strcpy(ofile, "");
                        }

                        pbp_day_tag = i; /* force day tag for ambiguity collector */

                        /* ensure B2b SSR restarts from the correct daily file */
                        reset_B2b_reader(i);

                        stat = execses_b(day_start, day_end, ti, &popt_run, sopt, fopt, 1,
                                         (const char **)ifile, index, nf, ofile, rov, base);

                        if (i==0 && stat==0) day_ok0 = 1;
                        if (i==1 && stat==0) day_ok1 = 1;
                        if (stat==1) break; /* user abort */
                    }

                    pbp_collect_flag = 0;
                    pbp_day_tag = -1;

                    if (stat!=1 && day_ok0 && day_ok1) {

                        printf("\n[PBP] Solving WL/NL and building fixed DD set ...\n");

                        if (ppp_ar_48h(&popt_ar, NULL, NULL) > 0 && popt->armode_pbp >= 3) {

                            gtime_t sess_start = day_aligned_start;
                            gtime_t sess_end   = timeadd(day_aligned_start, 2.0 * 86400.0 - DTTOL);
                            if (timediff(sess_end, te_auto) > 0.0) sess_end = te_auto;

                            printf("\n[PBP] Pass2 (48h float accumulation for final NEQ): %s to %s\n",
                                   time_str(sess_start,0), time_str(sess_end,0));

                            for (j=k=nf=0;j<n;j++) {
                                ttte = sess_end;
                                if (strstr(infile[j],"b2b") || strstr(infile[j],"B2b")) ttte = timeadd(ttte, 86400.0);
                                else if (strstr(infile[j],"brdc")) ttte = timeadd(ttte, 7200.0);
                                nf += reppaths(infile[j], ifile+nf, MAXINFILE-nf, sess_start, ttte, "", "");
                                while (k<nf) index[k++]=j;
                                if (nf>=MAXINFILE) { trace(2,"too many input files. truncated\n"); break; }
                            }

                            reppath(outfile, ofile, sess_start, "", "");
                            remove(ofile);
                            rtkclosestat();

                            pbp_day_tag      = -1;
                            pbp_collect_flag = 0;
                            pbp_apply_flag   = 0;
                            pbp_resolve_flag = 0;
                            pbp_neq_accum_flag = 0;
                            (void)pbp_neq_init(sess_start, sess_end, ti, &popt_run);
                            pbp_neq_accum_flag = 1;

                            reset_B2b_reader(0);
                            stat = execses_b(sess_start, sess_end, ti, &popt_run, sopt, fopt, 1,
                                             (const char **)ifile, index, nf, ofile, rov, base);
                            pbp_neq_accum_flag = 0;

                            if (stat==0 && pbp_finalize_final_neq()) {
                                printf("\n[PBP] Pass3 (write fixed clock series from solved NEQ)\n");
                                remove(ofile);
                                rtkclosestat();
                                pbp_resolve_flag = 1;
                                reset_B2b_reader(0);
                                stat = execses_b(sess_start, sess_end, ti, &popt_run, sopt, fopt, 1,
                                                 (const char **)ifile, index, nf, ofile, rov, base);
                                pbp_resolve_flag = 0;
                            }

                            pbp_day_tag = -1;
                            if (stat==0) processed_days = 2;
                        }
                        else {
                            /* only collect/solve, or AR failed */
                            processed_days = 2;
                        }
                    }
                    else {
                        /* fallback: original per-day processing */
                        printf("[PBP] Warning: pass1 failed, fallback to per-day processing.\n");
                        for (i=0;i<num_days;i++) {
                            /* (reuse the original per-day loop below if you want) */
                        }
                    }
                }

                /* ============================================================= */
                /* Original per-day processing (when not PBP) */
                /* ============================================================= */
                else {

                    for (i=0; i<num_days; i++) {

                        gtime_t day_start, day_end;
                        int current_week;
                        double current_sow;

                        current_sow  = (doy_start + i) * 86400.0;
                        current_week = week_start;

                        while (current_sow >= 604800.0) {
                            current_sow -= 604800.0;
                            current_week++;
                        }

                        day_start = gpst2time(current_week, current_sow);
                        day_end   = timeadd(day_start, 86400.0 - DTTOL);

                        if (timediff(day_start, ts_auto) < 0.0) day_start = ts_auto;
                        if (timediff(day_end,   te_auto) > 0.0) day_end   = te_auto;

                        printf("\nProcessing natural day %d: %s to %s\n", i+1,
                               time_str(day_start,0), time_str(day_end,0));

                        for (j=k=nf=0;j<n;j++) {

                            ext=strrchr(infile[j],'.');

                            if (ext&&(!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) {
                                strcpy(ifile[nf++],infile[j]);
                            }
                            else {
                                ttte=day_end;
                                if (ext&&(!strcmp(ext,".sp3")||!strcmp(ext,".SP3")||
                                          !strcmp(ext,".eph")||!strcmp(ext,".EPH"))) {
                                    ttte=timeadd(ttte,3600.0);
                                }
                                else if (strstr(infile[j],"brdc")) {
                                    ttte=timeadd(ttte,7200.0);
                                }
                                nf+=reppaths(infile[j],ifile+nf,MAXINFILE-nf,day_start,ttte,"","");
                            }
                            while (k<nf) index[k++]=j;

                            if (nf>=MAXINFILE) {
                                trace(2,"too many input files. truncated\n");
                                break;
                            }
                        }

                        reppath(outfile,ofile,day_start,"","");

                        /* ensure B2b SSR starts from the correct daily file */
                        reset_B2b_reader(i);

                        stat=execses_b(day_start,day_end,ti,popt,sopt,fopt,1,
                                       (const char **)ifile,index,nf,ofile,rov,base);

                        if (stat==0) {
                            processed_days++;
                            printf("Day %d processing completed\n", processed_days);
                        }
                        else {
                            printf("Day %d processing failed\n", i+1);
                            if (stat==1) break;
                        }
                    }
                }
            } /* end day-align scope */
        }

        /* obs_temp.n == 0 */
        else {
            if (obs_temp.data) free(obs_temp.data);
            stat=execses_b(ts,te,ti,popt,sopt,fopt,1,infile,index,n,outfile,rov,base);
            if (stat==0) {
                processed_days++;
                printf("Day %d processing completed\n", processed_days);
            }
        }

        for (i=0;i<MAXINFILE;i++) free(ifile[i]);
    }

    /* close processing session */
    closeses(&navs,&pcvss,&pcvsr);

    /* Output processing summary */
    if (processed_days == 0) {
        printf("\nNo complete day processed\n");
    } else if (processed_days == 1) {
        printf("\n1 day only\n");
    } else if (processed_days == 2) {
        printf("\n2 days complete\n");
    } else {
        printf("\n%d days complete\n", processed_days);
    }

    processed_days = 0;
    return stat;
}
