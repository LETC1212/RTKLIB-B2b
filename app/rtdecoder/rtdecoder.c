#define _POSIX_C_SOURCE 199506
#include <stdlib.h>
#include <signal.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <netdb.h>
#endif
#include <errno.h>
#include "rtklib.h"
#include "rtdecoder.h"
#include "B2b.h"
#include "vt.h"


#define PRGNAME     "B2bPPP"            /* program name */
#define CMDPROMPT   "B2bPPP> "          /* command prompt */
#define MAXCON      32                  /* max number of consoles */
#define MAXARG      10                  /* max number of args in a command */
#define MAXCMD      256                 /* max length of a command */
#define MAXSTR      1024                /* max length of a stream */
#define OPTSDIR     "../conf/"          /* default config directory */
#define OPTSFILE    "B2bPPP.conf"       /* default config file */
#define NAVIFILE    "rtkrcv.nav"        /* navigation save file */
#define STATFILE    "../log/B2bPPP_%Y_%m_%d.stat"  /* solution status file */
#define TRACEFILE   "./B2bPPP_%Y_%m_%d.trace" /* debug trace file */
#define B2bTRACEFILE   "./B2bPPP_%Y_%m_%d.B2bssr" /* debug trace file */
#define LOGFILE     "../log/B2bPPP_%Y_%m_%d.log"   /* Deamon log file */
#define INTKEEPALIVE 1000               /* keep alive interval (ms) */
#define MIN_INT_RESET   30000   /* mininum interval of reset command (ms) */

#define ESC_CLEAR   "\033[H\033[2J"     /* ansi/vt100 escape: erase screen */
#define ESC_RESET   "\033[0m"           /* ansi/vt100: reset attribute */
#define ESC_BOLD    "\033[1m"           /* ansi/vt100: bold */

#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))


#ifdef _WIN32
typedef rtklib_thread_t pthread_t;
typedef int socklen_t;
#endif

/* type definitions ----------------------------------------------------------*/

/* global variables ----------------------------------------------------------*/
static rtksvr_t svr;                    /* rtk server struct */
static stream_t moni;                   /* monitor stream */

static int intflg       =0;             /* interrupt flag (2:shutdown) */

static char passwd[MAXSTR]="admin";     /* login password */
static int timetype     =0;             /* time format (0:gpst,1:utc,2:jst,3:tow) */
static int soltype      =0;             /* sol format (0:dms,1:deg,2:xyz,3:enu,4:pyl) */
static int solflag      =2;             /* sol flag (1:std+2:age/ratio/ns) */
static int strtype[]={                  /* stream types */
    STR_SERIAL,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE,STR_NONE
};
static char strpath[8][MAXSTR]={"","","","","","","",""}; /* stream paths */
static int strfmt[]={                   /* stream formats */
    STRFMT_UBX,STRFMT_RTCM3,STRFMT_SP3,SOLF_LLH,SOLF_NMEA
};
static char rcvopt[3][256]={""};        /* Receiver options */
static int svrcycle     =10;            /* server cycle (ms) */
static int timeout      =10000;         /* timeout time (ms) */
static int reconnect    =10000;         /* reconnect interval (ms) */
static int nmeacycle    =5000;          /* nmea request cycle (ms) */
static int buffsize     =32768;         /* input buffer size (bytes) */
static int navmsgsel    =0;             /* navigation message select */
static char proxyaddr[256]="";          /* http/ntrip proxy */
static int nmeareq      =0;             /* nmea request type (0:off,1:lat/lon,2:single) */
static double nmeapos[] ={0,0,0};       /* nmea position (lat/lon/height) (deg,m) */
static char rcvcmds[3][MAXSTR]={""};    /* receiver commands files */
#ifdef RTKSHELLCMDS
static char startcmd[MAXSTR]="";        /* start command */
static char stopcmd [MAXSTR]="";        /* stop command */
#endif
static int modflgr[256] ={0};           /* modified flags of receiver options */
static int modflgs[256] ={0};           /* modified flags of system options */
static int moniport     =0;             /* monitor port */
static int keepalive    =0;             /* keep alive flag */
static int start        =0;             /* auto start */
static int fswapmargin  =30;            /* file swap margin (s) */
static char sta_name[256]="";           /* station name */

static prcopt_t prcopt;                 /* processing options */
static solopt_t solopt[2]={{0}};        /* solution options */
static filopt_t filopt  ={""};          /* file options */

PPPGlobal_t PPP_Glo = {0};

/* type definitions ----------------------------------------------------------*/

typedef struct {                       /* console type */
    int state;                         /* state (0:stop,1:run) */
    vt_t *vt;                          /* virtual terminal */
    pthread_t thread;                  /* console thread */
} con_t;


/* help text -----------------------------------------------------------------*/
static const char *usage[]={
    "usage: B2bPPP [-s][-p port][-d dev][-o file][-w pwd][-r level][-t level][-sta sta]",
    "options",
    "  -s         start RTK server on program startup",
    "  -nc        start RTK server on program startup with no console",
    "  -p port    port number for telnet console",
    "  -m port    port number for monitor stream",
    "  -d dev     terminal device for console",
    "  -o file    processing options file",
    "  -w pwd     login password for remote console (\"\": no password)",
    "  -r level   output solution status file (0:off,1:states,2:residuals)",
    "  -t level   debug trace level (0:off,1-5:on)",
    "  -sta sta   station name for receiver dcb",
    "  --deamon   detach from the console",
    "  --version  print the version and exit"
};
static const char *helptxt[]={
    "start                 : start rtk server",
    "stop                  : stop rtk server",
    "restart               : restart rtk sever",
    "solution [cycle]      : show solution",
    "status [cycle]        : show rtk status",
    "satellite [-n] [cycle]: show satellite status",
    "observ [-n] [cycle]   : show observation data",
    "navidata [cycle]      : show navigation data",
    "stream [cycle]        : show stream status",
    "ssr [cycle]           : show ssr corrections",
    "error                 : show error/warning messages",
    "option [opt]          : show option(s)",
    "set opt [val]         : set option",
    "load [file]           : load options from file",
    "save [file]           : save options to file",
    "log [file|off]        : start/stop log to file",
    "help|? [path]         : print help",
    "exit|ctr-D            : logout console (only for telnet)",
    "shutdown              : shutdown rtk server",
    ""
};
static const char *pathopts[]={         /* path options help */
    "stream path formats",
    "serial   : port[:bit_rate[:byte[:parity(n|o|e)[:stopb[:fctr(off|on)[#port]]]]]]]",
    "file     : path[::T[::+offset][::xspeed]]",
    "tcpsvr   : :port",
    "tcpcli   : addr:port",
    "ntripsvr : [passwd@]addr:port/mntpnt[:str]",
    "ntripcli : user:passwd@addr:port/mntpnt",
    "ntripcas : user:passwd@:[port]/mpoint[:srctbl]",
    "ftp      : user:passwd@addr/path[::T=poff,tint,off,rint]",
    "http     : addr/path[::T=poff,tint,off,rint]",
    ""
};
/* receiver options table ----------------------------------------------------*/
#define TIMOPT  "0:gpst,1:utc,2:jst,3:tow"
#define CONOPT  "0:dms,1:deg,2:xyz,3:enu,4:pyl"
#define FLGOPT  "0:off,1:std+2:age/ratio/ns"
#define ISTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,6:ntripcli,7:ftp,8:http"
#define OSTOPT  "0:off,1:serial,2:file,3:tcpsvr,4:tcpcli,5:ntripsvr,9:ntripcas,11:udpcli"
#define FMTOPT  "0:rtcm2,1:rtcm3,2:oem4,4:ubx,5:swift,6:hemis,7:skytraq,8:javad,9:nvs,10:binex,11:rt17,12:sbf,14,15:sp3"
#define NMEOPT  "0:off,1:latlon,2:single"
#define SOLOPT  "0:llh,1:xyz,2:enu,3:nmea,4:stat"
#define MSGOPT  "0:all,1:rover,2:base,3:corr"

static opt_t rcvopts[]={
    {"console-passwd",  2,  (void *)passwd,              ""     },
    {"console-timetype",3,  (void *)&timetype,           TIMOPT },
    {"console-soltype", 3,  (void *)&soltype,            CONOPT },
    {"console-solflag", 0,  (void *)&solflag,            FLGOPT },
    
    {"inpstr1-type",    3,  (void *)&strtype[0],         ISTOPT },
    {"inpstr2-type",    3,  (void *)&strtype[1],         ISTOPT },
    {"inpstr3-type",    3,  (void *)&strtype[2],         ISTOPT },
    {"inpstr1-path",    2,  (void *)strpath [0],         ""     },
    {"inpstr2-path",    2,  (void *)strpath [1],         ""     },
    {"inpstr3-path",    2,  (void *)strpath [2],         ""     },
    {"inpstr1-format",  3,  (void *)&strfmt [0],         FMTOPT },
    {"inpstr2-format",  3,  (void *)&strfmt [1],         FMTOPT },
    {"inpstr3-format",  3,  (void *)&strfmt [2],         FMTOPT },
    {"inpstr1-rcvopt",  2,  (void *)rcvopt[0],           ""     },
    {"inpstr2-rcvopt",  2,  (void *)rcvopt[1],           ""     },
    {"inpstr3-rcvopt",  2,  (void *)rcvopt[2],           ""     },
    {"inpstr2-nmeareq", 3,  (void *)&nmeareq,            NMEOPT },
    {"inpstr2-nmealat", 1,  (void *)&nmeapos[0],         "deg"  },
    {"inpstr2-nmealon", 1,  (void *)&nmeapos[1],         "deg"  },
    {"inpstr2-nmeahgt", 1,  (void *)&nmeapos[2],         "m"    },
    {"outstr1-type",    3,  (void *)&strtype[3],         OSTOPT },
    {"outstr2-type",    3,  (void *)&strtype[4],         OSTOPT },
    {"outstr1-path",    2,  (void *)strpath [3],         ""     },
    {"outstr2-path",    2,  (void *)strpath [4],         ""     },
    {"outstr1-format",  3,  (void *)&strfmt [3],         SOLOPT },
    {"outstr2-format",  3,  (void *)&strfmt [4],         SOLOPT },
    {"logstr1-type",    3,  (void *)&strtype[5],         OSTOPT },
    {"logstr2-type",    3,  (void *)&strtype[6],         OSTOPT },
    {"logstr3-type",    3,  (void *)&strtype[7],         OSTOPT },
    {"logstr1-path",    2,  (void *)strpath [5],         ""     },
    {"logstr2-path",    2,  (void *)strpath [6],         ""     },
    {"logstr3-path",    2,  (void *)strpath [7],         ""     },
    
    {"misc-svrcycle",   0,  (void *)&svrcycle,           "ms"   },
    {"misc-timeout",    0,  (void *)&timeout,            "ms"   },
    {"misc-reconnect",  0,  (void *)&reconnect,          "ms"   },
    {"misc-nmeacycle",  0,  (void *)&nmeacycle,          "ms"   },
    {"misc-buffsize",   0,  (void *)&buffsize,           "bytes"},
    {"misc-navmsgsel",  3,  (void *)&navmsgsel,          MSGOPT },
    {"misc-proxyaddr",  2,  (void *)proxyaddr,           ""     },
    {"misc-fswapmargin",0,  (void *)&fswapmargin,        "s"    },
    
#ifdef RTKSHELLCMDS
    {"misc-startcmd",   2,  (void *)startcmd,            ""     },
    {"misc-stopcmd",    2,  (void *)stopcmd,             ""     },
#endif
    
    {"file-cmdfile1",   2,  (void *)rcvcmds[0],          ""     },
    {"file-cmdfile2",   2,  (void *)rcvcmds[1],          ""     },
    {"file-cmdfile3",   2,  (void *)rcvcmds[2],          ""     },
    
    {"",0,NULL,""}
};


/* discard space characters at tail ------------------------------------------*/
static void chop(char *str)
{
    char *p;
    for (p=str+strlen(str)-1;p>=str&&!isgraph((int)*p);p--) *p='\0';
}

/* update ephemeris ----------------------------------------------------------*/
static void update_eph(rtksvr_t *svr, nav_t *nav, int ephsat, int ephset,
                       int index)
{
    eph_t *eph1,*eph2,*eph3;
    geph_t *geph1,*geph2,*geph3;
    int prn;
    
    if (satsys(ephsat,&prn)!=SYS_GLO) {
            if (!svr->navsel||svr->navsel==index+1) {
            /* svr->nav.eph={current_set1,current_set2,prev_set1,prev_set2} */
            eph1=nav->eph+ephsat-1+MAXSAT*ephset;         /* received */
            eph2=svr->nav.eph+ephsat-1+MAXSAT*ephset;     /* current */
            eph3=svr->nav.eph+ephsat-1+MAXSAT*(2+ephset); /* previous */

            trace(22, "update_eph: ephsat=%d ephset=%d index=%d\n", ephsat, ephset, index);
            trace(22, "eph1 (received): iode=%d, toe=%s, ttr=%s\n", eph1->iode, time_str(eph1->toe, 0), time_str(eph1->ttr, 0));
            trace(22, "eph2 (current): iode=%d, toe=%s, ttr=%s\n", eph2->iode, time_str(eph2->toe, 0), time_str(eph2->ttr, 0));
            trace(22, "eph3 (previous): iode=%d, toe=%s, ttr=%s\n", eph3->iode, time_str(eph3->toe, 0), time_str(eph3->ttr, 0));

                if (eph2->ttr.time==0||
                    (eph1->iode!=eph3->iode&&eph1->iode!=eph2->iode)||
                    (timediff(eph1->toe,eph3->toe)!=0.0&&
                 timediff(eph1->toe,eph2->toe)!=0.0)||
                (timediff(eph1->toc,eph3->toc)!=0.0&&
                 timediff(eph1->toc,eph2->toc)!=0.0)) {
                *eph3=*eph2; /* current ->previous */
                *eph2=*eph1; /* received->current */
                }
            }
            svr->nmsg[index][1]++;
        }
    }

static void update_B2bssr(rtksvr_t *svr, int index)
{
    int i;
    char satid[8], B2btime_str[64] = {0}, t0_str[6][64];
    raw_t *raw = &svr->raw[index];  // Direct pointer access for efficiency
    nav_t *nav = &svr->nav;         // Pointer to global navigation data

    time2str(raw->time, B2btime_str, 3);  // Initialize time string

    /* Output message header information */
    trace(22,"\n================ B2b Update @ %s ================\n", B2btime_str);
    if (raw->num_PPPB2BINF01 != 0)
        trace(22,"[INF01] GeoPRN=%d Count=%d\n", raw->geoprn, raw->num_PPPB2BINF01);
    if (raw->num_PPPB2BINF02 != 0)
        trace(22,"[INF02] Orbit Update Count=%d\n", raw->num_PPPB2BINF02);
    if (raw->num_PPPB2BINF03 != 0)
        trace(22,"[INF03] CodeBias Update Count=%d\n", raw->num_PPPB2BINF03);
    if (raw->num_PPPB2BINF04 != 0)
        trace(22,"[INF04] Clock Update Count=%d\n", raw->num_PPPB2BINF04);

    for (i = 0; i < MAXSAT; i++) {
        satno2id(i, satid);
        trace(22,"\n-- Processing SAT: %s --\n", satid);

        /* Skip satellites without updates */
        if (!raw->nav.B2bssr[i].update) {
            trace(22,"No update flag, skip\n");
            continue;
        }

        /* Time validity check */
        if (timediff(raw->time, svr->rtk.sol.time) < -1E-3) {
            trace(22,"Early data (B2b_time < sol_time), skip\n");
            raw->nav.B2bssr[i].update = 0;
            continue;
        }

        /* Message type dispatching */
        if (raw->num_PPPB2BINF02 != 0) {  // Orbit message
            /* Orbit data duplication check */
            if (fabs(timediff(raw->nav.B2bssr[i].t0[0], nav->B2bssr[i].t0[0])) < 1E-3) {
                if (!checkout_B2beph(&raw->nav.B2bssr[i], &nav->B2bssr[i])) {
                    printf("WARNING: Orbit data changed without t0 update!\n");
                    trace(22,"WARNING: Orbit data changed without t0 update!\n");
                }
                raw->nav.B2bssr[i].update = 0;
                trace(22,"Duplicate orbit data, skip\n");
                continue;
            }
            
            /* Calculate update interval */
            raw->nav.B2bssr[i].udi[0] = timediff(raw->nav.B2bssr[i].t0[0], 
                                               nav->B2bssr[i].t0[0]);
            if(raw->nav.B2bssr[i].udi[0] > 86400) raw->nav.B2bssr[i].udi[0] = 0;
            
            /* IODCOR consistency check */
            if (raw->nav.B2bssr[i].iodcorr[0] != raw->nav.B2bssr[i].iodcorr[1]) {
                trace(22,"ERROR: Orbit IODC mismatch (%d vs %d)\n",
                      raw->nav.B2bssr[i].iodcorr[0], 
                      raw->nav.B2bssr[i].iodcorr[1]);
                raw->nav.B2bssr[i].update = 0;
                continue;
            }
            /* B2b IODN vs Eph IODC consistency check */
            if(i>0){
                if (raw->nav.B2bssr[i].iodn != nav->eph[i-1].iodc && 
                    raw->nav.B2bssr[i].iodn != nav->eph[i-1+2*MAXSAT].iodc) {
                trace(22,"ERROR: B2b Orbit IODN mismatch EPH IODC (%d vs %d %d)\n",
                      raw->nav.B2bssr[i].iodn, 
                      nav->eph[i-1].iodc, nav->eph[i-1+2*MAXSAT].iodc);
                raw->nav.B2bssr[i].update = 0;
                continue;
                }
            }
            trace(22,"Orbit Updated: udi=%.1fs\n", raw->nav.B2bssr[i].udi[0]);
        }
        else if (raw->num_PPPB2BINF03 != 0) {  // Code bias message
            /* Code bias duplication check */
            if (fabs(timediff(raw->nav.B2bssr[i].t0[1], nav->B2bssr[i].t0[1])) < 1E-3) {
                if (!checkout_B2bcbia(&raw->nav.B2bssr[i], &nav->B2bssr[i])) {
                    printf("WARNING: CodeBias changed without t0 update!\n");
                    trace(22,"WARNING: CodeBias changed without t0 update!\n");
                }
                raw->nav.B2bssr[i].update = 0;
                trace(22,"Duplicate CodeBias, skip\n");
                continue;
            }
            
            raw->nav.B2bssr[i].udi[1] = timediff(raw->nav.B2bssr[i].t0[1],
                                                nav->B2bssr[i].t0[1]);
            // Set to 0 if received messages < 2
            if(raw->nav.B2bssr[i].udi[1] > 86400) raw->nav.B2bssr[i].udi[1] = 0;

            /* IODC consistency check */
            if (raw->nav.B2bssr[i].iodcorr[0] != raw->nav.B2bssr[i].iodcorr[1]) {
                trace(22,"ERROR: Orbit IODC mismatch (%d vs %d)\n",
                      raw->nav.B2bssr[i].iodcorr[0], 
                      raw->nav.B2bssr[i].iodcorr[1]);
                raw->nav.B2bssr[i].update = 0;
                continue;
            }
            /* B2b IODN vs Eph IODC consistency check */
            if(i>0){
                if (raw->nav.B2bssr[i].iodn != nav->eph[i-1].iodc && 
                    raw->nav.B2bssr[i].iodn != nav->eph[i-1+2*MAXSAT].iodc) {
                trace(22,"ERROR: B2b Orbit IODN mismatch EPH IODC (%d vs %d %d)\n",
                      raw->nav.B2bssr[i].iodn, 
                      nav->eph[i-1].iodc, nav->eph[i-1+2*MAXSAT].iodc);
                raw->nav.B2bssr[i].update = 0;
                continue;
                }
            }
            trace(22,"CodeBias Updated: udi=%.1fs\n", raw->nav.B2bssr[i].udi[1]);
        }
        else if (raw->num_PPPB2BINF04 != 0) {  // Clock message
            /* Clock duplication check */
            if (fabs(timediff(raw->nav.B2bssr[i].t0[2], nav->B2bssr[i].t0[2])) < 1E-3) {
                if (!checkout_B2bclk(&raw->nav.B2bssr[i], &nav->B2bssr[i])) {
                    trace(22,"WARNING: Clock changed without t0 update!\n");
                }
                raw->nav.B2bssr[i].update = 0;
                trace(22,"Duplicate Clock, skip\n");
                continue;
            }
            
            raw->nav.B2bssr[i].udi[2] = timediff(raw->nav.B2bssr[i].t0[2],
                                                nav->B2bssr[i].t0[2]);
            if(raw->nav.B2bssr[i].udi[2] > 86400) raw->nav.B2bssr[i].udi[2] = 0;
            /* IODC consistency check */
            if (raw->nav.B2bssr[i].iodcorr[0] != raw->nav.B2bssr[i].iodcorr[1]) {
                trace(22,"ERROR: Clock IODC mismatch (%d vs %d)\n",
                      raw->nav.B2bssr[i].iodcorr[0], 
                      raw->nav.B2bssr[i].iodcorr[1]);
                raw->nav.B2bssr[i].update = 0;
                continue;
            }
            /* B2b IODN vs Eph IODC consistency check */
            if(i>0){
                if (raw->nav.B2bssr[i].iodn != nav->eph[i-1].iodc && 
                    raw->nav.B2bssr[i].iodn != nav->eph[i-1+2*MAXSAT].iodc) {
                trace(22,"ERROR: B2b Clock IODN mismatch EPH IODC (%d vs %d %d)\n",
                      raw->nav.B2bssr[i].iodn, 
                      nav->eph[i-1].iodc, nav->eph[i-1+2*MAXSAT].iodc);
                raw->nav.B2bssr[i].update = 0;
                continue;
                }
            }
            trace(22,"Clock Updated: udi=%.1fs\n", raw->nav.B2bssr[i].udi[2]);
        }

        /* Execute data update */
        nav->B2bssr[i] = raw->nav.B2bssr[i];
        raw->nav.B2bssr[i].update = 0;  // Clear update flag
        svr->nmsg[index][8]++; // B2bssr counter update

        /* Detailed data logging */
        for (int j=0;j<3;j++) 
            time2str(nav->B2bssr[i].t0[j], t0_str[j], 3);
            
        trace(22,"Final Parameters:\n");
        trace(22,"  t0_orbit: %s\n", t0_str[0]);
        trace(22,"  t0_codeb: %s\n", t0_str[1]);
        trace(22,"  t0_clock: %s\n", t0_str[2]);
        trace(22,"  deph: %8.3f %8.3f %8.3f m\n", 
              nav->B2bssr[i].deph[0], nav->B2bssr[i].deph[1], nav->B2bssr[i].deph[2]);
        trace(22,"  dclk: %8.3f %8.3f %8.3f m\n",
              nav->B2bssr[i].dclk[0], nav->B2bssr[i].dclk[1], nav->B2bssr[i].dclk[2]);
    }

    /* Reset message counters */
    raw->num_PPPB2BINF01 = raw->num_PPPB2BINF02 = 0;
    raw->num_PPPB2BINF03 = raw->num_PPPB2BINF04 = 0;
}



/* update rtk server struct --------------------------------------------------*/
static void update_svr(rtksvr_t *svr, int ret, obs_t *obs, nav_t *nav,
                       int ephsat, int ephset, sbsmsg_t *sbsmsg, int index,
                       int iobs)
{
    tracet(4,"updatesvr: ret=%d ephsat=%d ephset=%d index=%d\n",ret,ephsat,
           ephset,index);
    if (ret==2) { /* ephemeris */
        update_eph(svr,nav,ephsat,ephset,index);
    }
    else if (ret==20) { /* B2bssr message */
        update_B2bssr(svr, index);
    }
    else if (ret==-1) { /* error */
        svr->nmsg[index][9]++;
    }
}

/* decode receiver raw/rtcm data ---------------------------------------------*/
static int decoderaw(rtksvr_t *svr, int index)
{
    obs_t *obs;
    nav_t *nav;
    sbsmsg_t *sbsmsg=NULL;
    int i,ret,ephsat,ephset,fobs=0;
    
    tracet(4,"decoderaw: index=%d\n",index);
    
    rtksvrlock(svr);
    
    for (i=0;i<svr->nb[index];i++) {
        
        /* input rtcm/receiver raw data from stream */
        if (svr->format[index]==STRFMT_UNICORE || svr->format[index]==STRFMT_SINO) {
			ret=input_raw(svr->raw+index,svr->format[index],svr->buff[index][i]);
            // obs=&svr->raw[index].obs;
            nav=&svr->raw[index].nav;
            ephsat=svr->raw[index].ephsat;
            ephset=svr->raw[index].ephset;
            // sbsmsg=&svr->raw[index].sbsmsg;
        }
        else {
            printf("Only support Unicore and Sino Format! \n");
            return 0;
        }
#if 0 /* record for receiving tick for debug */
        if (ret==1) {
            trace(0,"%d %10d T=%s NS=%2d\n",index,tickget(),
                  time_str(obs->data[0].time,0),obs->n);
        }
#endif
        /* update rtk server */
        if (ret>0) {
            update_svr(svr,ret,obs,nav,ephsat,ephset,sbsmsg,index,fobs);
        }
        /* observation data received */
        if (ret==1) {
            if (fobs<MAXOBSBUF) fobs++; else svr->prcout++;
        }
    }
    svr->nb[index]=0;
    
    rtksvrunlock(svr);
    
    return fobs;
}

/* rtk server thread ---------------------------------------------------------*/
#ifdef WIN32
static DWORD WINAPI decodesvrthread(void *arg)
#else
static void *decodesvrthread(void *arg)
#endif
{
    rtksvr_t *svr=(rtksvr_t *)arg;
    obs_t obs;
    obsd_t data[MAXOBS*2];
    sol_t sol={{0}};
    double tt;
    uint32_t tick,ticknmea,tick1hz,tickreset;
    uint8_t *p,*q;
    char msg[128];
    int i,j,n,cycle,cputime;
    int fobs[3]={0};

	double ep[6], dtt = 0;
	sol_t  sol1 = { { 0 } };
	int y1, mon1, day1, h1, doy1;
	int mjd2021 = 59215, pppnum = 30;
	char filename[20];
	FILE *fp;

    tracet(3,"decodesvrthread:\n");
    
    svr->state=1; obs.data=data;
    svr->tick=tickget();
    ticknmea=tick1hz=svr->tick-1000;
    tickreset=svr->tick-MIN_INT_RESET;
    
    for (cycle=0;svr->state;cycle++) {
        tick=tickget();
        for (i=0;i<3;i++) {
            p=svr->buff[i]+svr->nb[i]; q=svr->buff[i]+svr->buffsize;
            
            /* read receiver raw/rtcm data from input stream */
            if ((n=strread(svr->stream+i,p,q-p))<=0) {
                continue;
            }
            /* write receiver raw/rtcm data to log stream */
            strwrite(svr->stream+i+5,p,n);
            svr->nb[i]+=n;
            
            /* save peek buffer */
            rtksvrlock(svr);
            n=n<svr->buffsize-svr->npb[i]?n:svr->buffsize-svr->npb[i];
            memcpy(svr->pbuf[i]+svr->npb[i],p,n);
            svr->npb[i]+=n;
            rtksvrunlock(svr);
        }
		for (i=0;i<3;i++) {
            /* decode receiver raw/rtcm data */
            fobs[i]=decoderaw(svr,i);
            if (1==i&&svr->rtcm[1].staid>0) sol.refstationid=svr->rtcm[1].staid; 
            
        }
        /* averaging single base pos */
        // if (fobs[1]>0&&svr->rtk.opt.refpos==POSOPT_SINGLE) {
        //     if ((svr->rtk.opt.maxaveep<=0||svr->nave<svr->rtk.opt.maxaveep)&&
        //         pntpos(svr->obs[1][0].data,svr->obs[1][0].n,&svr->nav,
        //                &svr->rtk.opt,&sol,NULL,NULL,msg)) {
        //         svr->nave++;
        //         for (i=0;i<3;i++) {
        //             svr->rb_ave[i]+=(sol.rr[i]-svr->rb_ave[i])/svr->nave;
        //         }
        //     }
        //     for (i=0;i<3;i++) svr->rtk.opt.rb[i]=svr->rb_ave[i];
        // }
        // for (i=0;i<fobs[0];i++) { /* for each rover observation data */
        //     obs.n=0;
        //     for (j=0;j<svr->obs[0][i].n&&obs.n<MAXOBS*2;j++) {
        //         obs.data[obs.n++]=svr->obs[0][i].data[j];
        //     }
        //     for (j=0;j<svr->obs[1][0].n&&obs.n<MAXOBS*2;j++) {
        //         obs.data[obs.n++]=svr->obs[1][0].data[j];
		// 	}
        //     /* if cpu overload, increment obs outage counter and break */
        //     if ((int)(tickget()-tick)>=svr->cycle) {
        //         svr->prcout+=fobs[0]-i-1;
        //     }
        // }
        
        if ((cputime=(int)(tickget()-tick))>0) svr->cputime=cputime;
        
        /* sleep until next cycle */
        sleepms(svr->cycle-cputime);
    }
    for (i=0;i<MAXSTRRTK;i++) strclose(svr->stream+i);
    for (i=0;i<3;i++) {
        svr->nb[i]=svr->npb[i]=0;
        free(svr->buff[i]); svr->buff[i]=NULL;
        free(svr->pbuf[i]); svr->pbuf[i]=NULL;
        free_raw (svr->raw +i);
        // free_rtcm(svr->rtcm+i);
    }
    for (i=0;i<2;i++) {
        svr->nsb[i]=0;
        free(svr->sbuf[i]); svr->sbuf[i]=NULL;
    }
    return 0;
}



static int decodesvrstart(rtksvr_t *svr, int cycle, int buffsize, int *strs,
                       const char **paths, int *formats, int navsel, const char **cmds,
                       const char **cmds_periodic, const char **rcvopts, int nmeacycle,
                       int nmeareq, const double *nmeapos, prcopt_t *prcopt,
                       solopt_t *solopt, stream_t *moni, char *errmsg)
{
    gtime_t time,time0={0};
    int i,j,rw;
    
    tracet(3,"decodesvrstart: cycle=%d buffsize=%d navsel=%d nmeacycle=%d nmeareq=%d\n",
           cycle,buffsize,navsel,nmeacycle,nmeareq);
    
    if (svr->state) {
        sprintf(errmsg,"server already started");
        return 0;
    }
    strinitcom();
    svr->cycle=cycle>1?cycle:1;
    svr->nmeacycle=nmeacycle>1000?nmeacycle:1000;
    svr->nmeareq=nmeareq;
    for (i=0;i<3;i++) svr->nmeapos[i]=nmeapos[i];
    svr->buffsize=buffsize>4096?buffsize:4096;
    for (i=0;i<3;i++) svr->format[i]=formats[i];
    svr->navsel=navsel;
    svr->nsbs=0;
    svr->nsol=0;
    svr->prcout=0;
    rtkfree(&svr->rtk);
    rtkinit(&svr->rtk,prcopt);
    
    for (i=0;i<3;i++) { /* input/log streams */
        svr->nb[i]=svr->npb[i]=0;
        if (!(svr->buff[i]=(uint8_t *)malloc(buffsize))||
            !(svr->pbuf[i]=(uint8_t *)malloc(buffsize))) {
            tracet(1,"rtksvrstart: malloc error\n");
            sprintf(errmsg,"rtk server malloc error");
            return 0;
        }
        for (j=0;j<10;j++) svr->nmsg[i][j]=0;
        for (j=0;j<MAXOBSBUF;j++) svr->obs[i][j].n=0;
        // strcpy(svr->cmds_periodic[i],!cmds_periodic[i]?"":cmds_periodic[i]);
        
        /* initialize receiver raw and rtcm control */
        init_raw(svr->raw+i,formats[i]);
        // init_rtcm(svr->rtcm+i);
        
        /* set receiver and rtcm option */
        strcpy(svr->raw [i].opt,rcvopts[i]);
        // strcpy(svr->rtcm[i].opt,rcvopts[i]);
        
        /* connect dgps corrections */
        // svr->rtcm[i].dgps=svr->nav.dgps;
    }
    for (i=0;i<2;i++) { /* output peek buffer */
        if (!(svr->sbuf[i]=(uint8_t *)malloc(buffsize))) {
            tracet(1,"decodesvrstart: malloc error\n");
            sprintf(errmsg,"decode server malloc error");
            return 0;
        }
    }
    /* set solution options */
    for (i=0;i<2;i++) {
        svr->solopt[i]=solopt[i];
    }

    /* update navigation data */
    for (i=0;i<MAXSAT*4 ;i++) svr->nav.eph [i].ttr=time0;
    for (i=0;i<NSATGLO*2;i++) svr->nav.geph[i].tof=time0;
    for (i=0;i<NSATSBS*2;i++) svr->nav.seph[i].tof=time0;

    
    /* open input streams */
    for (i=0;i<8;i++) {
        rw=i<3?STR_MODE_R:STR_MODE_W;
        if (strs[i]!=STR_FILE) rw|=STR_MODE_W;
        if (!stropen(svr->stream+i,strs[i],rw,paths[i])) {
            sprintf(errmsg,"str%d open error path=%s",i+1,paths[i]);
            for (i--;i>=0;i--) strclose(svr->stream+i);
            return 0;
        }
        /* set initial time for rtcm and raw */
        if (i<3) {
            time=utc2gpst(timeget());
            svr->raw [i].time=strs[i]==STR_FILE?strgettime(svr->stream+i):time;
            svr->rtcm[i].time=strs[i]==STR_FILE?strgettime(svr->stream+i):time;
        }
    }
    /* sync input streams */
    strsync(svr->stream,svr->stream+1);
    strsync(svr->stream,svr->stream+2);
    
    /* write start commands to input streams */
    for (i=0;i<3;i++) {
        if (!cmds[i]) continue;
        strwrite(svr->stream+i,(unsigned char *)"",0); /* for connect */
        sleepms(100);
        strsendcmd(svr->stream+i,cmds[i]);
    }

    /* create rtk server thread */
#ifdef WIN32
    if (!(svr->thread=CreateThread(NULL,0,decodesvrthread,svr,0,NULL))) {
#else
    if (pthread_create(&svr->thread,NULL,decodesvrthread,svr)) {
#endif
        for (i=0;i<MAXSTRRTK;i++) strclose(svr->stream+i);
        sprintf(errmsg,"thread create error\n");
        return 0;
    }
    return 1;
}

/* print usage ---------------------------------------------------------------*/
static void printusage(void)
{
    int i;
    for (i=0;i<(int)(sizeof(usage)/sizeof(*usage));i++) {
        fprintf(stderr,"%s\n",usage[i]);
    }
    exit(0);
}
/* start rtdecoder server ----------------------------------------------------------*/
static int startsvr(vt_t *vt)
{
    static sta_t sta[MAXRCV]={{""}};
    double pos[3],npos[3];
    char s1[3][MAXRCVCMD]={"","",""},*cmds[]={NULL,NULL,NULL};
    char s2[3][MAXRCVCMD]={"","",""},*cmds_periodic[]={NULL,NULL,NULL};
    char *ropts[]={rcvopt[0],rcvopt[1],rcvopt[2]};
    char *paths[]={
        strpath[0],strpath[1],strpath[2],strpath[3],strpath[4],strpath[5],
        strpath[6],strpath[7]
    };
    char errmsg[2048]="";
    int i,stropt[8]={0};
    
    trace(3,"startdecode:\n");
    

    for (i=0;*rcvopts[i].name;i++) modflgr[i]=0;
    for (i=0;*sysopts[i].name;i++) modflgs[i]=0;
    
    /* set stream options */
    stropt[0]=timeout;
    stropt[1]=reconnect;
    stropt[2]=1000;
    stropt[3]=buffsize;
    stropt[4]=fswapmargin;
    strsetopt(stropt);
    
    if (strfmt[2]==8) strfmt[2]=STRFMT_SP3;
    
    /* set ftp/http directory and proxy */
    strsetdir(filopt.tempdir);
    strsetproxy(proxyaddr);
    
    solopt[0].posf=strfmt[3];
    solopt[1].posf=strfmt[4];
    
    /* start decode server */
    if (!decodesvrstart(&svr,svrcycle,buffsize,strtype,(const char **)paths,strfmt,navmsgsel,
                     (const char **)cmds,(const char **)cmds_periodic,(const char **)ropts,nmeacycle,nmeareq,npos,&prcopt,
                     solopt,&moni,errmsg)) {
        trace(2,"decode server start error (%s)\n",errmsg);
        printf("decode server start error (%s)\n",errmsg);
        return 0;
    }
    return 1;
}

/* stop decode server -----------------------------------------------------------*/
static void stopsvr()
{
    char s[3][MAXRCVCMD]={"","",""},*cmds[]={NULL,NULL,NULL};
    int i;
    
    trace(3,"stopsvr:\n");
    
    if (!svr.state) return;
    /* stop decode server */
    rtksvrstop(&svr,(const char **)cmds);
    
    
    printf("stop decode server\n");
}

/* print time ----------------------------------------------------------------*/
static void prtime(vt_t *vt, gtime_t time)
{
    double tow;
    int week;
    char tstr[64]="";
    
    if (timetype==1) {
        time2str(gpst2utc(time),tstr,2);
    }
    else if (timetype==2) {
        time2str(timeadd(gpst2utc(time),9*3600.0),tstr,2);
    }
    else if (timetype==3) {
        tow=time2gpst(time,&week); sprintf(tstr,"  %04d %9.2f",week,tow);
    }
    else time2str(time,tstr,1);
    vt_printf(vt,"%s ",tstr);
}
/* print solution ------------------------------------------------------------*/
static void prsolution(vt_t *vt, const sol_t *sol, const double *rb)
{
    const char *solstr[]={"------","FIX","FLOAT","SBAS","DGPS","SINGLE","PPP",""};
    double pos[3]={0},Qr[9],Qe[9]={0},dms1[3]={0},dms2[3]={0},bl[3]={0};
    double enu[3]={0},pitch=0.0,yaw=0.0,len;
    int i;
    
    trace(4,"prsolution:\n");
    
    if (sol->time.time==0||!sol->stat) return;
    prtime(vt,sol->time);
    vt_printf(vt,"(%-6s)",solstr[sol->stat]);
    
    if (norm(sol->rr,3)>0.0&&norm(rb,3)>0.0) {
        for (i=0;i<3;i++) bl[i]=sol->rr[i]-rb[i];
    }
    len=norm(bl,3);
    Qr[0]=sol->qr[0];
    Qr[4]=sol->qr[1];
    Qr[8]=sol->qr[2];
    Qr[1]=Qr[3]=sol->qr[3];
    Qr[5]=Qr[7]=sol->qr[4];
    Qr[2]=Qr[6]=sol->qr[5];
    
    if (soltype==0) {
        if (norm(sol->rr,3)>0.0) {
            ecef2pos(sol->rr,pos);
            covenu(pos,Qr,Qe);
            deg2dms(pos[0]*R2D,dms1,4);
            deg2dms(pos[1]*R2D,dms2,4);
            if (solopt[0].height==1) pos[2]-=geoidh(pos); /* geodetic */
        }       
        vt_printf(vt," %s:%2.0f %02.0f %07.4f",pos[0]<0?"S":"N",fabs(dms1[0]),dms1[1],dms1[2]);
        vt_printf(vt," %s:%3.0f %02.0f %07.4f",pos[1]<0?"W":"E",fabs(dms2[0]),dms2[1],dms2[2]);
        vt_printf(vt," H:%8.3f",pos[2]);
        if (solflag&1) {
            vt_printf(vt," (N:%6.3f E:%6.3f U:%6.3f)",SQRT(Qe[4]),SQRT(Qe[0]),SQRT(Qe[8]));
        }
    }
    else if (soltype==1) {
        if (norm(sol->rr,3)>0.0) {
            ecef2pos(sol->rr,pos);
            covenu(pos,Qr,Qe);
            if (solopt[0].height==1) pos[2]-=geoidh(pos); /* geodetic */
        }       
        vt_printf(vt," %s:%11.8f",pos[0]<0.0?"S":"N",fabs(pos[0])*R2D);
        vt_printf(vt," %s:%12.8f",pos[1]<0.0?"W":"E",fabs(pos[1])*R2D);
        vt_printf(vt," H:%8.3f",pos[2]);
        if (solflag&1) {
            vt_printf(vt," (E:%6.3f N:%6.3f U:%6.3fm)",SQRT(Qe[0]),SQRT(Qe[4]),SQRT(Qe[8]));
        }
    }
    else if (soltype==2) {
        vt_printf(vt," X:%12.3f",sol->rr[0]);
        vt_printf(vt," Y:%12.3f",sol->rr[1]);
        vt_printf(vt," Z:%12.3f",sol->rr[2]);
        if (solflag&1) {
            vt_printf(vt," (X:%6.3f Y:%6.3f Z:%6.3f)",SQRT(Qr[0]),SQRT(Qr[4]),SQRT(Qr[8]));
        }
    }
    else if (soltype==3) {
        if (len>0.0) {
            ecef2pos(rb,pos);
            ecef2enu(pos,bl,enu);
            covenu(pos,Qr,Qe);
        }       
        vt_printf(vt," E:%12.3f",enu[0]);
        vt_printf(vt," N:%12.3f",enu[1]);
        vt_printf(vt," U:%12.3f",enu[2]);
        if (solflag&1) {
            vt_printf(vt," (E:%6.3f N:%6.3f U:%6.3f)",SQRT(Qe[0]),SQRT(Qe[4]),SQRT(Qe[8]));
        }
    }
    else if (soltype==4) {
        if (len>0.0) {
            ecef2pos(rb,pos);
            ecef2enu(pos,bl,enu);
            covenu(pos,Qr,Qe);
            pitch=asin(enu[2]/len);
            yaw=atan2(enu[0],enu[1]); if (yaw<0.0) yaw+=2.0*PI;
        }
        vt_printf(vt," P:%12.3f",pitch*R2D);
        vt_printf(vt," Y:%12.3f",yaw*R2D);
        vt_printf(vt," L:%12.3f",len);
        if (solflag&1) {
            vt_printf(vt," (E:%6.3f N:%6.3f U:%6.3f)",SQRT(Qe[0]),SQRT(Qe[4]),SQRT(Qe[8]));
        }
    }
    if (solflag&2) {
        vt_printf(vt," A:%4.1f R:%5.1f N:%2d",sol->age,sol->ratio,sol->ns);
    }
    vt_printf(vt,"\n");
}
/* print status --------------------------------------------------------------*/
static void prstatus(vt_t *vt)
{
    rtk_t rtk;
    const char *svrstate[]={"stop","run"},*type[]={"stream_1","stream_2",""};
    const char *sol[]={"-","fix","float","SBAS","DGPS","single","PPP-FLOAT",""};
    const char *mode[]={
         "single","DGPS","kinematic","static","static-start","moving-base","fixed",
         "PPP-kinema","PPP-static"
    };

    const char *strfmt_str[] = {
    "rtcm2", "rtcm3", "oem4",    "", "ubx", "sbp", "cres", "stq", "javad", "nvs", "binex",
    "rt17",   "sept",     "",    "", "rinex", "sp3", "rnxclk", "sbas", "nmea", "sino", "unicore"
    };

    const char *ephopt_str[] = {
        "brdc", "prec", "sbas", "brdc+ssrapc", "brdc+ssrcom", "brdc+b2b"
    };
    // STRFMT_NMEA
    gtime_t eventime={0};
    const char *freq[]={"-","L1","L1+L2","L1+L2+E5b","L1+L2+E5b+L5","",""};
    rtcm_t rtcm[3];
    pthread_t thread;
    int i,j,n,cycle,state,rtkstat,nsat0,nsat1,prcout,rcvcount,tmcount,timevalid,nave;
    int cputime,nb[3]={0},nmsg[3][10]={{0}},raw_nmsg[2][16]={0};
    char tstr[64],tmstr[64],s[1024],*p;
    double runtime,rt[3]={0},dop[4]={0},rr[3],bl1=0.0,bl2=0.0;
    double azel[MAXSAT*2],pos[3],vel[3],*del;
    
    trace(4,"prstatus:\n");
    
    rtksvrlock(&svr);
    rtk=svr.rtk;
    thread=svr.thread;
    cycle=svr.cycle;
    state=svr.state;
    rtkstat=svr.rtk.sol.stat;
    nsat0=svr.obs[0][0].n;
    nsat1=svr.obs[1][0].n;
    rcvcount = svr.raw[0].obs.rcvcount;
    tmcount = svr.raw[0].obs.tmcount;
    cputime=svr.cputime;
    prcout=svr.prcout;
    nave=svr.nave;
    for (i=0;i<3;i++) nb[i]=svr.nb[i];
    for (i=0;i<3;i++) for (j=0;j<10;j++) {
        nmsg[i][j]=svr.nmsg[i][j];
    }
    for(i=0;i<16;i++){
        raw_nmsg[0][i] = svr.raw[0].raw_nmsg[i];   //base
        raw_nmsg[1][i] = svr.raw[1].raw_nmsg[i];   //corr
    }
    
    if (svr.state) {
        runtime=(double)(tickget()-svr.tick)/1000.0;
        rt[0]=floor(runtime/3600.0); runtime-=rt[0]*3600.0;
        rt[1]=floor(runtime/60.0); rt[2]=runtime-rt[1]*60.0;
    }
    for (i=0;i<3;i++) rtcm[i]=svr.rtcm[i];
    if (svr.raw[0].obs.data != NULL) {
        timevalid = svr.raw[0].obs.data[0].timevalid;
        eventime = svr.raw[0].obs.data[0].eventime;
    }
    time2str(eventime,tmstr,9);
    rtksvrunlock(&svr);
    
    for (i=n=0;i<MAXSAT;i++) {
        if (rtk.opt.mode==PMODE_SINGLE&&!rtk.ssat[i].vs) continue;
        if (rtk.opt.mode!=PMODE_SINGLE&&!rtk.ssat[i].vsat[0]) continue;
        azel[  n*2]=rtk.ssat[i].azel[0];
        azel[1+n*2]=rtk.ssat[i].azel[1];
        n++;
    }
    dops(n,azel,0.0,dop);

    // const char *rover_fmt = (strfmt[0] >= 0 && strfmt[0] < 22) ? strfmt_str[strfmt[0]] : "Unknown";
    const char *base_fmt = (strfmt[1] >= 0 && strfmt[1] < 22) ? strfmt_str[strfmt[0]] : "Unknown";
    const char *corr_fmt = (strfmt[2] >= 0 && strfmt[2] < 22) ? strfmt_str[strfmt[1]] : "Unknown";
    const char *ephopt = (rtk.opt.sateph >= 0 && rtk.opt.sateph <= 5) ? ephopt_str[rtk.opt.sateph] : "Unknown";
    
    vt_printf(vt,"\n%s%-28s: %s%s\n",ESC_BOLD,"Parameter","Value",ESC_RESET);
    vt_printf(vt,"%-28s: %s %s\n","B2bPPP version",VER_RTKLIB,PATCH_LEVEL);
    vt_printf(vt,"%-28s: %d\n","B2bPPP server thread",thread);
    vt_printf(vt,"%-28s: %s\n","B2bPPP server state",svrstate[state]);
    vt_printf(vt,"%-28s: %d\n","processing cycle (ms)",cycle);
    // vt_printf(vt,"%-28s: %s\n","positioning mode",mode[rtk.opt.mode]);
    // vt_printf(vt,"%-28s: %s\n","satpositioning ephopt:",ephopt);
    // vt_printf(vt,"%-28s: %s\n","frequencies",freq[rtk.opt.nf]);
    // vt_printf(vt,"%-28s: %02.0f:%02.0f:%04.1f\n","accumulated time to run",rt[0],rt[1],rt[2]);
    // vt_printf(vt,"%-28s: %d\n","cpu time for a cycle (ms)",cputime);
    // vt_printf(vt,"%-28s: %d\n","missing obs data count",prcout);
    // vt_printf(vt,"%-28s: %d,%d\n","bytes in input buffer",nb[0],nb[1]);    
    // vt_printf(vt,"%-28s: %s\n","Base stream format",base_fmt);
    // vt_printf(vt,"%-28s: %s\n","Corr stream format",corr_fmt);
    for (i=0;i<2;i++) {
        sprintf(s,"# of input data %s",type[i]);
        vt_printf(vt,"%-28s: nav(%d),B2bssr(%d),err(%d)\n",
                s,nmsg[i][1],nmsg[i][8],nmsg[i][9]);}

    for (i = 0; i < 2; i++) {
        if (i == 0) {
            // 添加 base 的数据流格式
            char param_name[64];
            sprintf(param_name, "# of nav messages (%s)", base_fmt);
            vt_printf(vt, "%-28s: gps_lnav(%d), bds_cnav1(%d)\n",
                    param_name, raw_nmsg[i][8], raw_nmsg[i][9]);
        }
        if (i == 1) {
            // 添加 corr 的数据流格式
            char param_name[64];
            sprintf(param_name, "# of B2bssr messages (%s)", corr_fmt);
            vt_printf(vt, "%-28s: mask(%d), eph_ura(%d), diff_cbia(%d), clk(%d), others(%d)\n",
                    param_name, raw_nmsg[i][0], raw_nmsg[i][1], raw_nmsg[i][2], raw_nmsg[i][3], raw_nmsg[i][7]);
        }
    }
}
/* print satellite -----------------------------------------------------------*/
static void prsatellite(vt_t *vt, int nf)
{
    rtk_t rtk;
    double az,el;
    char id[32];
    int i,j,fix,frq[]={1,2,5,7,8,6};
    
    trace(4,"prsatellite:\n");
    
    rtksvrlock(&svr);
    rtk=svr.rtk;
    rtksvrunlock(&svr);
    if (nf<=0||nf>NFREQ) nf=NFREQ;
    vt_printf(vt,"\n%s%3s %2s %5s %4s",ESC_BOLD,"SAT","C1","Az","El");
    for (j=0;j<nf;j++) vt_printf(vt," L%d"    ,frq[j]);
    for (j=0;j<nf;j++) vt_printf(vt,"  Fix%d" ,frq[j]);
    for (j=0;j<nf;j++) vt_printf(vt,"  P%dRes",frq[j]);
    for (j=0;j<nf;j++) vt_printf(vt,"   L%dRes",frq[j]);
    for (j=0;j<nf;j++) vt_printf(vt,"  Sl%d"  ,frq[j]);
    for (j=0;j<nf;j++) vt_printf(vt,"  Lock%d",frq[j]);
    for (j=0;j<nf;j++) vt_printf(vt," Rj%d"   ,frq[j]);
    vt_printf(vt,"%s\n",ESC_RESET);
    
    for (i=0;i<MAXSAT;i++) {
        if (rtk.ssat[i].azel[1]<=0.0) continue;
        satno2id(i+1,id);
        vt_printf(vt,"%3s %2s",id,rtk.ssat[i].vs?"OK":"-");
        az=rtk.ssat[i].azel[0]*R2D; if (az<0.0) az+=360.0;
        el=rtk.ssat[i].azel[1]*R2D;
        vt_printf(vt," %5.1f %4.1f",az,el);
        for (j=0;j<nf;j++) vt_printf(vt," %2s",rtk.ssat[i].vsat[j]?"OK":"-");
        for (j=0;j<nf;j++) {
            fix=rtk.ssat[i].fix[j];
            vt_printf(vt," %5s",fix==1?"FLOAT":(fix==2?"FIX":(fix==3?"HOLD":"-")));
        }
        for (j=0;j<nf;j++) vt_printf(vt,"%7.3f",rtk.ssat[i].resp[j]);
        for (j=0;j<nf;j++) vt_printf(vt,"%8.4f",rtk.ssat[i].resc[j]);
        for (j=0;j<nf;j++) vt_printf(vt," %4d",rtk.ssat[i].slipc[j]);
        for (j=0;j<nf;j++) vt_printf(vt," %6d",rtk.ssat[i].lock [j]);
        for (j=0;j<nf;j++) vt_printf(vt," %3d",rtk.ssat[i].rejc [j]);
        vt_printf(vt,"\n");
    }
}
/* print observation data ----------------------------------------------------*/
static void probserv(vt_t *vt, int nf)
{
    obsd_t obs[MAXOBS*2];
    char tstr[64],id[32];
    int i,j,n=0,frq[]={1,2,5,7,8,6,9};
    
    trace(4,"probserv:\n");
    
    rtksvrlock(&svr);
    for (i=0;i<svr.obs[0][0].n&&n<MAXOBS*2;i++) {
        obs[n++]=svr.obs[0][0].data[i];
    }
    for (i=0;i<svr.obs[1][0].n&&n<MAXOBS*2;i++) {
        obs[n++]=svr.obs[1][0].data[i];
    }
    rtksvrunlock(&svr);
    
    if (nf<=0||nf>NFREQ) nf=NFREQ;
    vt_printf(vt,"\n%s%-22s %3s %s",ESC_BOLD,"      TIME(GPST)","SAT","R");
    for (i=0;i<nf;i++) vt_printf(vt,"        P%d(m)" ,frq[i]);
    for (i=0;i<nf;i++) vt_printf(vt,"       L%d(cyc)",frq[i]);
    for (i=0;i<nf;i++) vt_printf(vt,"  D%d(Hz)"      ,frq[i]);
    for (i=0;i<nf;i++) vt_printf(vt," S%d"           ,frq[i]);
    vt_printf(vt," LLI%s\n",ESC_RESET);
    for (i=0;i<n;i++) {
        time2str(obs[i].time,tstr,2);
        satno2id(obs[i].sat,id);
        vt_printf(vt,"%s %3s %d",tstr,id,obs[i].rcv);
        for (j=0;j<nf;j++) vt_printf(vt,"%13.3f",obs[i].P[j]);
        for (j=0;j<nf;j++) vt_printf(vt,"%14.3f",obs[i].L[j]);
        for (j=0;j<nf;j++) vt_printf(vt,"%8.1f" ,obs[i].D[j]);
        for (j=0;j<nf;j++) vt_printf(vt,"%3.0f" ,obs[i].SNR[j]*SNR_UNIT);
        for (j=0;j<nf;j++) vt_printf(vt,"%2d"   ,obs[i].LLI[j]);
        vt_printf(vt,"\n");
    }
}
/* print navigation data -----------------------------------------------------*/
static void prnavidata(vt_t *vt)
{
    eph_t eph[MAXSAT];
    geph_t geph[MAXPRNGLO];
    double ion[8],utc[8];
    gtime_t time;
    char id[32],s1[64],s2[64],s3[64];
    int i,valid,prn;
    
    trace(4,"prnavidata:\n");
    
    rtksvrlock(&svr);
    // time=svr.rtk.sol.time;
    time=svr.raw[0].time;
    for (i=0;i<MAXSAT;i++) eph[i]=svr.nav.eph[i];
    for (i=0;i<MAXPRNGLO;i++) geph[i]=svr.nav.geph[i];
    for (i=0;i<8;i++) ion[i]=svr.nav.ion_gps[i];
    for (i=0;i<8;i++) utc[i]=svr.nav.utc_gps[i];
    rtksvrunlock(&svr);
    
    vt_printf(vt,"\n%s%3s %3s %3s %3s %3s %3s %3s %19s %19s %19s %3s %3s%s\n",
              ESC_BOLD,"SAT","S","IOD","IOC","FRQ","A/A","SVH","Toe","Toc",
              "Ttr/Tof","L2C","L2P",ESC_RESET);
    for (i=0;i<MAXSAT;i++) {
        if (!(satsys(i+1,&prn)&(SYS_GPS|SYS_GAL|SYS_QZS|SYS_CMP))||
            eph[i].sat!=i+1) continue;
        valid=eph[i].toe.time!=0&&!eph[i].svh&&
              fabs(timediff(time,eph[i].toe))<=MAXDTOE;
        satno2id(i+1,id);
        if (eph[i].toe.time!=0) time2str(eph[i].toe,s1,0); else strcpy(s1,"-");
        if (eph[i].toc.time!=0) time2str(eph[i].toc,s2,0); else strcpy(s2,"-");
        if (eph[i].ttr.time!=0) time2str(eph[i].ttr,s3,0); else strcpy(s3,"-");
        vt_printf(vt,"%3s %3s %3d %3d %3d %3d %03X %19s %19s %19s %3d %3d\n",
                id,valid?"OK":"-",eph[i].iode,eph[i].iodc,0,eph[i].sva,
                eph[i].svh,s1,s2,s3,eph[i].code,eph[i].flag);
    }
    for (i=0;i<MAXSAT;i++) {
        if (!(satsys(i+1,&prn)&SYS_GLO)||geph[prn-1].sat!=i+1) continue;
        valid=geph[prn-1].toe.time!=0&&!geph[prn-1].svh&&
              fabs(timediff(time,geph[prn-1].toe))<=MAXDTOE_GLO;
        satno2id(i+1,id);
        if (geph[prn-1].toe.time!=0) time2str(geph[prn-1].toe,s1,0); else strcpy(s1,"-");
        if (geph[prn-1].tof.time!=0) time2str(geph[prn-1].tof,s2,0); else strcpy(s2,"-");
        vt_printf(vt,"%3s %3s %3d %3d %3d %3d  %02X %19s %19s %19s %3d %3d\n",
                id,valid?"OK":"-",geph[prn-1].iode,0,geph[prn-1].frq,
                geph[prn-1].age,geph[prn].svh,s1,"-",s2,0,0);
    }
    // vt_printf(vt,"ION: %9.2E %9.2E %9.2E %9.2E %9.2E %9.2E %9.2E %9.2E\n",
    //         ion[0],ion[1],ion[2],ion[3],ion[4],ion[5],ion[6],ion[7]);
    // vt_printf(vt,"UTC: %9.2E %9.2E %9.2E %9.2E  LEAPS: %.0f\n",utc[0],utc[1],utc[2],
    //         utc[3],utc[4]);
}
/* print error/warning messages ----------------------------------------------*/
static void prerror(vt_t *vt)
{
    int n;
    
    trace(4,"prerror:\n");
    
    rtksvrlock(&svr);
    if ((n=svr.rtk.neb)>0) {
        svr.rtk.errbuf[n]='\0';
        vt_puts(vt,svr.rtk.errbuf);
        svr.rtk.neb=0;
    }
    rtksvrunlock(&svr);
}
/* print stream --------------------------------------------------------------*/
static void prstream(vt_t *vt)
{
    const char *ch[]={
        "input stream1","input stream2","input corr","output sol1","output sol2",
        "log stream1","log stream2","log corr","monitor"
    };
    const char *type[]={
        "-","serial","file","tcpsvr","tcpcli","ntrips","ntripc","ftp",
        "http","ntripcas","udpsvr","udpcli","membuf"
    };
    const char *fmt[]={
    "rtcm2", "rtcm3", "oem4",    "", "ubx", "sbp", "cres", "stq", "javad", "nvs", "binex",
    "rt17",   "sept",     "",    "", "rinex", "sp3", "rnxclk", "sbas", "nmea", "sino", "unicore"
    };
    const char *sol[]={"llh","xyz","enu","nmea","stat","-"};
    stream_t stream[9];
    int i,format[9]={0};
    
    trace(4,"prstream:\n");
    
    rtksvrlock(&svr);
    for (i=0;i<8;i++) stream[i]=svr.stream[i];
    for (i=0;i<3;i++) format[i]=svr.format[i];
    for (i=3;i<5;i++) format[i]=svr.solopt[i-3].posf;
    stream[8]=moni;
    format[8]=SOLF_LLH;
    rtksvrunlock(&svr);
    
    vt_printf(vt,"\n%s%-12s %-8s %-5s %s %10s %7s %10s %7s %-24s %s%s\n",ESC_BOLD,
              "Stream","Type","Fmt","S","In-byte","In-bps","Out-byte","Out-bps",
              "Path","Message",ESC_RESET);
    for (i=0;i<9;i++) {
        if(i==2||i==3||i==4||i==7||i==8) continue;
        vt_printf(vt,"%-12s %-8s %-5s %s %10d %7d %10d %7d %-24.24s %s\n",
            ch[i],type[stream[i].type],i<3?fmt[format[i]]:(i<5||i==8?sol[format[i]]:"-"),
            stream[i].state<0?"E":(stream[i].state?"C":"-"),
            stream[i].inb,stream[i].inr,stream[i].outb,stream[i].outr,
            stream[i].path,stream[i].msg);
    }
}
/* print ssr correction ------------------------------------------------------*/
/* print B2bSSR correction ---------------------------------------------------*/
static void prssr(vt_t *vt)
{
    static char buff[128*MAXSAT];
    gtime_t time;
    B2bssr_t b2bssr[MAXSAT];
    int i, valid;
    char tstr[64], id[32], *p = buff;
    
    /* Retrieve B2bSSR data */
    rtksvrlock(&svr);
    time = svr.raw[1].time;  // Time when B2b was last received
    for (i = 0; i < MAXSAT; i++) {
        b2bssr[i] = svr.nav.B2bssr[i];
    }
    rtksvrunlock(&svr);
    
    /* Output header */
    p += sprintf(p, "\n%s%3s %3s %3s %3s %3s %19s %6s %6s %6s %8s %6s %6s%s\n",
                 ESC_BOLD, "SAT", "S", "UDI", "IODSSR", "URA", "T0", "D0-A", "D0-C", "D0-R",
                 "C0", "C1", "C2", ESC_RESET);
    
    /* Traverse satellites and output correction information */
    for (i = 1; i < MAXSAT; i++) {
        if (!b2bssr[i].t0[0].time) continue;  /* Check if eph type timestamp exists */
        satno2id(i, id);  /* Get satellite ID */
        valid = fabs(timediff(time, b2bssr[i].t0[0])) <= 1800.0;  /* Determine validity */
        time2str(b2bssr[i].t0[0], tstr, 0);  /* Format timestamp */
        
        p += sprintf(p, "%3s %3s %3.0f %3d %6d %19s %6.3f %6.3f %6.3f %8.3f %6.3f %6.4f\n",
                     id, valid ? "OK" : "-",                  /* Satellite ID and status */
                     b2bssr[i].udi[0],                        /* Update interval for eph */
                     b2bssr[i].iodssr[0],                     /* IODSSR for eph */
                     b2bssr[i].ura,                           /* URA index */
                     tstr,                                    /* Timestamp */
                     b2bssr[i].deph[0], b2bssr[i].deph[1], b2bssr[i].deph[2],    /* Orbit correction (m) */
                     b2bssr[i].dclk[0],                       /* Clock correction c0 (m) */
                     b2bssr[i].dclk[1] * 1E3,                 /* Clock correction c1 (mm/s) */
                     b2bssr[i].dclk[2] * 1E3);                /* Clock correction c2 (mm/s^2) */
    }
    vt_puts(vt, buff);  /* Output to terminal */
}




/* start command -------------------------------------------------------------*/
static void cmd_start(char **args, int narg, vt_t *vt)
{
    trace(3,"cmd_start:\n");
    
    if (!startsvr(vt)) return;
    vt_printf(vt,"rtk server start\n");
}
/* stop command --------------------------------------------------------------*/
static void cmd_stop(char **args, int narg, vt_t *vt)
{
    trace(3,"cmd_stop:\n");
    
    stopsvr(vt);
    vt_printf(vt,"rtk server stop\n");
}
/* restart command -----------------------------------------------------------*/
static void cmd_restart(char **args, int narg, vt_t *vt)
{
    trace(3,"cmd_restart:\n");
    
    stopsvr(vt);
    if (!startsvr(vt)) return;
    vt_printf(vt,"rtk server restart\n");
}
/* solution command ----------------------------------------------------------*/
static void cmd_solution(char **args, int narg, vt_t *vt)
{
    int i,cycle=0;
    
    trace(3,"cmd_solution:\n");
    
    if (narg>1) cycle=(int)(atof(args[1])*1000.0);
    
    if (cycle>0) svr.nsol=0;
    
    while (!vt_chkbrk(vt)) {
        rtksvrlock(&svr);
        for (i=0;i<svr.nsol;i++) prsolution(vt,&svr.solbuf[i],svr.rtk.rb);
        svr.nsol=0;
        rtksvrunlock(&svr);
        if (cycle>0) sleepms(cycle); else return;
    }
}
/* status command ------------------------------------------------------------*/
static void cmd_status(char **args, int narg, vt_t *vt)
{
    int cycle=0;
    
    trace(3,"cmd_status:\n");
    
    if (narg>1) cycle=(int)(atof(args[1])*1000.0);
    
    while (!vt_chkbrk(vt)) {
        if (cycle>0) vt_printf(vt,ESC_CLEAR);
        prstatus(vt);
        if (cycle>0) sleepms(cycle); else return;
    }
    vt_printf(vt,"\n");
}
/* satellite command ---------------------------------------------------------*/
static void cmd_satellite(char **args, int narg, vt_t *vt)
{
    int i,nf=2,cycle=0;
    
    trace(3,"cmd_satellite:\n");
    
    for (i=1;i<narg;i++) {
        if (sscanf(args[i],"-%d",&nf)<1) cycle=(int)(atof(args[i])*1000.0);
    }
    while (!vt_chkbrk(vt)) {
        if (cycle>0) vt_printf(vt,ESC_CLEAR);
        prsatellite(vt,nf);
        if (cycle>0) sleepms(cycle); else return;
    }
    vt_printf(vt,"\n");
}
/* observ command ------------------------------------------------------------*/
static void cmd_observ(char **args, int narg, vt_t *vt)
{
    int i,nf=2,cycle=0;
    
    trace(3,"cmd_observ:\n");
    
    for (i=1;i<narg;i++) {
        if (sscanf(args[i],"-%d",&nf)<1) cycle=(int)(atof(args[i])*1000.0);
    }
    while (!vt_chkbrk(vt)) {
        if (cycle>0) vt_printf(vt,ESC_CLEAR);
        probserv(vt,nf);
        if (cycle>0) sleepms(cycle); else return;
    }
    vt_printf(vt,"\n");
}
/* navidata command ----------------------------------------------------------*/
static void cmd_navidata(char **args, int narg, vt_t *vt)
{
    int cycle=0;
    
    trace(3,"cmd_navidata:\n");
    
    if (narg>1) cycle=(int)(atof(args[1])*1000.0);
    
    while (!vt_chkbrk(vt)) {
        if (cycle>0) vt_printf(vt,ESC_CLEAR);
        prnavidata(vt);
        if (cycle>0) sleepms(cycle); else return;
    }
    vt_printf(vt,"\n");
}
/* error command -------------------------------------------------------------*/
static void cmd_error(char **args, int narg, vt_t *vt)
{
    trace(3,"cmd_error:\n");
    
    rtksvrlock(&svr);
    svr.rtk.neb=0;
    rtksvrunlock(&svr);
    
    while (!vt_chkbrk(vt)) {
        prerror(vt);
        sleepms(100);
    }
    vt_printf(vt,"\n");
}
/* stream command ------------------------------------------------------------*/
static void cmd_stream(char **args, int narg, vt_t *vt)
{
    int cycle=0;
    
    trace(3,"cmd_stream:\n");
    
    if (narg>1) cycle=(int)(atof(args[1])*1000.0);
    
    while (!vt_chkbrk(vt)) {
        if (cycle>0) vt_printf(vt,ESC_CLEAR);
        prstream(vt);
        if (cycle>0) sleepms(cycle); else return;
    }
    vt_printf(vt,"\n");
}
/* ssr command ---------------------------------------------------------------*/
static void cmd_ssr(char **args, int narg, vt_t *vt)
{
    int cycle=0;
    
    trace(3,"cmd_ssr:\n");
    
    if (narg>1) cycle=(int)(atof(args[1])*1000.0);
    
    while (!vt_chkbrk(vt)) {
        if (cycle>0) vt_printf(vt,ESC_CLEAR);
        prssr(vt);
        if (cycle>0) sleepms(cycle); else return;
    }
    vt_printf(vt,"\n");
}
/* option command ------------------------------------------------------------*/
static void cmd_option(char **args, int narg, vt_t *vt)
{
    char buff[MAXSTR],*p;
    int i,n;
    
    trace(3,"cmd_option:\n");
    
    for (i=0;*rcvopts[i].name;i++) {
        if (narg>=2&&!strstr(rcvopts[i].name,args[1])) continue;
        p=buff;
        p+=sprintf(p,"%-18s =",rcvopts[i].name);
        p+=opt2str(rcvopts+i,p);
        if (*rcvopts[i].comment) {
            if ((n=(int)(buff+30-p))>0) p+=sprintf(p,"%*s",n,"");
            p+=sprintf(p," # (%s)",rcvopts[i].comment);
        }
        vt_printf(vt,"%s%s\n",modflgr[i]?"*":" ",buff);
    }
    for (i=0;*sysopts[i].name;i++) {
        if (narg>=2&&!strstr(sysopts[i].name,args[1])) continue;
        p=buff;
        p+=sprintf(p,"%-18s =",sysopts[i].name);
        p+=opt2str(sysopts+i,p);
        if (*sysopts[i].comment) {
            if ((n=(int)(buff+30-p))>0) p+=sprintf(p,"%*s",n,"");
            p+=sprintf(p," # (%s)",sysopts[i].comment);
        }
        vt_printf(vt,"%s%s\n",modflgs[i]?"*":" ",buff);
    }
}
/* set command ---------------------------------------------------------------*/
static void cmd_set(char **args, int narg, vt_t *vt)
{
    opt_t *opt;
    int *modf;
    char buff[MAXSTR];
    
    trace(3,"cmd_set:\n");
    
    if (narg<2) {
        vt_printf(vt,"specify option type\n");
        return;
    }
    if ((opt=searchopt(args[1],rcvopts))) {
        modf=modflgr+(int)(opt-rcvopts);
    }
    else if ((opt=searchopt(args[1],sysopts))) {
        modf=modflgs+(int)(opt-sysopts);
    }
    else {
        vt_printf(vt,"no option type: %s\n",args[1]);
        return;
    }
    if (narg<3) {
        vt_printf(vt,"%s",opt->name);
        if (*opt->comment) vt_printf(vt," (%s)",opt->comment);
        vt_printf(vt,": ");
        if (!vt_gets(vt,buff,sizeof(buff))||vt->brk) return;
    }
    else strcpy(buff,args[2]);
    
    chop(buff);
    if (!str2opt(opt,buff)) {
        vt_printf(vt,"invalid option value: %s %s\n",opt->name,buff);
        return;
    }
    getsysopts(&prcopt,solopt,&filopt);
    solopt[1]=solopt[0];
    
    vt_printf(vt,"option %s changed.",opt->name);
    if (strncmp(opt->name,"console",7)) {
        *modf=1;
        vt_printf(vt," restart to enable it");
    }
    vt_printf(vt,"\n");
}
/* load command --------------------------------------------------------------*/
static void cmd_load(char **args, int narg, vt_t *vt)
{
    char file[MAXSTR]="";
    
    trace(3,"cmd_load:\n");
    
    if (narg>=2) {
        strcpy(file,args[1]);
    }
    else {
        sprintf(file,"%s/%s",OPTSDIR,OPTSFILE);
    }
    resetsysopts();
    if (!loadopts(file,sysopts)) {
        vt_printf(vt,"no options file: %s\n",file);
        return;
    }
    getsysopts(&prcopt,solopt,&filopt);
    solopt[1]=solopt[0];
    
    if (!loadopts(file,rcvopts)) {
        vt_printf(vt,"no options file: %s\n",file);
        return;
    }
    vt_printf(vt,"options loaded from %s. restart to enable them\n",file);
}
// /* save command --------------------------------------------------------------*/
static void cmd_save(char **args, int narg, vt_t *vt)
{
    // char file[MAXSTR]="",comment[256],s[64];
    
    // trace(3,"cmd_save:\n");
    
    // if (narg>=2) {
    //     strcpy(file,args[1]);
    // }
    // else {
    //     sprintf(file,"%s/%s",OPTSDIR,OPTSFILE);
    // }
    // if (!confwrite(vt,file)) return;
    // time2str(utc2gpst(timeget()),s,0);
    // sprintf(comment,"%s options (%s, v.%s %s)",PRGNAME,s,VER_RTKLIB,PATCH_LEVEL);
    // setsysopts(&prcopt,solopt,&filopt);
    // if (!saveopts(file,"w",comment,rcvopts)||!saveopts(file,"a",NULL,sysopts)) {
    //     vt_printf(vt,"options save error: %s\n",file);
    //     return;
    // }
    vt_printf(vt,"save is invalid in decoder!\n");
}
// /* log command ---------------------------------------------------------------*/
static void cmd_log(char **args, int narg, vt_t *vt)
{
    trace(3,"cmd_log:\n");
    
    // if (narg<2) {
    //     vt_printf(vt,"specify log file\n");
    //     return;
    // }
    // if (!strcmp(args[1],"off")) {
    //     vt_closelog(vt);
    //     vt_printf(vt,"log off\n");
    //     return;
    // } 
    // if (!confwrite(vt,args[1])) return;
    
    // if (!vt_openlog(vt,args[1])) {
    //     vt_printf(vt,"log open error: %s\n",args[1]);
    //     return;
    // }
    vt_printf(vt,"log is invalid in decoder!\n");
}
/* help command --------------------------------------------------------------*/
static void cmd_help(char **args, int narg, vt_t *vt)
{
    char str[]="path";
    int i;
    
    if (narg<2) {
        vt_printf(vt,"%s (ver.%s %s)\n",PRGNAME,VER_RTKLIB,PATCH_LEVEL);
        for (i=0;*helptxt[i];i++) vt_printf(vt,"%s\n",helptxt[i]);
    }
    else if (strstr(str,args[1])==str) {
        for (i=0;*pathopts[i];i++) vt_printf(vt,"%s\n",pathopts[i]);
    }
    else {
        vt_printf(vt,"unknown help: %s\n",args[1]);
    }
}
/* console thread ------------------------------------------------------------*/
static void *con_thread(void *arg)
{
    const char *cmds[]={
        "start","stop","restart","solution","status","satellite","observ",
        "navidata","stream","ssr","error","option","set","load", "save","log",
        "help","?","exit","shutdown",""
    };
    con_t *con=(con_t *)arg;
    int i,j,narg;
    char buff[MAXCMD],*args[MAXARG],*p;
    
    trace(3,"console_thread:\n");
    
    vt_printf(con->vt,"\n%s** %s ver.%s %s console (h:help) **%s\n",ESC_BOLD,
              PRGNAME,VER_RTKLIB,PATCH_LEVEL,ESC_RESET);
    
    // if (!login(con->vt)) {
    //     vt_close(con->vt);
    //     con->state=0;
    //     return NULL;
    // }
 
    /* auto start if option set */
    if (start&1) { /* start with console */
        cmd_start(NULL,0,con->vt);
        start=0;
    }
    
    while (con->state) {
        
        /* output prompt */
        if (!vt_puts(con->vt,CMDPROMPT)) {
            con->state = 0;
            break;
        }
        
        /* input command */
        if (!vt_gets(con->vt,buff,sizeof(buff))) break;
        
        /* parse command */
        narg=0;
        char *r;
        for (p=strtok_r(buff," \t\n",&r);p&&narg<MAXARG;p=strtok_r(NULL," \t\n",&r)) {
            args[narg++]=p;
        }
        if (narg==0) continue;
        
        for (i=0,j=-1;*cmds[i];i++) {
            if (strstr(cmds[i],args[0])==cmds[i]) j=i;
        }
        switch (j) {
            case  0: cmd_start    (args,narg,con->vt); break;
            case  1: cmd_stop     (args,narg,con->vt); break;
            case  2: cmd_restart  (args,narg,con->vt); break;
            case  3: cmd_solution (args,narg,con->vt); break;
            case  4: cmd_status   (args,narg,con->vt); break;
            case  5: cmd_satellite(args,narg,con->vt); break;
            case  6: cmd_observ   (args,narg,con->vt); break;
            case  7: cmd_navidata (args,narg,con->vt); break;
            case  8: cmd_stream   (args,narg,con->vt); break;
            case  9: cmd_ssr      (args,narg,con->vt); break;
            case 10: cmd_error    (args,narg,con->vt); break;
            case 11: cmd_option   (args,narg,con->vt); break;
            case 12: cmd_set      (args,narg,con->vt); break;
            case 13: cmd_load     (args,narg,con->vt); break;
            case 14: cmd_save     (args,narg,con->vt); break;
            case 15: cmd_log      (args,narg,con->vt); break;
            case 16: cmd_help     (args,narg,con->vt); break;
            case 17: cmd_help     (args,narg,con->vt); break;
            case 18: /* exit */
                if (con->vt->type) con->state=0;
                break;
            case 19: /* shutdown */
                if (!strcmp(args[0],"shutdown")) {
                    vt_printf(con->vt,"rtk server shutdown ...\n");
                    sleepms(1000);
                    intflg=1;
                    con->state=0;
                }
                break;
            default:
                vt_printf(con->vt,"unknown command: %s.\n",args[0]);
                break;
        }
    }
    vt_close(con->vt);
    con->vt = NULL;
    return NULL;
}
/* open console --------------------------------------------------------------*/
static con_t *con_open(int sock, const char *dev)
{
    con_t *con;
    
    trace(3,"con_open: sock=%d dev=%s\n",sock,dev);
    
    if (!(con=(con_t *)malloc(sizeof(con_t)))) return NULL;
    
    if (!(con->vt=vt_open(sock,dev))) {
        free(con);
        return NULL;
    }
    /* start console thread */
    con->state=1;
#ifndef WIN32
    if (pthread_create(&con->thread,NULL,con_thread,con)) {
        free(con);
        return NULL;
    }
    return con;
#else
    return NULL;
#endif
}
/* close console -------------------------------------------------------------*/
static void con_close(con_t *con)
{
    trace(3,"con_close:\n");
    
    if (!con) return;
    con->state=0;
    //pthread_join(con->thread,NULL);
    free(con);
}




static void load_rtdecoder(const char *filename){
    FILE *file = fopen(filename, "r");
    char line[256];
    if (!file) {
        perror("Failed to open config file");
        return;
    }
    if (PPP_Glo.RT_flag == 0) {
        perror("Please set RT_flag = 1");
        return;
    }

    int infile_count = 0; // Counter for the infile array

    while (fgets(line, sizeof(line), file)) {
        char *start = line;
        while (isspace((unsigned char)*start)) start++;

        char *comment = strchr(start, '#');
        if (comment) {
            *comment = '\0';
        }

        char *key = strtok(start, "=");
        char *value = strtok(NULL, "\n");

        if (key && value) {
            trim(key);
            trim(value);
    /*--------------------------------- Process stream options ---------------------------------*/ 
            if (strstr(key, "strtype") != NULL) {
                int a = strlen(key);
                int index = key[strlen(key) - 2] - '0';
                if (index >= 0 && index < 8) strtype[index] = atoi(value);
            }
            else if (strcmp(key, "navmsgsel") == 0)      navmsgsel = atoi(value);
            else if (strstr(key, "strpath") != NULL) {
                int a = strlen(key);
                int index = key[strlen(key) - 2] - '0';
                if (index >= 0 && index < 8) strcpy(strpath[index], value);
            }
            else if (strstr(key, "strfmt") != NULL) {
                int a = strlen(key);
                int index = key[strlen(key) - 2] - '0';
                if (index >= 0 && index < 5) strfmt[index] = atoi(value);
            }
    /*--------------------------------- Process misc options ---------------------------------*/
            else if (strcmp(key, "svrcycle") == 0)  svrcycle = atoi(value);
            else if (strcmp(key, "timeout") == 0)   timeout = atoi(value);
            else if (strcmp(key, "reconnect") == 0) reconnect = atoi(value);
            else if (strcmp(key, "nmeacycle") == 0) nmeacycle = atoi(value);
            else if (strcmp(key, "buffsize") == 0)  buffsize = atoi(value);
    /*--------------------------------- Process station options ---------------------------------*/
            else if (strcmp(key, "sta_name") == 0)  strcpy(sta_name, value);
    /*------------------------------------------------------------------------------------------*/ 
        }
    }

    fclose(file);
}

static void print_rt_options() {
    FILE *file = fopen("rt_options_log.txt", "w");
    if (!file) {
        perror("Failed to open rt_options_log file");
        return;
    }

    fprintf(file, "System Options:\n");
    fprintf(file, "Interrupt Flag: %d\n", intflg);
    fprintf(file, "Login Password: %s\n", passwd);
    fprintf(file, "Time Type: %d (0:gpst, 1:utc, 2:jst, 3:tow)\n", timetype);
    fprintf(file, "Solution Type: %d (0:dms, 1:deg, 2:xyz, 3:enu, 4:pyl)\n", soltype);
    fprintf(file, "Solution Flag: %d (1:std+2:age/ratio/ns)\n", solflag);

    fprintf(file, "Stream Types: ");
    for (int i = 0; i < 8; i++) {
        fprintf(file, "%d ", strtype[i]);
    }
    fprintf(file, "\n");

    fprintf(file, "Stream Paths: ");
    for (int i = 0; i < 8; i++) {
        fprintf(file, "\"%s\" ", strpath[i]);
    }
    fprintf(file, "\n");

    fprintf(file, "Stream Formats: ");
    for (int i = 0; i < 5; i++) {
        fprintf(file, "%d ", strfmt[i]);
    }
    fprintf(file, "\n");

    fprintf(file, "Receiver Options: ");
    for (int i = 0; i < 3; i++) {
        fprintf(file, "\"%s\" ", rcvopt[i]);
    }
    fprintf(file, "\n");

    fprintf(file, "Server Cycle: %d ms\n", svrcycle);
    fprintf(file, "Timeout Time: %d ms\n", timeout);
    fprintf(file, "Reconnect Interval: %d ms\n", reconnect);
    fprintf(file, "NMEA Cycle: %d ms\n", nmeacycle);
    fprintf(file, "Buffer Size: %d bytes\n", buffsize);
    fprintf(file, "Navigation Message Select: %d\n", navmsgsel);
    fprintf(file, "Proxy Address: %s\n", proxyaddr);
    fprintf(file, "NMEA Request Type: %d (0:off, 1:lat/lon, 2:single)\n", nmeareq);

    fprintf(file, "NMEA Position: %.6f, %.6f, %.2f (lat, lon, height)\n", nmeapos[0], nmeapos[1], nmeapos[2]);

    fprintf(file, "Receiver Command Files: ");
    for (int i = 0; i < 3; i++) {
        fprintf(file, "\"%s\" ", rcvcmds[i]);
    }
    fprintf(file, "\n");

    fprintf(file, "Modified Flags (Receiver Options): ");
    for (int i = 0; i < 256; i++) {
        fprintf(file, "%d ", modflgr[i]);
    }
    fprintf(file, "\n");

    fprintf(file, "Modified Flags (System Options): ");
    for (int i = 0; i < 256; i++) {
        fprintf(file, "%d ", modflgs[i]);
    }
    fprintf(file, "\n");

    fprintf(file, "Monitor Port: %d\n", moniport);
    fprintf(file, "Keep Alive Flag: %d\n", keepalive);
    fprintf(file, "Auto Start: %d\n", start);
    fprintf(file, "File Swap Margin: %d s\n", fswapmargin);
    fprintf(file, "Station Name: %s\n", sta_name);

    fclose(file);
}


int main(int argc, char **argv)
{
    con_t *con[MAXCON]={0};
    int i,port=0,outstat=0,trace=22,sock=0;
    char *dev="/dev/tty",file[MAXSTR]="";
    int deamon=0;
    
    for (i=1;i<argc;i++) {
        if      (!strcmp(argv[i],"-s")) start|=1; /* console */
        else if (!strcmp(argv[i],"-nc")) start|=2; /* no console */
        else if (!strcmp(argv[i],"-o")&&i+1<argc) strcpy(file,argv[++i]);
        else if (!strcmp(argv[i], "--version")) {
            fprintf(stderr, "rtppp B2bPPP %s %s\n", VER_RTKLIB, PATCH_LEVEL);
            exit(0);
        }
        else printusage();
    }
    if (trace>0) {
        // traceopen(TRACEFILE);
        // B2b_traceopen(B2bTRACEFILE);
        // tracelevel(trace);
    }
    B2b_tracelevel(22);
    B2b_traceopen(B2bTRACEFILE);
    /* initialize rtk server and monitor port */
    rtksvrinit(&svr);
    strinit(&moni);
    
    /* load options file */
    load_config(file, &prcopt,solopt,&filopt);
    load_rtdecoder(file);

    if (start&2) { /* Start without console */
        startsvr(NULL);
    } 
     else {
        /* open device for local console */
        if (!(con[0]=con_open(0,dev))) {
            fprintf(stderr,"console open error dev=%s\n",dev);
            // if (moniport>0) closemoni();
            if (outstat>0) rtkclosestat();
            traceclose();
            B2b_traceclose();
            return -1;
        }
    }

    while (!intflg) {
        sleepms(100);
    }
    /* stop rtk server */
    stopsvr(NULL);
    
    
    traceclose();
    return 0;
}
