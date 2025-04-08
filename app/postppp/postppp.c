/*------------------------------------------------------------------------------
* rnx2rtkp.c : read rinex obs/nav files and compute receiver positions
*
*          Copyright (C) 2007-2016 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:55:16 $
* history : 2007/01/16  1.0 new
*           2007/03/15  1.1 add library mode
*           2007/05/08  1.2 separate from postpos.c
*           2009/01/20  1.3 support rtklib 2.2.0 api
*           2009/12/12  1.4 support glonass
*                           add option -h, -a, -l, -x
*           2010/01/28  1.5 add option -k
*           2010/08/12  1.6 add option -y implementation (2.4.0_p1)
*           2014/01/27  1.7 fix bug on default output time format
*           2015/05/15  1.8 -r or -l options for fixed or ppp-fixed mode
*           2015/06/12  1.9 output patch level in header
*           2016/09/07  1.10 add option -sys
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include "rtklib.h"

#define PROGNAME    "rnx2rtkp"          /* program name */
#define MAXFILE     16                  /* max number of input files */
#define MAXSTR      1024                /* max length of a stream */

PPPGlobal_t PPP_Glo = {0};
char *infile[MAXFILE] = {0};
char outfile[MAXSTR] = {0};   //{ "abpo0320.pos" };
int n_file;


/* help text -----------------------------------------------------------------*/
static const char *help[]={
"",
" usage: rnx2rtkp [option]... file file [...]",
"",
" Read RINEX OBS/NAV/GNAV/HNAV/CLK, SP3, SBAS message log files and ccompute ",
" receiver (rover) positions and output position solutions.",
" The first RINEX OBS file shall contain receiver (rover) observations. For the",
" relative mode, the second RINEX OBS file shall contain reference",
" (base station) receiver observations. At least one RINEX NAV/GNAV/HNAV",
" file shall be included in input files. To use SP3 precise ephemeris, specify",
" the path in the files. The extension of the SP3 file shall be .sp3 or .eph.",
" All of the input file paths can include wild-cards (*). To avoid command",
" line deployment of wild-cards, use \"...\" for paths with wild-cards.",
" Command line options are as follows ([]:default). A maximum number of", 
" input files is currently set to 16. With -k option, the",
" processing options are input from the configuration file. In this case,",
" command line options precede options in the configuration file.",
"",
" -?        print help",
" -k file   input options from configuration file [off]",
" -o file   set output file [stdout]",
" -ts ds ts start day/time (ds=y/m/d ts=h:m:s) [obs start time]",
" -te de te end day/time   (de=y/m/d te=h:m:s) [obs end time]",
" -ti tint  time interval (sec) [all]",
" -p mode   mode (0:single,1:dgps,2:kinematic,3:static,4:static-start,",
"                 5:moving-base,6:fixed,7:ppp-kinematic,8:ppp-static,9:ppp-fixed) [2]",
" -m mask   elevation mask angle (deg) [15]",
" -sys s[,s...] nav system(s) (s=G:GPS,R:GLO,E:GAL,J:QZS,C:BDS,I:IRN) [G|R]",
" -f freq   number of frequencies for relative mode (1:L1,2:L1+L2,3:L1+L2+L5) [2]",
" -v thres  validation threshold for integer ambiguity (0.0:no AR) [3.0]",
" -b        backward solutions [off]",
" -c        forward/backward combined solutions [off]",
" -i        instantaneous integer ambiguity resolution [off]",
" -h        fix and hold for integer ambiguity resolution [off]",
" -bl bl,std baseline distance and stdev",
" -e        output x/y/z-ecef position [latitude/longitude/height]",
" -a        output e/n/u-baseline [latitude/longitude/height]",
" -n        output NMEA-0183 GGA sentence [off]",
" -g        output latitude/longitude in the form of ddd mm ss.ss' [ddd.ddd]",
" -t        output time in the form of yyyy/mm/dd hh:mm:ss.ss [sssss.ss]",
" -u        output time in utc [gpst]",
" -d col    number of decimals in time [3]",
" -s sep    field separator [' ']",
" -r x y z  reference (base) receiver ecef pos (m) [average of single pos]",
"           rover receiver ecef pos (m) for fixed or ppp-fixed mode",
" -l lat lon hgt reference (base) receiver latitude/longitude/height (deg/m)",
"           rover latitude/longitude/height for fixed or ppp-fixed mode",
" -y level  output soltion status (0:off,1:states,2:residuals) [0]",
" -x level  debug trace level (0:off) [0]"
};
/* show message --------------------------------------------------------------*/
extern int showmsg(const char *format, ...)
{
    va_list arg;
    va_start(arg,format); vfprintf(stderr,format,arg); va_end(arg);
    fprintf(stderr,"\r");
    return 0;
}
// extern void settspan(gtime_t ts, gtime_t te) {}
// extern void settime(gtime_t time) {}

/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    int i;
    for (i=0;i<(int)(sizeof(help)/sizeof(*help));i++) fprintf(stderr,"%s\n",help[i]);
    exit(0);
}

static void load_post(const char *filename){
    FILE *file = fopen(filename, "r");
    char line[256];
    if (!file) {
        perror("Failed to open config file");
        return;
    }
    if (PPP_Glo.RT_flag == 1) {
        perror("Please set RT_flag = 0");
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
    /*--------------------------------- Process infile options ---------------------------------*/ 
            if (strstr(key, "infile") != NULL && infile_count < MAXFILE) {
                infile[infile_count] = (char *)malloc(strlen(value) + 1);
                if (infile[infile_count]) {
                    strcpy(infile[infile_count], value);
                    infile_count++;
                }
            }
            else if (strstr(key, "outfile") != NULL) {
                strcpy(outfile, value);
            }
    /*------------------------------------------------------------------------------------------*/ 
        }
    }

    n_file = infile_count;  // Set the number of input files

    fclose(file);
}

static void extract_station_and_doy(const char *filename, char *result) {
    // 获取文件名（去掉路径部分）
    const char *base_name = strrchr(filename, '/');
    if (base_name == NULL) {
        base_name = filename;
    } else {
        base_name++;
    }

    // 找到最后一个点号的位置
    const char *dot_pos = strrchr(base_name, '.');

    // 提取年份并转换为四位数
    char two_digit_year[3];
    strncpy(two_digit_year, dot_pos + 1, 2);
    two_digit_year[2] = '\0';

    char four_digit_year[5];
    sprintf(four_digit_year, "20%s", two_digit_year);

    // 提取站点名（前4个字符）
    char station_name[5];
    strncpy(station_name, base_name, 4);
    station_name[4] = '\0';

    // 提取DOY（第5到7个字符）
    char doy[4];
    strncpy(doy, base_name + 4, 3);
    doy[3] = '\0';

    // 组合站点名、年份和DOY
    sprintf(result, "%s%s%s", station_name, four_digit_year, doy);
}


/* rnx2rtkp main -------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char *dev="",file[MAXSTR]="";
    gtime_t ts={0},te={0};
    double tint=0.0,es[]={2000,1,1,0,0,0},ee[]={2000,12,31,23,59,59},pos[3];
    int i,j,ret,n;
    char *p;

    for (int i=1;i<argc;i++) {

        if (!strcmp(argv[i],"-o")&&i+1<argc) strcpy(file,argv[++i]);
        // else printhelp();
    }

    prcopt_t prcopt=prcopt_default;
    solopt_t solopt=solopt_default;
    filopt_t filopt={""};
    
// #ifdef _DEBUG
//     printf("Start with DEBUG MODE! \n");
//     load_config(file, &prcopt, &solopt, &filopt);
//     load_post(file);
// #else
    // printf("Start with RELEASE MODE! \n");
    // load_config(file, &prcopt, &solopt, &filopt, infile, outfile, &n);
    load_config(file, &prcopt, &solopt, &filopt);
    load_post(file);
    PPP_Glo.b2b_flag=0;
    if(prcopt.sateph==5) {
        PPP_Glo.b2b_flag=1;
    }
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-in") && i + 1 < argc) {
            i++;
            while (i < argc && argv[i][0] != '-') {
                infile[n++] = argv[i++];
            }
            i--;
        } 
        else if (!strcmp(argv[i], "-out") && i + 1 < argc) {
            strcpy(outfile, argv[++i]);
        }else if (!strcmp(argv[i], "-s") && i + 1 < argc) {
            prcopt.navsys = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-m") && i + 1 < argc) {
            prcopt.mode = atoi(argv[++i]);
        }
    }
// #endif
    // strncpy(infile[0],prcopt.sationname,4);  
    char result[20] = {0};
    extract_station_and_doy(infile[0], result);
    strcpy(prcopt.sationname,result);

    ret = postpos(ts, te, tint, 0.0, &prcopt, &solopt, &filopt, (const char **)infile, n_file, outfile, "", "");
    
    if (!ret) fprintf(stderr,"%40s\r","");
    return ret;
}
