/*------------------------------------------------------------------------------
* postdecoder.c : test decoding for various receiver-specific custom data formats
*
* Built by Liu@APM(SDUST), Jiang@SEU(SDUST), Li@SEU(SDUST), Zhu@SEU(SDUST) for RTCM3 decoder
*
* version : $Revision:$ $Date:$
* history : 14 Sep 2023  1.0 created file for initial decoding tests by Liu@APM, Jiang@SEU, Li@SEU, Zhu@SEU
*           05 Oct 2023  1.1 added RTCM3 SSR format decoding support by Liu@APM
*           19 Nov 2024  1.2 added Sinan/Unicore B2b custom format decoding by Liu@APM,Xu@CUGB
*
* description : This program tests decoding for receiver manufacturers' custom data formats.
*               Input design simulates the RTKLIB framework for enhanced compatibility.
*-----------------------------------------------------------------------------*/
#include "postdecoder.h"
#include "rtklib.h"

#define TRACEFILE   "./B2bPPP_%Y_%m_%d.trace"    /* debug trace file */
#define B2bTRACEFILE "./B2bPPP_%Y_%m_%d.B2bssr"   /* B2b SSR debug file */

// static nav_t navs={0};          /* navigation data structure */


PPPGlobal_t PPP_Glo = {0};        /* PPP global parameters */


int main(int argc, char **argv)
{
    int i;
    char output_path[1024];        // Buffer for output file path
    const char *B2bfilepath = NULL; // Input file path
    const char *format_str = NULL;  // Format type identifier
    const char *output_specified = NULL; // User-specified output path

    // Parse command-line arguments
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-U") == 0) {
            format_str = "U";
        } else if (strcmp(argv[i], "-S") == 0) {
            format_str = "S";
        } else if (strcmp(argv[i], "-in") == 0 && i + 1 < argc) {
            B2bfilepath = argv[++i];
        } else if (strcmp(argv[i], "-out") == 0 && i + 1 < argc) {
            output_specified = argv[++i];
        } else {
            printf("Invalid option: %s\n", argv[i]);
            printf("Usage: postdecoder [-U|-S] -in B2b_filepath [-out B2bSSR_resultpath]\n");
            return -1;
        }
    }

    // Validate required parameters
    if (format_str == NULL || B2bfilepath == NULL) {
        printf("Missing required options. Usage: postdecoder [-U|-S] -in B2b_filepath [-out B2bSSR_resultpath]\n");
        return -1;
    }

    // Determine data format
    int format;
    if (format_str[0] == 'S' || format_str[0] == 's') {
        format = STRFMT_SINO;
    } else if (format_str[0] == 'U' || format_str[0] == 'u') {
        format = STRFMT_UNICORE;
    } else {
        printf("Invalid format type. Use '-U' for UNICORE or '-S' for SINO.\n");
        return -1;
    }

    // Configure output path
    if (output_specified != NULL) {
        strcpy(output_path, output_specified); // Use user-defined path
    } else {
        strcpy(output_path, B2bfilepath);      // Append default suffix
        strcat(output_path, ".B2bSSR");
    }

    raw_t raw0 = {0};            // Raw data container
    nav_t navs = {0};            // Navigation data container

    if (init_raw(&raw0, format) != 1) {
        printf("Error initializing raw data handler!\n");
        return -1;
    }

    FILE *B2bfp = fopen(B2bfilepath, "rb");
    if (B2bfp == NULL) {
        printf("Error: Failed to open: %s\n", B2bfilepath);
        return -1;
    }

    // tracelevel(100);             // Enable full tracing
    B2b_tracelevel(22);
    B2b_traceopen(output_path);  // Initialize tracing with output path

    printf("Decoding started...\n");
    while (1) {
        if (input_rawf(&raw0, format, B2bfp) < -1) break;
        for (i = 0; i < MAXSAT; i++) {
            if (!raw0.nav.B2bssr[i].update) continue;
            if (raw0.num_PPPB2BINF02 != 0) {
                raw0.nav.B2bssr[i].udi[0] = timediff(raw0.nav.B2bssr[i].t0[0], navs.B2bssr[i].t0[0]);
                if (raw0.nav.B2bssr[i].udi[0] > 86400) raw0.nav.B2bssr[i].udi[0] = 0;
            }
            if (raw0.num_PPPB2BINF03 != 0) {
                raw0.nav.B2bssr[i].udi[1] = timediff(raw0.nav.B2bssr[i].t0[1], navs.B2bssr[i].t0[1]);
                if (raw0.nav.B2bssr[i].udi[1] > 86400) raw0.nav.B2bssr[i].udi[1] = 0;
            }
            if (raw0.num_PPPB2BINF04 != 0) {
                raw0.nav.B2bssr[i].udi[2] = timediff(raw0.nav.B2bssr[i].t0[2], navs.B2bssr[i].t0[2]);
                if (raw0.nav.B2bssr[i].udi[2] > 86400) raw0.nav.B2bssr[i].udi[2] = 0;
            }
            navs.B2bssr[i] = raw0.nav.B2bssr[i];
            raw0.nav.B2bssr[i].update = 0;
        }
        raw0.num_PPPB2BINF01 = 0;
        raw0.num_PPPB2BINF02 = 0;
        raw0.num_PPPB2BINF03 = 0;
        raw0.num_PPPB2BINF04 = 0;
    }

    traceclose();                // Close trace files
    B2b_traceclose();

    printf("Decoding completed!\n");

    return 0;
}