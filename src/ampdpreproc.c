#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <unistd.h>
#include "filters.h"

#define V_MIN 9
#define V_MAJ 0
#define MAX_PATH 1024

#define DEF_OUTPUT_PATH "ampdprerpoc.out"
#define DEF_SMOOTH_WINDOW 2
#define DEF_SAMPLING_RATE 200

static int verbose_flag;
static int smooth_only_flag;

void printf_version(){

    printf("ampdpreproc, v%d.%d\n",V_MAJ, V_MIN);
}

void printf_help(){
    
    printf("\nampdpreproc\n"
           "===========\n"
           "Prepare a single data series for peak counting with the AMPD \n"
           "routine. This utility program applies simple time-domeain filters\n"
           "and smoothing to uniformly sammpled data ina  file input.\n"
           "The input file should contain one data value each line.\n"
           "The output is a new file of the same format.\n\n"
           "Usage:\n"
           "$ ampdpreproc --[args]=[vals]\n"
           "-f --infile=[INFILE]\n"
           "-o --outfile=[OUTFILE]\n"
           "-s --sampling-rate=[SAMPLING_RATE]\n"
           "-h --highpass-cutoff=[HIGPASS_CUTOFF]\n"
           "-l --lowpass-cutoff=[LOWPASS_CUTOFF]\n"
           "--verbose\n"
           "--help\n"
           );
}
/* Count occurrences of a character ina  file. Useful for counting lines*/
int count_char(char *path, char cc);
int get_fprec_from_str(char *str);

int main(int argc, char **argv){

    int opt;
    int index;
    FILE *verbose_fs = stdout; // output verbose stuff here

    double sampling_rate = DEF_SAMPLING_RATE;  // sampling rate in Hz
    double highpass_cutoff = 0.0;     // cutoff frequency in Hz for high pass filter
    double lowpass_cutoff = 0.0;    // cutoff frequency in Hz for lowpass filter
    float *data;                // input date, later processed as well
    int n;                      // number of data points
    int fprec;                  // precision of data

    char infile[MAX_PATH] = {0};
    char outfile[MAX_PATH] = DEF_OUTPUT_PATH;
    FILE *fp;

    // moving avg
    int w = DEF_SMOOTH_WINDOW;

    if(argc == 1){
        printf_help();
        return 0;
    }

    while(1){
        /*
         * Option parsing
         */
        static struct option long_options[] = {
            {"help",no_argument, 0, 0},
            {"verbose",no_argument, &verbose_flag, 'v'},
            {"version",no_argument, 0, 1},
            {"infile",required_argument, 0, 'f'},
            {"outfile",required_argument, 0, 'o'},
            {"sampling-rate",required_argument, 0, 's'},
            {"lowpass-freq", required_argument, 0, 'l'},
            {"highpass-freq", required_argument, 0, 'h'},
            {"smooth-only",no_argument, &smooth_only_flag, 2},
            {0,0,0,0}
        };

        opt = getopt_long(argc, argv, "vh:l:s:t:f:o:",long_options, &index);

        if(opt == -1)
            break;

        switch (opt){

            case 0:
                printf_help();
                return 0;
            case 1:
                printf_version();
                break;
            case 's':
                sampling_rate = atof(optarg);
                break;
            case 'l':
                lowpass_cutoff = atof(optarg);
                break;
            case 'h':
                highpass_cutoff = atof(optarg);
                break;
            case 'v':
                verbose_flag = 1;
                break;
            case 'f':
                strcpy(infile, optarg);
                break;
            case 'o':
                strcpy(outfile, optarg);
                break;
            case 2:
                smooth_only_flag = 1;
                break;

        }
    }
    // input check
    if(strcmp(infile,"")==0){
        fprintf(stderr, "No input file given, exiting...\n");
        exit(1);
    }
    if(verbose_flag == 1){
        fprintf(verbose_fs, "verbose=%d\n",verbose_flag);
        fprintf(verbose_fs, "sampling_rate=%lf\n",sampling_rate);
        fprintf(verbose_fs, "lowpass_cutoff=%lf\n",lowpass_cutoff);
        fprintf(verbose_fs, "highpass_cutoff=%lf\n",highpass_cutoff);
        fprintf(verbose_fs, "smooth_only=%d\n",smooth_only_flag);
        fprintf(verbose_fs, "infile=%s\n",infile);
        fprintf(verbose_fs, "outfile=%s\n",outfile);
    }
    /* count data length*/
    n = count_char(infile, '\n');
    data = malloc(sizeof(float) * n);
    /* reading input data */
    fp = fopen(infile, "r");
    if(fp == NULL){
        fprintf(stderr, "Cannot open file on path '%s'\n",infile);
        exit(1);
    }
    char buf[32];
    int i = 0;
    while(fgets(buf, 32, fp)){
        if(i==0){
            // check precision, and use this later for output
            fprec = get_fprec_from_str(buf);
        }
        sscanf(buf, "%f\n",&data[i]);
        i++;
    }
    fclose(fp);
    /* Apply smoothing*/
    movingavg(data, n, w);
    /* Apply filters*/
    if(highpass_cutoff != 0.0)
        tdhpfilt(data, n, sampling_rate, highpass_cutoff);
    if(lowpass_cutoff != 0.0)
        tdlpfilt(data, n, sampling_rate, lowpass_cutoff);

    /* Save filtered, smoothed data*/

    fp = fopen(outfile, "w+");
    if(fp == NULL){
        fprintf(stderr, "cannot open file for writing '%s'\n",outfile);
        exit(1);
    }
    for(i=0; i<n; i++){
        fprintf(fp, "%.*f\n",fprec, data[i]);
    }
    fclose(fp);
    free(data);
    return 0;
}

/**
 * Count occurrences of a character in a file.
 *
 */
int count_char(char *path, char cc){
    int count = 0;
    FILE *fp;
    char c;
    fp = fopen(path, "r");
    if(fp == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    for (c = getc(fp); c != EOF; c = getc(fp))
        if(c == cc)
            count++;
    fclose(fp);
    return count;
}
/**
 * count numbers after decimal separator
 * example input: "74.45673", output=5
 */
int get_fprec_from_str(char *str){

    int prec;
    int i, n, count = 0;
    for(i=0; i<strlen(str); i++){
        if(str[i] == '.'){
            n = i;
            break;
        }
    }
    prec = strlen(str)-(n+1);
    return prec;
}

