/*
 * ampd.c
 *
 */

#include "ampd.h"

#define TEST_LENGTH 10000

/**
 * Print description and general usage
 */
void printf_help(){

    printf(
    "AMPD\n"
    "================================================================\n"
    "Peak detection algorithm for quasiperiodic data. Main usage of this "
    "implementation is detection of peaks in rat physiological data: "
    "respiration and pulsoxymmetry waveforms.\n\n"

    "Reference paper:\n"
    "An Efficient Algorithm for Automatic Peak Detection in Noisy "
    "Periodic and Quasi-Periodic Signals\n"
    "DOI:10.3390/a5040588\n\n"

    "This program takes a file input which contains a single timeseries "
    "of quasi-periodic data. The output is a file containing the indices "
    "of the peaks as calculated.\n\n"

    "Usage from linux command line:\n"
    " $ ampd -f [input file] -o [output file]\n"
    "Optional arguments:\n"
    "\t-v : verbose"
    "\t-h : print help"
    "\t-a : output intermediary data"
            );
}
/**
 * Least-squares linear regression.
 *
 */
int linreg(float *data, int n, double *a, double *b){

    return 0;
}

/**
 * Subtract linear trendline from data.
 *
 */
int linear_detrend(float *data, int n, double a, double b){

    return 0;
}

int main(int argc, char **argv){

    int opt;
    int verbose = VERBOSE_DEFAULT;
    int output_all = OUTPUT_ALL_DEFAULT;
    char infile[MAX_PATH_LEN];
    char outfile[MAX_PATH_LEN];

    float *data;
    int datalen = TEST_LENGTH; // load only partial data for dev
    int i, j;

    struct Mtx lms;

    // parse options
    while((opt = getopt(argc, argv, "hvf:o:")) != -1){
        switch(opt){
            case 'h':
                printf_help();
                return 0;
            case 'v':
                verbose = 1;
                break;
            case 'o':
                strcpy(outfile, optarg);
                break;
            case 'f':
                strcpy(infile, optarg);
                break;
            case 'a':
                output_all = 1;
                break;

        }
    }
}
