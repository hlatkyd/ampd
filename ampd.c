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
 * Count occurrences of a character in a file
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
 * Load a part of the full timeseries data into memory from file.
 * File should only contain one float value on each line.
 */
int fetch_data(char *path, float *data, int n, int ind){

    FILE *fp;
    int i;
    size_t len = 0;
    ssize_t read;
    char *tok;
    char *line = NULL;

    fp = fopen(path, "r");
    if(fp == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    for(i=0; i<n; i++){

    }
    fclose(fp);
}

/**
 * Least-squares linear regression. Data is assumed to be uniformly spaced.
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
    char outfile[MAX_PATH_LEN] = {0};

    float *data;
    int datalen = TEST_LENGTH; // load only partial data for dev
    int i, j;

    struct Mtx lms;

    // aux output files
    char aux_dir[MAX_PATH_LEN] = {0};
    char outfile_def[] = "ampd.out";
    // parse options
    while((opt = getopt(argc, argv, "hvf:o:a")) != -1){
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
    // setting to defaults if arguments were not given
    if(strcmp(aux_dir,"")==0){
        getcwd(aux_dir, sizeof(aux_dir));
    }
    if(strcmp(outfile,"")==0){
        snprintf(outfile, sizeof(outfile), "%s/%s",aux_dir, outfile_def);
    }

    if(output_all == 1){
        ;
    }
    if(verbose == 1){
        printf("infile: %s\n", infile);
        printf("outfile: %s\n", outfile);
        printf("aux_dir: %s\n", aux_dir);
    }
    return 0;
}
