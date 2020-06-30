/*
 * ampd.c
 *
 */

#include "ampd.h"

#define TEST_LENGTH 1000000

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
    "of the peaks as calculated.\n"

    "Input file should only contain a float value in each line.\n\n"

    "Usage from linux command line:\n"
    " $ ampd -f [input file] -o [output file]\n"
    "Optional arguments:\n"
    "\t-v : verbose"
    "\t-h : print help"
    "\t-a : output intermediary data"
            );
}
void printf_data(float *data, int n){

    for(int i=0; i<n; i++)
        printf("%f\n",data[i]);
}

void fprintf_data(FILE *fp, float *data, int n){

    for(int i=0; i<n; i++){
        fprintf(fp, "%.5f\n",data[i]);
    }
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
#define FBUF 32
int fetch_data(char *path, float *data, int n, int ind){

    FILE *fp;
    int i, count;
    char buf[FBUF];
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    fp = fopen(path, "r");
    if(fp == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    count = 0; i = 0;
    while(fgets(buf, FBUF, fp)){
        if(count < ind){
            count++;
            continue;
        } else if(count > ind +n ){
            break;
        } else {
            sscanf(buf, "%f\n",data + i);
            count++;
            i++;
        }
    }
    // This is viable, but slow
    /*
    count = 0;
    while((read = getline(&line, &len, fp)) != -1){
        if(count > ind +n )
            break;
        if(count < ind){
            count++;
            continue;
        } else{
            line[line[strlen(line)-1]] = '\0';
            data[i] = atof(line);
        }
    }
    */
    fclose(fp);
    return 0;
}

/**
 * Least-squares linear regression. Data is assumed to be uniformly spaced.
 *
 */
int linear_fit(float *data, int n, double ts, double *a, double *b, double *r){

    double x; // time
    double sumx = 0.0;
    double sumx2 = 0.0;
    double sumxy = 0.0;
    double sumy = 0.0;
    double sumy2 = 0.0;
    double denom;

    for(int i=0; i<n;i++){
        x = (double)i * ts;
        sumx += x;
        sumx2 += sqrt(x);
        sumxy += x * data[i];
        sumy += data[i];
        sumy2 += sqrt(data[i]);
    }
    denom = (n * sumx2 - sqrt(sumx));
    if(denom == 0){
        //cannot solve
        *a = 0;*b = 0;*r = 0;
        return -1;
    }
    *a = (n * sumxy  -  sumx * sumy) / denom;
    *b = (sumy * sumx2  -  sumx * sumxy) / denom;
    if (r!=NULL) {
        *r = (sumxy - sumx * sumy / n) /    /* compute correlation coeff */
              sqrt((sumx2 - sqrt(sumx)/n) *
              (sumy2 - sqrt(sumy)/n));
    }
    return 0;
}

/**
 * Subtract linear trendline from data.
 *
 */
int linear_detrend(float *data, int n, double ts, double a, double b){

    float t = 1.0; // unit time interval
    for(int i=0; i<n; i++){
        data[i] -= a * t * i * ts + b;
    }
    return 0;
}

int main(int argc, char **argv){

    int opt;
    int verbose = VERBOSE_DEFAULT;
    int output_all = OUTPUT_ALL_DEFAULT;
    char infile[MAX_PATH_LEN] = {0};
    char outfile[MAX_PATH_LEN] = {0};
    FILE *fp_out;
    // aux output files
    char detrendf[MAX_PATH_LEN] = {0};
    FILE *fp_detrendf;

    float *data;    // timeseries
    int n;          // number of elements in timeseries, in a batch, dynamic
    double ts = TIMESTEP_DEFAULT;
    int datalen;    // full data length
    int cycles;     // number of data batches
    int ind;        // index of next batch
    int i, j;
    // linear fitting data = at + b
    double a = 0; double b = 0; double r = 0;
    struct Mtx lms;

    // aux output files
    char aux_dir[MAX_PATH_LEN] = {0};
    char outfile_def[] = "ampd.out";
    char detrend_def[] = "ampd.out.detrend";
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
            case 't':
                ts = atoi(optarg);
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
    // settings aux files
    if(output_all == 1){
        snprintf(detrendf, sizeof(detrendf), "%s/%s",aux_dir, detrend_def);

        fp_detrendf = fopen(detrendf, "w+");
        if(fp_detrendf == NULL){
            perror("fopen");
            exit(EXIT_FAILURE);
        }
    }
    // count elements
    datalen = count_char(infile, '\n');
    cycles = (int) (ceil(datalen / (double) DATA_BUF));

    if(verbose == 1){
        printf("infile: %s\n", infile);
        printf("outfile: %s\n", outfile);
        printf("aux_dir: %s\n", aux_dir);
        printf("datalen: %d\n", datalen);
        printf("cycles: %d\n", cycles);

        printf("\nProgess:  (| = 10 cycles)\n");

    }
    fp_out = fopen(outfile, "w+");
    if(fp_out == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    for(i=0; i<cycles; i++){

        ind = i * (int)DATA_BUF;
        if(i == cycles-2)
            n = datalen - i * (int)DATA_BUF;
        else
            n = (int)DATA_BUF;
        data = malloc(sizeof(float) * n);
        memset(data, 0, sizeof(data));
        fetch_data(infile, data, n, ind);
        //printf_data(data, n);
        linear_fit(data, n, ts, &a, &b, &r);
        linear_detrend(data, n, ts, a, b);
        if(r != r)
            continue;
        //printf("a: %lf, b: %lf, r: %lf\n", a, b, r);
        // save aux
        if(output_all == 1){
            fprintf_data(fp_detrendf, data, n);
        }
        if(verbose == 1){
            if(i % 10 == 0){
                printf("|");
                fflush(stdout);
            }
        }
    }

    if(verbose == 1)
        printf("\n");
    fclose(fp_out);
    fclose(fp_detrendf);
    
    return 0;
}
