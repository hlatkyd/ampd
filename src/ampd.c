/*
 * ampd.c
 *
 */

#include "ampd.h"

#define TEST_LENGTH 1000000 // TODO wtf is this?

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

    "Input file should only contain a float value in each line.\n"
    "Main output file contains the indices of peaks, while aux output "
    "directory contains various intermediate data for error checking\n\n"

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
    return;
}

void fprintf_data(FILE *fp, float *data, int n){

    for(int i=0; i<n; i++){
        fprintf(fp, "%.5f\n",data[i]);
    }
    return;
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
 * Extract basename from a path, and omit extension as well.
 * Return 0 on success, -1 on error. Extracted string is put
 * into bname.
 */
void extract_raw_filename(char *path, char *bname, int buf){

    char *tmp;
    memset(bname, 0, buf);
    tmp = basename(path);
    if(strlen(tmp) > buf){
        fprintf(stderr, "extract_raw_filename: insufficent buffer\n");
        exit(EXIT_FAILURE);
    }
    for(int i=0; i<strlen(tmp); i++){
        if(tmp[i] != '.')
            bname[i] = tmp[i];
        else
            break;
    }
    return;
}

/**
 * Save list into file, one value per line.
 */
void save_data(float *data, int n, char *path){

    FILE *fp;
    fp = fopen(path, "w");
    printf("path=%s\n",path);
    if(fp == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    for(int i=0; i<n; i++){
        fprintf(fp, "%.3f\n",data[i]);
    }
    printf("HERE\n");
    fclose(fp);
    return;
}
/**
 * Save matrix to a tab delimited file, withoud any headers.
 */
void save_mtx(struct Mtx *mtx, char *path){

    FILE *fp;
    int i, j;
    fp = fopen(path, "w");
    if(fp == NULL){
        perror("fopen, exiting...");
        exit(EXIT_FAILURE);
    }
    printf("rows, cols: %d, %d\n", mtx->rows, mtx->cols);
    for(i=0; i<mtx->rows; i++){
        for(j=0; j<mtx->cols; j++){
            if(j == mtx->cols-1)
                fprintf(fp,"%.3f\n",mtx->data[i][j]);
            else
                fprintf(fp,"%.3f\t",mtx->data[i][j]);
        }
    }
    fclose(fp);
    return;
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
    data = malloc(sizeof(float) * n);
    memset(data, 0, sizeof(data));
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

    for(int i=0; i<n; i++){
        data[i] -= a * i * ts + b;
    }
    return 0;
}

/**
 * Allocate memory for local maxima scalogram matrix
 *
 */
struct Mtx* malloc_mtx(int rows, int cols){

    struct Mtx *mtx = malloc(sizeof(struct Mtx));
    mtx->rows = rows;
    mtx->cols = cols;
    mtx->data = malloc((mtx->rows*sizeof(float *)));
    for(int i=0; i<mtx->rows;i++){
        mtx->data[i] = malloc(mtx->cols * sizeof(float));
    }
    return mtx;
}
/**
 * Caclulate local maxima scalogram.
 * Moving window w_k = {2k | k=1,2...L} where L = ceil(n/2)
 * For the timeseries x_i the matrix elements:
 *  m_k_i = 0 if (x_i-1 > x_i-k-1 && x_i-1 > x_i+1-1)
 *  m_k_i = r + alpha othervise
 * where r is double rand[0,1]
 *
 */
void calc_lms(struct Mtx *lms, float *data){

    int i, k;
    int n, l;
    float rnd;
    float a = ALPHA;
    n = lms->cols;
    l = lms->rows;
    printf("cols=%d, rows=%d\n",n, l);
    for(i = 0; i<n; i++){
        for(k=0; k<l; k++){
            rnd = (float) rand() / (float)RAND_MAX;
            //printf("i=%d, k=%d\n",i,k);
            if(i < k){
                lms->data[k][i] = rnd + a;
                continue;
            }
            if(i>n-k+1){
                lms->data[k][i] = rnd + a;
                continue;
            }
            if(data[i-1] > data[i-k-1] && data[i-1] > data[i+k-1])
                lms->data[k][i] = 0.0;
            else
                lms->data[k][i] = rnd + a;
        }
    }
    return;
}
/**
 * Calculate the vector gamma, by summing the LMS row-wise.
 *
 */
void row_sum_lms(struct Mtx *lms, double *gamma){
    ;
}

int concat_peaks(int *sum_peaks, int *peaks, int ind){

    return 0;
}
/**
 * Main routine for peak detection on a dataseries
 *
 * @param data      Input dataseries, previosly loaded from file
 * @param n         Length of dataseries
 * @param lms       Local maxima scalogram matrix
 * @param rlm       Rscaled local maxima scalogram
 * @param gamma     Vector coming from the row-wise summation of the LMS matrix.
 * @param sigma     Vector from the column-wise standard deviation of the rescaled
 *                  LMS matrix
 * @param peaks     Vector containing the indices of peaks corresponding to the 
 *                  original input dataseries
 *
 * @return          Number of peaks if successful, -1 on error.
 */
int ampd(float *data, int n, struct Mtx *lms, struct Mtx *rlms, 
         double *gamma, double *sigma, int *peaks){

    double a = 0; double b = 0; double r = 0;
    double ts = TIMESTEP_DEFAULT;
    int l;
    int n_peaks;
    // LMS matrix is N x L, so make L:
    l = (int)ceil(n / 2) - 1;
    linear_fit(data, n, ts, &a, &b, &r);
    if(r != r) // return if fitting was not working
        fprintf(stderr, "ampd: linear fit error r=%lf\n",r);
        return -1;
    linear_detrend(data, n, ts, a, b);
    lms = malloc_mtx(l, n);
    calc_lms(lms, data);

    return n_peaks;

}
/**
 * Main handles the command line input parsing and output file management.
 */
int main(int argc, char **argv){

    int opt;
    int verbose = VERBOSE;
    int output_all = OUTPUT_ALL;

    char infile[MAX_PATH_LEN] = {0};
    char infile_basename[MAX_PATH_LEN]; // input, without dir and extension
    char outfile[MAX_PATH_LEN] = {0}; // main output file with indices of peaks
    char cwd[MAX_PATH_LEN]; // current directory
    FILE *fp_out;
    // aux output files
    char detrend_path[MAX_PATH_LEN];
    char lms_path[MAX_PATH_LEN];
    char rlms_path[MAX_PATH_LEN];
    char gamma_path[MAX_PATH_LEN];
    char sigma_path[MAX_PATH_LEN];
    char peaks_path[MAX_PATH_LEN];

    // batch processing
    float *data;        // timeseries
    double ts;          // timestep
    int n;              // number of elements in timeseries, in a batch, dynamic
    int ind;            // index of next batch
    int datalen;        // full data length
    int cycles;         // number of data batches
    int n_peaks;        // main output, number of peaks
    int sum_n_peaks;    // summed peak number from all batches
    int *sum_peaks;      // concatenated vector from peak indices
    // ampd routine pointers
    struct Mtx *lms;
    struct Mtx *rlms;
    double *gamma;
    double *sigma;
    int *peaks;

    // main output file, containing only the peak indices
    char outfile_def[] = "ampd_out_peaks"; //

    // aux output files
    char aux_dir[MAX_PATH_LEN] = {0};
    char aux_dir_def[] = "ampd.out"; // full default is cwd plus this

    char detrend_def[] = "ampd.out.detrend";
    char lms_def[] = "ampd.out.lms";
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
    /* Setting up output paths.
     * Main output path is given as command line argument, if not it is the cwd
     * General format of aux output files:
     * aux_root/data_ID/batch_[i]/out_intermediate_files
     * aux_root is the current working dir by default, but can be set with a
     * command line argument. data_ID contains the name of the file where the
     * input data is coming from.
     */
    getcwd(cwd, sizeof(cwd));
    extract_raw_filename(infile, infile_basename, sizeof(infile_basename));
    // setting to defaults if arguments were not given
    if(strcmp(aux_dir,"")==0){
        snprintf(aux_dir, sizeof(aux_dir), "%s/%s",cwd, aux_dir_def);
    }
    if(strcmp(outfile,"")==0){
        snprintf(outfile, sizeof(outfile), "%s/%s",aux_dir, outfile_def);
    }
    // settings aux files
    if(output_all == 1){
        snprintf(detrend_path,sizeof(detrend_path),"%s/%s",aux_dir,detrend_def);
        snprintf(lms_path,sizeof(lms_path),"%s/%s",aux_dir, lms_def);
    }

    if(strcmp(outfile, "") == 0){

    }
    // setting remaining variables for processing
    sum_n_peaks = 0;
    datalen = count_char(infile, '\n');
    cycles = (int) (ceil(datalen / (double) DATA_BUF));
    if(ts == 0)
        ts = TIMESTEP_DEFAULT;

    fp_out = fopen(outfile, "w");
    if(fp_out == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    if(verbose == 1){
        printf("infile: %s\n", infile);
        printf("outfile: %s\n", outfile);
        printf("aux_dir: %s\n", aux_dir);
        printf("datalen: %d\n", datalen);
        printf("cycles: %d\n", cycles);

        printf("\nProgess:  (| = 10 cycles)\n");

    }
    //for(i=0; i<cycles; i++){
    for(int i=0; i<1; i++){

        ind = i * (int)DATA_BUF;
        if(i == cycles-1)
            n = datalen - i * (int)DATA_BUF;
        else
            n = (int)DATA_BUF;
        fetch_data(infile, data, n, ind);
        // main ampd routine, contains malloc
        n_peaks = ampd(data, n, lms, rlms, gamma, sigma, peaks);
        sum_n_peaks += n_peaks;
        concat_peaks(sum_peaks, peaks, ind);
        // save aux
        if(output_all == 1){
            save_data(data, n, detrend_path); // save detrended data
            save_mtx(lms, lms_path);
        }

        free(data);
        free(sigma);
        free(gamma);
        free(peaks);
        free(lms->data);
        free(lms);
        //TODO rlms not needed at all probably
        free(rlms->data);
        free(rlms);
    }

    fclose(fp_out);
    fprintf(stdout, "%d\n", sum_n_peaks);

    return 0;
}
