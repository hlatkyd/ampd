/*
 * ampd.c
 *
 */

#include "ampd.h"

#define TESTING 1

static struct option long_options[] = 
{
    {"infile",required_argument, NULL, 'f'},
    {"help",optional_argument, NULL, 'h'},
    {"outfile",optional_argument, NULL, 'o'},
    {"auxdir",optional_argument, NULL, 'x'},
    {"timestep",optional_argument, NULL, 't'},
    {"overlap", optional_argument, NULL, 'p'},
    {"verbose", optional_argument, NULL, 'v'},
    {"output-all", optional_argument, NULL, 'a'},
    {"output-lms", optional_argument, NULL, 'm'},
    {"length", optional_argument, NULL, 'l'},
    {"rate",optional_argument,NULL, 'r'},
    {NULL, 0, NULL, 0}
};

/**
 * Print description and general usage
 */
void printf_help(){

    printf(
    "AMPD\n"
    "=========================================================================\n"
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
    "directory contains various intermediate data for error checking. "
    "The final peak count is sent to stdout as well."
    "\n"
    "For long data files the processing is done in batches to avoid too high "
    "resource usage. This improves accuracy as well, since data homogeneity is "
    "better preserved in smaller batches. These can overlap redundantly. "
    ""

    "\n\n"

    "Usage from linux command line:\n"
    " $ ampd -f [input file]\n"
    "Optional arguments:\n"
    "\t-o --outfile:\tpath to main output file\n"
    "\t-r --rate:\toutput peak-per-min\n"
    "\t-v --verbose:\tverbose\n"
    "\t-h --help:\tprint help\n"
    "\t-a --output-all:\toutput aux data, local maxima scalogram not included\n"
    "\t-m --output-lms:\toutput local maxima scalogram (high disk space usage)\n"
    "\t-x --auxdir:\taux data root dir, default is cwd\n"
    "\t-t --samplig-rate:\tsampling rate input data in Hz\n"
    "\t-p --overlap:\tmake batches overlapping in time domain\n"
    "\t-l --length:\tdata window length in seconds\n"
    "\n"
            );
}
/**
 * Main handles the command line input parsing and output file management.
 */
int main(int argc, char **argv){

    int opt;
    int verbose = VERBOSE;
    int output_all = OUTPUT_ALL;
    int output_lms = OUTPUT_LMS;
    int output_rate = OUTPUT_RATE;
    double overlap = OVERLAP_DEF;   //overlap ratio
    int n_ovlap;             // overlapping data points

    // timing
    clock_t begin, end;
    double time_spent;

    char infile[MAX_PATH_LEN] = {0};
    char infile_basename[MAX_PATH_LEN]; // input, without dir and extension
    char outfile[MAX_PATH_LEN] = {0}; // main output file with indices of peaks
    char cwd[MAX_PATH_LEN]; // current directory
    FILE *fp_out;
    // aux output files
    char raw_path[MAX_PATH_LEN];
    char detrend_path[MAX_PATH_LEN];
    char lms_path[MAX_PATH_LEN];
    char rlms_path[MAX_PATH_LEN];
    char gamma_path[MAX_PATH_LEN];
    char sigma_path[MAX_PATH_LEN];
    char peaks_path[MAX_PATH_LEN];
    char preproc_path[MAX_PATH_LEN];

    // batch processing
    float *data;        //  timeseries
    int n;              // number of elements in timeseries, in a batch, dynamic
    int ind;            // index of next batch
    int data_buf = DATA_BUF_DEF;
    int datalen;        // full data length
    int cycles;         // number of data batches
    int n_peaks;        // main output, number of peaks
    int sum_n_peaks;    // summed peak number from all batches
    int *sum_peaks;      // concatenated vector from peak indices
    double data_buf_sec = 0.0;
    double a, b; // linear fitting
    // ampd routine pointers
    struct ampd_param *param;
    double sampling_rate;
    int l;
    struct fmtx *lms;
    double *gamma;
    double *sigma;
    int *peaks;

    // main output file, containing only the peak indices
    char outfile_def[] = "ampd.out.peaks"; //

    // aux output paths
    char aux_dir[MAX_PATH_LEN] = {0};
    char batch_dir[MAX_PATH_LEN] = {0};
    char aux_dir_def[] = "ampd_out"; // full default is cwd plus this
    // parse options
    while((opt = getopt_long(argc,argv,"hmvf:o:ax:p:l:",long_options,NULL)) != -1){
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
            case 'm':
                output_lms = 1;
                break;
            case 't':
                sampling_rate = atoi(optarg);
                break;
            case 'x':
                strcpy(aux_dir, optarg);
                break;
            case 'p':
                overlap = atof(optarg);
                if(overlap > 1){
                    printf("error: overlap higher than 1.0 not permitted\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'l':
                data_buf_sec = atof(optarg);
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
    begin = clock();
    getcwd(cwd, sizeof(cwd));
    extract_raw_filename(infile, infile_basename, sizeof(infile_basename));
    // setting to defaults if arguments were not given
    if(strcmp(outfile,"")==0){
        snprintf(outfile, sizeof(outfile), "%s/%s",cwd, outfile_def);
    }
    if(strcmp(aux_dir,"")==0){
        snprintf(aux_dir, sizeof(aux_dir), "%s/%s",cwd, aux_dir_def);
    }
    if(data_buf_sec != 0.0){
        printf("here!\n");
        data_buf = (int)(data_buf_sec * sampling_rate);
    }
    // setting remaining variables for processing
    sum_n_peaks = 0;
    datalen = count_char(infile, '\n');
    n_ovlap = (int)((double)data_buf * overlap);
    cycles = (int) (ceil(datalen / (double) (data_buf - n_ovlap)));
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
        printf("sampling_rate: %.5lf\n",sampling_rate);
        printf("data_buf: %d\n",data_buf);
        printf("cycles: %d\n", cycles);
        printf("output-all: %d\n",output_all);

    }
    /*
     * Processing
     *
     */
    for(int i=0; i<cycles; i++){
        if(TESTING == 1){
            if(i>0)
                break;
        }

        // settings aux output paths
        snprintf(batch_dir, sizeof(batch_dir),"%s/batch_%d",aux_dir,i);
        snprintf(detrend_path,sizeof(detrend_path),"%s/detrend.dat",batch_dir);
        snprintf(raw_path,sizeof(raw_path),"%s/raw.dat",batch_dir);
        snprintf(lms_path,sizeof(lms_path),"%s/lms.dat",batch_dir);
        snprintf(rlms_path,sizeof(lms_path),"%s/rlms.dat",batch_dir);
        snprintf(gamma_path, sizeof(gamma_path),"%s/gamma.dat",batch_dir);
        snprintf(sigma_path, sizeof(sigma_path),"%s/sigma.dat",batch_dir);
        snprintf(peaks_path, sizeof(peaks_path),"%s/peaks.dat",batch_dir);
        snprintf(preproc_path, sizeof(preproc_path),"%s/smoothed.dat",batch_dir);

        ind = i * (int)(data_buf - n_ovlap);
        if(i == cycles-1){
            break; //TODO del test
            //n_init = datalen - i * (int)data_buf;
            n = datalen - i * (int)data_buf;
        }
        else {
            //n_init = (int)data_buf;
            n = (int)data_buf;
        }
        // init matrix and other array pointers
        l = (int)ceil(n/2)-1;
        lms = malloc_fmtx(l, n);
        gamma = malloc(sizeof(double)*l);
        sigma = malloc(sizeof(double)*n);
        peaks = malloc(sizeof(int)*n);

        // getting data
        data = malloc(sizeof(float)*n);
        fetch_data(infile, data, n, ind);
        if(output_all == 1)
            save_data(data, n, raw_path,"float"); // save raw data
        movingavg(data, n, 5);
        if(output_all == 1)
            save_data(data, n, preproc_path,"float"); // save smoothed data

        // main ampd routine, contains malloc
        n_peaks = ampdcpu(data, n, param, lms, gamma, sigma, peaks);
        //catch_false_pks(peaks, &n_peaks, ts, thresh);
        sum_n_peaks += (int)(n_peaks * (1.0-overlap));
        if(verbose == 1){
            printf("batch=%d, n_peaks=%d, sum_n_peaks=%d\n",i,n_peaks,sum_n_peaks);
        }

        //concat_peaks(sum_peaks, peaks, ind);
        // save aux
        if(output_all == 1){
            //save_data(data, n, raw_path); // save raw data
            save_data(data, n, detrend_path,"float"); // save detrended data
            save_data(sigma, n, sigma_path, "double");
            save_data(gamma, l, gamma_path, "double");
            save_data(peaks, n_peaks, peaks_path, "int"); // save peak indices
            if(output_lms == 1){
                save_fmtx(lms, lms_path);
            }
        }

        free(data);
        free(sigma);
        free(gamma);
        free(peaks);
        for(int j=0; j<l; j++)
            free(lms->data[j]);
        free(lms->data);
        free(lms);
        free(param);
    }
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC; 
    if(verbose == 1)
        printf("runtime = %lf\n",time_spent);
    fprintf(fp_out, "%d\n", sum_n_peaks);
    fprintf(stdout, "%d\n", sum_n_peaks);
    fclose(fp_out);
    return 0;
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
 * Function: mkpath
 * ----------------
 *  Recursively create directories. Return -1 on error, 0 on success.
 *  Use mode 0755 for usual read-exec access or 0777 for read-write-exec.
 */
int mkpath(char *file_path, mode_t mode){

    assert(file_path && *file_path);
    for (char* p = strchr(file_path + 1, '/'); p; p = strchr(p + 1, '/')) {
        *p = '\0';
        if (mkdir(file_path, mode) == -1) {
            if (errno != EEXIST) {
                *p = '/';
                return -1;
            }
        }
        *p = '/';
    }
    return 0;
}

/**
 * Save list into file, one value per line.
 * Creates necessary directories.
 * Type represents the type of data: ["float", "double", "int"]
 */
void save_data(void *indata, int n, char *path, char *type){

    FILE *fp;
    float *fdata; int *idata; double *ddata;

    if(mkpath(path, 0777) == -1){
        fprintf(stderr, "cannot make path %s\n",path);
        exit(EXIT_FAILURE);
    }
    fp = fopen(path, "w");
    if(fp == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    // conditional fprint
    if(strcmp(type, "float")){
        fdata = (float*)&indata;
        for(int i=0; i<n; i++){
            fprintf(fp, "%.3f\n",fdata[i]);
        }
    }
    if(strcmp(type, "double")){
        idata = (int*)&indata;
        for(int i=0; i<n; i++){
            fprintf(fp, "%.5lf\n",ddata[i]);
        }
    }

    if(strcmp(type, "int")){
        ddata = (double*)&indata;
        for(int i=0; i<n; i++){
            fprintf(fp, "%d\n",idata[i]);
        }
    }
    fclose(fp);
    return;
}

/**
 * Save matrix to a tab delimited file, without any headers.
 * Creates necessary directories.
 */
void save_fmtx(struct fmtx *mtx, char *path){

    FILE *fp;
    int i, j;
    if(mkpath(path, 0777) == -1){
        fprintf(stderr, "cannot make path %s\n",path);
        exit(EXIT_FAILURE);
    }

    fp = fopen(path, "w");
    if(fp == NULL){
        perror("fopen, exiting...");
        exit(EXIT_FAILURE);
    }
    //printf("rows, cols: %d, %d\n", mtx->rows, mtx->cols);
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
int fetch_data(char *path, float *data, int n, int ind ){

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
    fclose(fp);
    return 0;
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

void set_ampd_param(struct ampd_param *p){

}



