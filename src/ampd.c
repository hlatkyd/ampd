/*
 * ampd.c
 *
 */

#include "ampd.h"

#define TESTING 0
#define CONF_PATH "ampd.conf"

/* only for getopt */
#define ARG_OUTPUT_ALL 1
#define ARG_OUTPUT_LMS 2
#define ARG_OUTPUT_RATE 3
#define ARG_OVERLAP 4
#define ARG_PREPROCESS 5
#define ARG_HPFILT 6
#define ARG_LPFILT 7

static int verbose;
static int output_all;
static int output_lms;
static int output_rate;
static int preprocess;
static int smooth_flag; // TODO
static struct option long_options[] = 
{
    {"infile",required_argument, NULL, 'f'},
    {"outfile",required_argument, NULL, 'o'},
    {"auxdir",required_argument, NULL, 'a'},
    {"verbose", optional_argument, NULL, 'v'},
    {"help",no_argument, NULL, 'h'},
    {"datatype",optional_argument,NULL, 't'},
    {"sampling-rate",required_argument, NULL, 'r'},
    {"batch-length", required_argument, NULL, 'l'},
    {"overlap", required_argument, NULL, ARG_OVERLAP},
    {"hpfilt", optional_argument, NULL, ARG_HPFILT},
    {"lpfilt", optional_argument, NULL, ARG_LPFILT},
    {"smooth", no_argument, &smooth_flag, 1},
    {"output-all", no_argument, NULL, ARG_OUTPUT_ALL},
    {"output-lms", no_argument, NULL, ARG_OUTPUT_LMS},
    {"output-rate", no_argument, NULL, ARG_OUTPUT_RATE},
    {"preprocess", no_argument, NULL, ARG_PREPROCESS},
    {NULL, 0, NULL, 0}
};
static char optstring[] = "hv:f:o:a:l:r:t:";

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
    // TODO cleanup
    "Usage from linux command line:\n\n"
    " $ ampd -f [input file]\n\n"
    "Optional arguments:\n"
    "\t-o --outfile:\tpath to main output file\n"
    "\t-a --auxdir:\taux data root dir, default is cwd\n"
    "\t-v --verbose:\tverbose\n"
    "\t-h --help:\tprint help\n"
    "\t-r --samplig-rate:\tsampling rate input data in Hz\n"
    "\t-l --batch-length:\tdata window length in seconds\n"
    "\t--lpfilt:\tapply lowpass filter\n"
    "\t--hpfilt:\tapply highpass filter\n"
    "\t--overlap:\tmake batches overlapping in time domain\n"
    "\t--output-all:\toutput aux data, local maxima scalogram not included\n"
    "\t--output-lms:\toutput local maxima scalogram (high disk space usage)\n"
    "\t--output-rate:\toutput peak-per-min\n"
    "\n"
            );
}
/**
 * Main handles the command line input parsing and output file management.
 */
int main(int argc, char **argv){

    int opt;
    struct ampd_config *conf;
    char datatype[32] = {0};
    // timing
    clock_t begin, end;
    double time_spent;

    char infile[MAX_PATH_LEN] = {0};
    char infile_basename[MAX_PATH_LEN]; // input, without dir and extension
    char outfile[MAX_PATH_LEN] = {0}; // main output file with indices of peaks
    FILE *fp_out;        // main output file containing the peak indices
    char cwd[MAX_PATH_LEN]; // current directory
    /*
     * aux output files
     */
    char raw_path[MAX_PATH_LEN];
    char detrend_path[MAX_PATH_LEN];
    char lms_path[MAX_PATH_LEN];
    char rlms_path[MAX_PATH_LEN];
    char gamma_path[MAX_PATH_LEN];
    char sigma_path[MAX_PATH_LEN];
    char peaks_path[MAX_PATH_LEN];
    char preproc_path[MAX_PATH_LEN];
    char param_path[MAX_PATH_LEN];
    char batch_param_path[MAX_PATH_LEN];

    /*
     * batch processing
     */
    float *data;        //  timeseries
    int n;              // number of elements in timeseries, in a batch, dynamic
    int ind;            // index of next batch
    int data_buf;
    int datalen;        // full data length
    double batch_length;
    double overlap;          //overlap ratio
    int n_overlap;             // overlapping data points
    int cycles;         // number of data batches
    int sum_n_peaks;    // summed peak number from all batches
    int *sum_peaks;      // concatenated vector from peak indices
    double peaks_per_min;
    double data_buf_sec = 0.0;
    struct batch_param *bparam; // only for outputting batch utility parameters

    /* 
     * filtering
     */
    double cutoff_freq;
    double highpassfilt = 0.0;
    double lowpassfilt = 0.0;
    /*
     * ampd routine
     */
    struct ampd_param *param;
    double sampling_rate = 0.0;
    int l;
    struct fmtx *lms;
    double *gamma;
    double *sigma;
    int *peaks;
    int n_peaks;        // main output, number of peaks

    // main output file, containing only the peak indices
    char outfile_def[] = "ampd.out.peaks"; //

    // aux output paths
    char aux_dir[MAX_PATH_LEN] = {0};
    char batch_dir[MAX_PATH_LEN] = {0};
    char aux_dir_def[] = "ampd_out"; // full default is cwd plus this
    if(argc == 1){
        printf_help();
        return 0;
    }
    // load config
    char conf_path[] = CONF_PATH;
    conf = malloc(sizeof(struct ampd_config));
    load_config(conf_path, conf, NULL);
    // parse options
    while((opt = getopt_long(argc,argv,optstring,long_options,NULL)) != -1){
        switch(opt){
            case 'h':
                printf_help();
                return 0;
            case 'v':
                if(optarg != NULL){
                    verbose = atoi(optarg);
                }
                else
                    verbose = 1;
                break;
            case 't':
                if(optarg == NULL)
                    strcpy(datatype,"def");
                else
                    strcpy(datatype,optarg);
                break;
        }
    }
    optind = 0;
    while((opt = getopt_long(argc,argv,optstring,long_options,NULL)) != -1){
        switch(opt){
            case 'o':
                strcpy(outfile, optarg);
                break;
            case 'f':
                strcpy(infile, optarg);
                break;
            case 'r':
                sampling_rate = atof(optarg);
                break;
            case 'a':
                strcpy(aux_dir, optarg);
                break;
            case 'l':
                batch_length = atof(optarg);
                break;
            case ARG_OUTPUT_ALL:
                output_all = 1;
                break;
            case ARG_OUTPUT_LMS:
                output_lms = 1;
                break;
            case ARG_OVERLAP:
                overlap = atof(optarg);
                if(overlap > 1){
                    printf("error: overlap higher than 1.0 not permitted\n");
                    exit(EXIT_FAILURE);
                }
                //TODO remove once done
                fprintf(stderr,"Warning: overlap > 0 is not finished yet.\n");
                fprintf(stderr,"Defaulting to  overlap = 0.\n");
                overlap = 0;
                break;
            case ARG_PREPROCESS:
                // do filtering and smoothing
                preprocess = 1;
                break;
            case ARG_HPFILT:
                if(optarg == NULL)
                    highpassfilt = conf->hpfilt;
                else 
                    highpassfilt = atof(optarg);
                break;
            case ARG_LPFILT:
                if(optarg == NULL)
                    lowpassfilt = conf->lpfilt;
                else 
                    lowpassfilt = atof(optarg);
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
    param = malloc(sizeof(struct ampd_param));
    bparam = malloc(sizeof(struct batch_param));
    memset(bparam, 0, sizeof(bparam));
    set_ampd_param(param, "resp");
    if(sampling_rate != 0)
        param->sampling_rate = sampling_rate;
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
    data_buf = (int)(batch_length * param->sampling_rate);
    n_overlap = (int)((double)data_buf * overlap);
    cycles = (int) (ceil(datalen / (double) (data_buf - n_overlap)));
    fp_out = fopen(outfile, "w");
    if(fp_out == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    if(verbose > 0){
        printf("infile: %s\n", infile);
        printf("outfile: %s\n", outfile);
        printf("lowpassfilt=%lf\n",lowpassfilt);
        printf("highpassfilt=%lf\n",highpassfilt);
        printf("aux_dir: %s\n", aux_dir);
        printf("sampling_rate: %.5lf\n",param->sampling_rate);
        printf("batch_length: %lf\n",batch_length);
        printf("overlap: %lf, n_overlap: %d\n",overlap, n_overlap);
        printf("datalen: %d\n", datalen);
        printf("data_buf: %d\n",data_buf);
        printf("cycles: %d\n", cycles);
        printf("output-all: %d\n",output_all);
        printf("output-lms: %d\n",output_lms);

    }
    /* fill batch param */
    bparam->cycles = cycles;
    bparam->batch_length = batch_length;
    bparam->sampling_rate = sampling_rate;
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
        snprintf(param_path, sizeof(param_path),"%s/param.txt",batch_dir);
        snprintf(batch_param_path, sizeof(batch_param_path),"%s/bparam.txt",batch_dir);

        ind = i * (int)(data_buf - n_overlap);
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

        // fill some more batch param
        bparam->ind = i;
        bparam->n = n;
        bparam->l = l;
        // getting data
        data = malloc(sizeof(float)*n);
        if(verbose > 1){
            printf("\nfetchig data:\n");
            printf("infile=%s\n",infile);
            printf("data=%p\n",data);
            printf("n=%d, ind=%d\n",n,ind);
        }
        fetch_data(infile, data, n, ind);
        if(output_all == 1)
            save_data(data, n, raw_path,"float"); // save raw data

        /* Filter data depending on datatype. Pulsoxy signals benefit from high
         * pass filtering near the respiration frequency.
         *
         */
        if(lowpassfilt != 0.0){
            tdlpfilt(data, n, sampling_rate, 0.5);
        }
        if(highpassfilt != 0.0){
            tdhpfilt(data, n, sampling_rate, 2);
        }
        if(smooth_flag == 1){
            movingavg(data, n, 5);
        }
        if(output_all == 1)
            save_data(data, n, preproc_path,"float"); // save smoothed data

        // main ampd routine, contains malloc
        n_peaks = ampdcpu(data, n, param, lms, gamma, sigma, peaks);

        //catch_false_pks(peaks, &n_peaks, ts, thresh);
        sum_n_peaks += (int)(n_peaks * (1.0-overlap));
        if(verbose > 0){
            printf("batch=%d/%d, n_peaks=%d, sum_n_peaks=%d\n",
                    i,cycles, n_peaks,sum_n_peaks);
        }
        // calc peak rate
        bparam->n_peaks = n_peaks;
        peaks_per_min = (double)n_peaks / batch_length * 60.0;
        bparam->peaks_per_min = peaks_per_min; 

        //concat_peaks(sum_peaks, peaks, ind);
        //TODO
        // save aux
        if(output_all == 1){
            //save_data(data, n, raw_path); // save raw data
            save_data(data, n, detrend_path,"float"); // save detrended data
            save_data(sigma, n, sigma_path, "double");
            save_data(gamma, l, gamma_path, "double");
            //TODO peak indices not correct if overlap!=0 until concatenation is
            // not addressed
            save_data(peaks, n_peaks, peaks_path, "int"); // save peak indices
            save_ampd_param(param, param_path);
            save_batch_param(bparam, batch_param_path);
        }
        if(output_lms == 1){
            save_fmtx(lms, lms_path);
        }

        free(data);
        free(sigma);
        free(gamma);
        free(peaks);
        for(int j=0; j<l; j++)
            free(lms->data[j]);
        free(lms->data);
        free(lms);
    }
    free(param);
    free(bparam);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC; 
    if(verbose > 0)
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
    fp = fopen(path, "w+");
    if(fp == NULL){
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    // conditional fprint
    if(strcmp(type, "float")==0){
        fdata = (float*)indata;
        for(int i=0; i<n; i++){
            fprintf(fp, "%.3f\n",fdata[i]);
        }
    }
    if(strcmp(type, "double")==0){
        ddata = (double*)indata;
        for(int i=0; i<n; i++){
            fprintf(fp, "%.5lf\n",ddata[i]);
        }
    }
    if(strcmp(type, "int")==0){
        idata = (int*)indata;
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
 * Save ampd_param struct contents into file
 */
void save_ampd_param(struct ampd_param *p, char *path){

    FILE *fp;
    fp = fopen(path, "w+");
    if(fp == NULL){
        fprintf(stderr, "cannot open file for writing: %s\n",path);
        exit(1);
    }
    fprintf(fp, "sampling_rate=%lf\n",p->sampling_rate);
    fprintf(fp, "datatype=%s\n", p->datatype);
    fprintf(fp, "a=%lf\n",p->a);
    fprintf(fp, "rnd_factor=%lf\n",p->rnd_factor);
    fprintf(fp, "fit_a=%lf\n",p->fit_a);
    fprintf(fp, "fit_b=%lf\n",p->fit_b);
    fprintf(fp, "fit_r=%lf\n",p->fit_r);
    fprintf(fp, "lambda=%d\n",p->lambda);
    fprintf(fp, "sigma_thresh=%lf\n",p->sigma_thresh);
    fprintf(fp, "peak_thresh=%lf\n",p->peak_thresh);

    fclose(fp);

}
/**
 * Save utility parameters to batch dir
 */
void save_batch_param(struct batch_param *p, char *path){

    FILE *fp;
    fp = fopen(path, "w+");
    if(fp == NULL){
        fprintf(stderr, "can't open file for writing %s\n",path);
        exit(1);
    }
    fprintf(fp,"ind=%d\n",p->ind);
    fprintf(fp,"cycles=%d\n",p->cycles);
    fprintf(fp,"n=%d\n",p->n);
    fprintf(fp,"l=%d\n",p->l);
    fprintf(fp,"batch_length=%lf\n",p->batch_length);
    fprintf(fp,"sampling_rate=%lf\n",p->sampling_rate);
    fprintf(fp,"n_peaks=%d\n",p->n_peaks);
    fprintf(fp,"peaks_per_min=%lf\n",p->peaks_per_min);
    fclose(fp);
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

void set_ampd_param(struct ampd_param *p, char *type){

    strcpy(p->datatype,type);
    p->a = DEF_A;
    p->rnd_factor = DEF_RND_FACTOR;
    p->sampling_rate = DEF_SAMPLING_RATE;
    if(strcmp(type, "resp")==0){
        // respiration optimized
        p->sigma_thresh = RESP_SIGMA_THRESHOLD;
        p->peak_thresh = RESP_PEAK_THRESHOLD;

    }
    else if(strcmp(type, "puls")==0){
        // pulsoxy optimized
        p->sigma_thresh = PULS_SIGMA_THRESHOLD;
        p->peak_thresh = PULS_PEAK_THRESHOLD;
    }
    else {
        // default
        p->sigma_thresh = DEF_SIGMA_THRESHOLD;
        p->peak_thresh = DEF_PEAK_THRESHOLD;
    }

}
/**
 * Load default config from file
 */
void preload_config(char *path, struct ampd_config *conf){
}
/**
 * Load datatype specific custom settings from config file
 */
void load_config(char *path, struct ampd_config *conf, char *type){

}
